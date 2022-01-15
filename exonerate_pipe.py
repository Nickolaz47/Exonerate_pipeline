import os
import subprocess
from concurrent.futures import ProcessPoolExecutor
from argparse import ArgumentParser
from Bio import SeqIO


def multifasta2single(fasta_inp):
    """
    It gets the multifasta file, breaks into many singlefasta files, and writes
    them into geneset folder.
    :param fasta_inp: Fasta input
    :return:
    """
    # Storing the sequence in a list
    with open(fasta_inp, 'r') as inp:
        seqs = [record for record in SeqIO.parse(inp, 'fasta')]
    # Making the dir
    try:
        os.mkdir('geneset')
    except Exception as error:
        print(f'The error is {error}')
    os.chdir(os.path.join(os.getcwd(), 'geneset'))
    # Writing each sequence
    for seq in seqs:
        file_name = seq.id + '.fasta'
        with open(file_name, 'w') as out:
            SeqIO.write(seq, out, 'fasta')
    os.chdir(os.path.join(os.getcwd(), '..'))
    return


def exonerate(file, target, mode):
    """
    This function runs Exonerate.
    :param file: Fasta file (query)
    :param target: Fasta file to align the query
    :param mode: mRNA or ptn
    :return:
    """
    # Defining output
    output = file.split('.fasta')[0] + '.gff'
    # Defining the command line
    if mode == 'mrna':
        exonerate_line = f'exonerate -q {file} -t {target} ' \
                         f'--exhaustive False --showtargetgff -n 1' \
                         f' -m est2genome --softmasktarget True' \
                         f' -M 1000 --showalignment False' \
                         f' --showvulgar False --percent 90' \
                         f' --refine region --ryo "# %qi %ql %qab %qae' \
                         f' %qal %ti %tab %tae %pi\\n" > {output}'
    else:
        exonerate_line = f'exonerate -q {file} -t {target} ' \
                         f'--exhaustive False --showtargetgff -n 1' \
                         f' -m protein2genome --softmasktarget True' \
                         f' -M 1000 --showalignment False' \
                         f' --showvulgar False --percent 90' \
                         f' --refine region --ryo "# %qi %ql %qab %qae' \
                         f' %qal %ti %tab %tae %pi\\n" > {output}'
    # Running exonerate
    subprocess.run(exonerate_line, shell=True)
    return


def clean_geneset():
    """
    This function moves the results to a new folder and deletes the geneset
    folder.
    :return:
    """
    # Listing the folders
    folder = 'geneset'
    # Creanting the results folder
    os.mkdir('results')
    # Deleting the geneset folders and fasta files, and moving the gffs files
    # to results
    os.chdir(os.path.join(os.getcwd(), folder))
    files = os.listdir()
    for file in files:
        if file.endswith('.fasta') or os.path.islink(file):
            os.remove(file)
        else:
            os.replace(os.path.join(os.getcwd(), file),
                       os.path.join(os.getcwd(), f'../results/{file}'))
    os.chdir(os.path.join(os.getcwd(), '..'))
    os.rmdir(folder)
    return


def organize_results():
    """
    This function organizes the results in perfect, not perfect and not aligned.
    :return:
    """
    os.chdir(os.path.join(os.getcwd(), 'results'))
    # Creating the directories
    gffs = os.listdir()
    os.mkdir('perfect')
    os.mkdir('not_perfect')
    os.mkdir('not_aligned')
    for gff in gffs:
        with open(gff, 'r') as inp:
            inp = inp.readlines()
            inp = inp[-2].strip()
            if len(inp.split()) == 10:
                query_len = inp.split()[2]
                align_len = inp.split()[5]
                identity = inp.split()[9]
                if query_len == align_len and identity == '100.00':
                    os.replace(os.path.join(os.getcwd(), gff),
                               os.path.join(os.getcwd(), f'perfect/{gff}'))
                else:
                    os.replace(os.path.join(os.getcwd(), gff),
                               os.path.join(os.getcwd(), f'not_perfect/{gff}'))
            else:
                os.replace(os.path.join(os.getcwd(), gff),
                           os.path.join(os.getcwd(), f'not_aligned/{gff}'))
    os.chdir(os.path.join(os.getcwd(), '..'))
    return


def count_alignments():
    """
    This function counts the number of perfect, not perfect, and not aligned
    results.
    :return:
    """
    count_dic = {}
    os.chdir(os.path.join(os.getcwd(), 'results'))
    folders = os.listdir()
    # Counting the number of files in each folder
    for folder in folders:
        os.chdir(os.path.join(os.getcwd(), folder))
        count_dic[folder] = len(os.listdir())
        os.chdir(os.path.join(os.getcwd(), '..'))
    os.chdir(os.path.join(os.getcwd(), '..'))
    return count_dic


def pipe(query, target, mode, cores=4):
    """
    This function runs the complete pipeline.
    :param query: Fasta query
    :param target: Fasta target
    :param mode: Alignment mode (mRNA or ptn)
    :param cores: Number of cores to use
    :return:
    """
    # Running the pipeline
    print('Starting preprocessing...')
    multifasta2single(query)
    folder = 'geneset'
    # Creating a symlink for the target
    os.symlink(os.path.join(os.getcwd(), target),
               os.path.join(os.getcwd(), f'{folder}/{target}'))
    print('Running Exonerate...')
    # Entering into the folder
    os.chdir(os.path.join(os.getcwd(), folder))
    # Excluding the target file as query
    files = [file for file in os.listdir() if file != target]
    # Creating generators with the same lenght of the files list to run the
    # map function
    mode_gen = (mode for _ in range(len(files)))
    target_gen = (target for _ in range(len(files)))
    # Generating many Exonerate processes limited by the cores number
    with ProcessPoolExecutor(max_workers=cores) as executor:
        executor.map(exonerate, files, target_gen, mode_gen)
    # Getting out of the current folder
    os.chdir(os.path.join(os.getcwd(), '..'))
    print('Starting posprocessing...')
    # Cleaning the geneset folder
    clean_geneset()
    # Classifying the results
    organize_results()
    print('Finished.\n')
    # Counting the results
    count_dic = count_alignments()
    print('Alignments information:\n')
    for k, v in count_dic.items():
        print(f'{k} -> {v} alignments')
    print(f'total -> {sum(count_dic.values())} alignments')
    return


if __name__ == '__main__':
    # Defining the program
    parser = ArgumentParser(description='It runs the Exonerate pipeline.')
    parser.add_argument('-q', required=True, help='Fasta query')
    parser.add_argument('-t', required=True, help='Fasta subject')
    parser.add_argument('-c', required=True, help='Number of cores to use')
    parser.add_argument('-m', required=True, choices=['mrna', 'ptn'],
                        help='Mode to align the queries to target. mrna -> '
                             'nucleotide x nucleotide; ptn -> aminoacid x '
                             'nucleotide')
    # Parsing the arguments
    args = parser.parse_args()
    # Running the pipeline
    pipe(args.q, args.t, args.m, int(args.c))
