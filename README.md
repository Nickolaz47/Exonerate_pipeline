# Exonerate Pipeline

>The main objective of this pipeline is speed up the Exonerate analysis using core parallelism. Furthermore, the script
>classifies the alignments in not aligned, not perfect and perfect.

## Dependencies

- External programs
  - [Exonerate](https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate)
- Main Python libraries
  - biopython
  - numpy
  - the others are present in the requirements file

## Installation

After install the external program, download the repository:
```
git clone https://github.com/Nickolaz47/Exonerate_pipeline.git
```
Enter the folder and install the python libraries:
```
pip install -r requirements.txt
```

## Usage

This program has two modes: **mrna**, this mode is responsible for align nucleotides against nucleotides; **ptn**, where
aminoacids are aligned against nucleotides. Both inputs must be in fasta format.

- Usage:
```
python3 exonerate_pipe.py -q gene/protein fasta input -t genome fasta input  -c core numbers -m alignment mode
```

## Parameters

- **-q**: Gene/protein fasta query.
- **-t**: Genome fasta subject. 
- **-c**: Core numbers to use (default 4).
- **-m**: Alignment mode (mrna or ptn).
- **-h**: Show the help.

## Algorithm

The algorithm basically splits the fasta query to speed up the Exonerate analysis and classifies the GFFs. 

### 1. Exonerate

The Exonerate alignment is pretty fast and the command line was built to improve that feature. The main parameters 
responsible for this are: *--exhaustive False* to use alignment heuristics sparing time; *-n 1* to show only the best 
alignment (once in a while can show more than one great alignment); *--percent 90* to report alignments with 90% of the 
maximal score optainable for that query.

Exonerate command line:
```
exonerate -q file -t target --exhaustive False --showtargetgff -n 1 -m est2genome --softmasktarget True -M 1000 --showalignment False --showvulgar False --percent 90 --refine region --ryo "# %qi %ql %qab %qae %qal %ti %tab %tae %pi\\n" > output
```

### 2. Classification

To classify the alignments the parameter *--ryo "# %qi %ql %qab %qae %qal %ti %tab %tae %pi\\n"* is used for give some 
information about the query length, alignment lenght, etc. Through this data is possible to classify the alignment in: 
**perfect** 100% of identity and coverage; **not perfect** less than 100% of identity and coverage; **not aligned** no 
alignment was generated.
