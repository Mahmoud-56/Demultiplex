#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH --time=0-3




conda activate bgmp-velvet


/usr/bin/time -v zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz | python ./multiplex.py -s 101 -n Read1

/usr/bin/time -v zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz | python ./multiplex.py -s 8 -n Index1

/usr/bin/time -v zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz | python ./multiplex.py -s 8 -n Index2

/usr/bin/time -v zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz | python ./multiplex.py -s 101 -n Read2


#The stdout of the zcat command is the stdin for the python script

