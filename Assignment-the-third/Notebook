# Demultiplexing and Index Swapping – Assignment the First

### Part 1 


In Assignment the First, we will develop a strategy to de-multiplex samples to create 48 FASTQ files that contain acceptable index pairs (read1 and read2 for 24 different index pairs), two FASTQ files with index-hopped reads-pairs, and two FASTQ files undetermined (non-matching or low quality) index-pairs.

De-multiplexing is necessary for downstream analyses.

We submitted 24 indexed (dual matched) libraries for sequencing. The indexes are:

```B1	GTAGCGTA    A5	CGATCGAT    C1	GATCAAGG
B9	AACAGCGA    C9	TAGCCATG    C3	CGGTAATC
B3	CTCTGGAT    C4	TACCGGAT    A11	CTAGCTCA
C7	CACTTCAC    B2	GCTACTCT    A1	ACGATCAG
B7	TATGGCAC    A3	TGTTCCGT    B4	GTCCTAAG
A12	TCGACAAG    C10	TCTTCGAC    A2	ATCATGCG
C2	ATCGTGGT    A10	TCGAGAGT    B8	TCGGATTC
A7	GATCTTGC    B10	AGAGTCCA    A8	AGGATAGC
```

**Dual matched is when the index matches on both sides** 

We have 4 FASTQ input files on Talas: ```/projects/bgmp/shared/2017_sequencing/```

```
1294_S1_L008_R1_001.fastq.gz
1294_S1_L008_R2_001.fastq.gz
1294_S1_L008_R3_001.fastq.gz
1294_S1_L008_R4_001.fastq.gz
```
Initial data exploration includes the following: 
**Use zcat before commands to view the contents of the files without unzipping them**
- How big is the each file? ```ls -lah```
```
-rw-r-xr--+  1 coonrod  is.racs.pirg.bgmp  20G Jul 30  2018 1294_S1_L008_R1_001.fastq.gz
-rw-r-xr--+  1 coonrod  is.racs.pirg.bgmp 2.6G Jul 30  2018 1294_S1_L008_R2_001.fastq.gz
-rw-r-xr--+  1 coonrod  is.racs.pirg.bgmp 2.8G Jul 30  2018 1294_S1_L008_R3_001.fastq.gz
-rw-r-xr--+  1 coonrod  is.racs.pirg.bgmp  21G Jul 30  2018 1294_S1_L008_R4_001.fastq.gz 
```
- Number of lines? ```zcat <filename> | wc -l``` 
```1452986940 1294_S1_L008_R1_001.fastq.gz 
1452986940 1294_S1_L008_R2_001.fastq.gz
1452986940 1294_S1_L008_R3_001.fastq.gz
1452986940 1294_S1_L008_R4_001.fastq.gz
```
Index files will have sequences much shorter than files containing reads so just by looking at the files we can tell which one is which. Also the order should always be this 
R1 --> fw read
R3 --> Index 1
R3 --> Index 2 (has to be reverse complemented here but this is not always the case)
R4 --> rv read 

To find the Phred encoding, we look at the quality score lines, and see where on the ASCII table do the values fall. Here they are all Phred+33

To find the read length for each FASTQ file: 
```
zcat <filename> | grep -v "@"|head -n 1| wc -c 
Then the answer -1 to account for tab character 
```

To find the number of sequences in the index files that contain an undetermined base "N" I used the following bash command:

```zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz | awk 'NR%4==2' | sort| uniq -c | grep -c "N"```

The awk command looks a bit weird but it makes sense and I should use it this way easier than using greps. NR refres to the number of records followed by the normal modulus command. So it is going through every sequence line in the record. 

### Part 2
The second part included writing a pseudocode to understand the problem and envision how we will be solving it, step by step. It was helpful to look at peers pseudocode and it made me make slight changes to my own strategy. 


### Part 3
 Writing the actual code took some effort. I struggled with organizing my thoughts and going step by step. It was helpful to use argsparse, and split my code into functions that do each step. The main error i made here is opening the files and writing to the files within the loops, each time i was going to write to the file i need to open it, which would have taken forever. So, I wrote several functions; a function open all the files before the while loop, a function to write to the files using the join command since all my records are saved in lists, and a function to close all the files at the end. 

I also realized that I did not have a clear understanding of file handles at the start. I made a dictionary that has the file names and respective file handles, so when the files are opened, I can just write to each file according to its barcode. 

Another mistake was not accounting for valid low-quality barcodes. So my total reads at the beginning was less than what it should be. It was an easy fix; I added another else: block after counting the hopped reads, which counts what is left out, added that to my total, but these still ended up in the unknown files since we don't have to write them to their own file. 

finally I ran the following sbatch command:

```
#!/bin/bash
#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH -c 1
#SBATCH --time=0-3

/usr/bin/time -v python ./demultiplex.py -r1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -r2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -r3 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz -r4 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -q 30 -o ./output-FASTQ-q30_updated
```

slurm-8560412.out 

```
Done reading files
	Command being timed: "python ./demultiplex.py -r1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -r2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -r3 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz -r4 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -q 30 -o ./output-FASTQ-q30_updated"
	User time (seconds): 3070.31
	System time (seconds): 34.02
	Percent of CPU this job got: 98%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 52:16.32
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 248640
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 38104
	Voluntary context switches: 44466
	Involuntary context switches: 10978
	Swaps: 0
	File system inputs: 0
	File system outputs: 0
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
```


