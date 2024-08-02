# Assignment the First

## Part 1
1. Be sure to upload your Python script. Provide a link to it here:

| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz | Read1 |  101| Phred+33  |
| 1294_S1_L008_R2_001.fastq.gz | Index1 | 8 |  Phred+33 |
| 1294_S1_L008_R3_001.fastq.gz | Index2 | 101 |  Phred+33 |
| 1294_S1_L008_R4_001.fastq.gz | Read2 | 8 |  Phred+33 |

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.

    [R1 Distribution](https://github.com/Mahmoud-56/Demultiplex/blob/master/Assignment-the-first/quality_score_distribution_Read1.png)
     [R2 Distribution](https://github.com/Mahmoud-56/Demultiplex/blob/master/Assignment-the-first/quality_score_distribution_Index1.png)  
     [R3 Distribution](https://github.com/Mahmoud-56/Demultiplex/blob/master/Assignment-the-first/quality_score_distribution_Index2.png) 
     [R4 Distribution](https://github.com/Mahmoud-56/Demultiplex/blob/master/Assignment-the-first/quality_score_distribution_Read2.png)       
  

    Quality score cutoff:

    I will set the cutoff at 30 since by looking at the average reads graph, the average is slightly above 30. Also according to Illumina, Q30 quality score is considered the benchmark for high quality next-generation sequencing.

    https://www.illumina.com/documents/products/technotes/technote_Q-Scores.pdf


## Part 2
1. Define the problem
   
   We need to categorize each read pair from the forward and reverse files as a match, mismatch, or unknown. With 24 different indexes, we should end up with 48 match files (24 for each direction), 2 mismatch files, and 2 unknown files. The second index (R3) needs to be reverse complemented before the comparison with R2. Each header in the "matched" and "unmatched" files should include their respective index pair appended to the header. 
2. Describe output 

    The script should input the following:
- 24 fastq files for matching reads for the fw reads
- 24 fastq files for unmatched (hopped) reads for the rv reads
- 2 files for hopped reads (one for fw and one for rv)
- 2 files for unknown reads (one for fw and one for rv)

    The script will also output some statistics about our resutls including: 

- The number of matching index pairs 
- The number of hopped index pairs (percentage of hopped reads)
- The number of unknown index pairs 
3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).

    The 4 input files are here [4 input FASTQ files](../TEST-input_FASTQ)
   
    The 6 output files are here [>=6 expected output FASTQ files](../TEST-output_FASTQ)


4. Pseudocode
    
    [Pseudocode](https://github.com/Mahmoud-56/Demultiplex/blob/master/pseudocode) 
5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement




```def reverse_complement(sequence: str) -> str:
    '''Returns the reverse complement of a DNA sequence.'''
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    rev_comp = ''
    for base in sequence:
        rev_comp = complement[base] + rev_comp
    return rev_comp
    Input: 'ACC'
    
    Expected output: 'TGG' 
    assert reverse_complement('ATCG') == 'CGAT'
````
```
def create_new_header(r1_header: str, r4_header: str, r2_seq: str, r3_seq: str) -> (str, str):
   ''' Takes in the header of R1 and R4 and appends the index pair for that sequence to the end of the header in the fw R1 and the rv R4 files.'''


    index_pair = f"{r2_seq}-{r3_seq}"
    r1_new_header = f"{r1_header} {index_pair}"
    r4_new_header = f"{r4_header} {index_pair}"


    Returns r1_new_header and r4_new_header.
````
````
def read_fastq_files(r1_file, r2_file, r3_file, r4_file):

''Opens four zipped fastq files and reads one record at a time simultaneously from each file. Then it defines each line from each record with it's own variable '''

    return processed_records
````

```
def write_output_files(processed_records, output_file):
'''writes the processed records to their respective output file'''
````
