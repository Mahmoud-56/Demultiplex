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
    2. **YOUR ANSWER HERE**
    3. **YOUR ANSWER HERE**
    
## Part 2
1. Define the problem
   We need to categorize each read pair from the forward and reverse files as a match, mismatch, or unknown. With 24 different indexes, we should end up with 48 match files (24 for each direction), 2 mismatch files, and 2 unknown files. The second index (R3) needs to be reverse complemented before comparison.
2. Describe output 
3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).
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
````

    