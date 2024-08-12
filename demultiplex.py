#!/usr/bin/env python

import argparse
import bioinfo
import gzip

known_barcodes = {
    "GTAGCGTA", "CGATCGAT", "GATCAAGG", "AACAGCGA", "TAGCCATG", "CGGTAATC", "CTCTGGAT",
    "TACCGGAT", "CTAGCTCA", "CACTTCAC", "GCTACTCT", "ACGATCAG", "TATGGCAC", "TGTTCCGT",
    "GTCCTAAG", "TCGACAAG", "TCTTCGAC", "ATCATGCG", "ATCGTGGT", "TCGAGAGT", "TCGGATTC",
    "GATCTTGC", "AGAGTCCA", "AGGATAGC"
}

number_of_low_quality_reads = 0
number_of_matched_reads = 0
number_of_hopped_reads = 0
number_of_unk_reads = 0
number_of_reads_per_index = {}  # keys are each barcode in the set, value is how many times that high-quality matched read occurred.

def get_args():
    parser = argparse.ArgumentParser(description="Demultiplexing")
    parser.add_argument("-r1", "--r1_file", type=str, required=True)
    parser.add_argument("-r2", "--r2_file", type=str, required=True)
    parser.add_argument("-r3", "--r3_file", type=str, required=True)
    parser.add_argument("-r4", "--r4_file", type=str, required=True)
    parser.add_argument("-q", "--quality_threshold", help="quality score cutoff", type=int)
    parser.add_argument("-o", "--output_dir", type=str, required=True)
    return parser.parse_args()

args = get_args()

def reverse_complement(sequence: str):
    '''Returns the reverse complement of a DNA sequence.'''
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    rev_comp = ''
    for base in sequence:
        rev_comp = complement[base] + rev_comp
    return rev_comp

def write_fastq(output_file_handle, record):
    '''writes to output files, takes in a list and outputs in fastq format using join function'''
    output_file_handle.write("\n".join(record) + "\n")

def open_output_files(output_dir):
    '''Opens output FASTQ files for each barcode plus hopped and unk. 'a' refers to append mode
       Parameters:
            output_dir (str): The directory where the output files will be created.

        Returns:
            A dictionary where keys are filenames and values are file handles for writing.
       '''
    file_handles = {}
    for barcode in known_barcodes:
        file_handles[f"R1_{barcode}"] = open(f"{output_dir}/R1_{barcode}.fastq", 'a')
        file_handles[f"R2_{barcode}"] = open(f"{output_dir}/R2_{barcode}.fastq", 'a')
    file_handles["unknown_R1"] = open(f"{output_dir}/unknown_R1.fastq", 'a')
    file_handles["unknown_R2"] = open(f"{output_dir}/unknown_R2.fastq", 'a')
    file_handles["hopped_R1"] = open(f"{output_dir}/hopped_R1.fastq", 'a')
    file_handles["hopped_R2"] = open(f"{output_dir}/hopped_R2.fastq", 'a')
    return file_handles

def close_output_files(file_handles):
    '''Closes the files at the end'''
    for fh in file_handles.values():
        fh.close()

def demultiplex(r1_file, r2_file, r3_file, r4_file, output_dir, quality_threshold):
    '''This function reads through the four input FASTQ files containing the sequencing reads (2) and their
    corresponding barcodes. It classifies reads into matched, hopped, unknown, or low-quality
    categories based on barcodes and quality scores, then writes them to appropriate
    output files in the specified directory. Statistics of the demultiplexing process, including
    counts of each category, are maintained globally.'''

    global number_of_matched_reads, number_of_hopped_reads, number_of_unk_reads, number_of_low_quality_reads

    file_handles = open_output_files(output_dir)
    statistics_file = f"{output_dir}/Statistics.txt"

    with gzip.open(r1_file, "rt") as r1, gzip.open(r2_file, "rt") as r2, gzip.open(r3_file, "rt") as r3, gzip.open(r4_file, "rt") as r4:
        while True:
            # Read each record of each file
            r1_record = [r1.readline().strip() for _ in range(4)]
            r2_record = [r2.readline().strip() for _ in range(4)]
            r3_record = [r3.readline().strip() for _ in range(4)]
            r4_record = [r4.readline().strip() for _ in range(4)]

            if not r1_record[0] or not r2_record[0] or not r3_record[0] or not r4_record[0]:
                break

            r1_header = r1_record[0]
            r1_seq = r1_record[1]
            r1_plus = r1_record[2]
            r1_quality = r1_record[3]

            r2_seq = r2_record[1]

            r3_seq = reverse_complement(r3_record[1])
            
            r4_header = r4_record[0]
            r4_seq = r4_record[1]
            r4_plus = r4_record[2]
            r4_quality = r4_record[3]

            index_pair = f"{r2_seq}-{r3_seq}"
            r1_new_header = f"{r1_header} {index_pair}"
            r4_new_header = f"{r4_header} {index_pair}"

            # Check if R2 and R3 are known barcodes and pass the quality score threshold
            r2_quality_score = bioinfo.qual_score(r2_record[3])
            r3_quality_score = bioinfo.qual_score(r3_record[3])

            if r2_seq in known_barcodes and r3_seq in known_barcodes:  # Check if both barcodes are valid
                if r2_quality_score >= quality_threshold and r3_quality_score >= quality_threshold:
                    if r2_seq == r3_seq:
                        number_of_matched_reads += 1
                        if r2_seq not in number_of_reads_per_index:
                            number_of_reads_per_index[r2_seq] = 0
                        number_of_reads_per_index[r2_seq] += 1

                        write_fastq(file_handles[f"R1_{r2_seq}"], [r1_new_header, r1_seq, "+", r1_quality])
                        write_fastq(file_handles[f"R2_{r2_seq}"], [r4_new_header, r4_seq, "+", r4_quality])
                    else:
                        number_of_hopped_reads += 1
                        write_fastq(file_handles["hopped_R1"], [r1_new_header, r1_seq, r1_plus, r1_quality])
                        write_fastq(file_handles["hopped_R2"], [r4_new_header, r4_seq, r4_plus, r4_quality])
                
                else:
                    number_of_low_quality_reads += 1
            
            else:
                    number_of_unk_reads += 1
                    write_fastq(file_handles["unknown_R1"], [r1_new_header, r1_seq, r1_plus, r1_quality])
                    write_fastq(file_handles["unknown_R2"], [r4_new_header, r4_seq, r4_plus, r4_quality])


    # Writing statistics file
    total_number_of_reads = number_of_matched_reads + number_of_hopped_reads + number_of_unk_reads + number_of_low_quality_reads
    perc_matched_reads = (number_of_matched_reads / total_number_of_reads) * 100
    perc_of_hopped_reads = (number_of_hopped_reads / total_number_of_reads) * 100
    perc_of_unk_reads = (number_of_unk_reads/total_number_of_reads) * 100
    percentage_of_low_quality = (number_of_low_quality_reads/total_number_of_reads) * 100 


    with open(statistics_file, 'w') as stats:
        '''opens the statistics file for writing, and writes all the calculated values'''
        stats.write(f"Total number of reads: {total_number_of_reads}\n")
        stats.write(f"Number of matched reads: {number_of_matched_reads}, percentage of matched reads: {perc_matched_reads}%\n")
        stats.write(f"Number of hopped reads: {number_of_hopped_reads}, percentage of hopped reads: {perc_of_hopped_reads}%\n")
        stats.write(f"Number of unknown reads: {number_of_unk_reads},  percentage of unknown reads: {perc_of_unk_reads}%\n")
        stats.write(f"Number of low quality reads: {number_of_low_quality_reads}, percentage of low quality:{ percentage_of_low_quality}%\n")
        stats.write("Number of reads per index:\n")
        for index, count in number_of_reads_per_index.items():
            stats.write(f"{index}: {count}, Percentage: {count/(number_of_matched_reads)*100}%\n")

    close_output_files(file_handles) #closes all the files after writing 

if __name__ == "__main__":
    demultiplex(args.r1_file, args.r2_file, args.r3_file, args.r4_file, args.output_dir, args.quality_threshold)

