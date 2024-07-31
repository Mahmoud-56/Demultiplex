#!/usr/bin/env python 

# import argparse
# import bioinfo
# import matplotlib.pyplot as plt

# def get_args():
#     parser = argparse.ArgumentParser(description="Quality score per base distribution")
#     parser.add_argument("-f", "--filename", help="What is the file name", type=str)
#     parser.add_argument("-s", "--list_size", help="Size of the list (number of nucleotides in the read)", type=int)
#     return parser.parse_args()

# args = get_args()

# #file = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz"
# #size = 8

# def init_list(size, value=0.0):
#     return [value] * size



# #def summ_qscore_per_nuc(file):
# my_list = init_list(args.list_size)
# num_lines = 0
# with open(args.file) as fh:
#         for i, line in enumerate(fh):
#             if i % 4 == 3:  # Only process the quality score lines
#                 num_lines += 1
#                 newline = line.strip()
#                 for j, ch in enumerate(newline):
#                     qual_score = bioinfo.convert_phred(ch)
#                     my_list[j] += qual_score
                
# for i, summ in enumerate(my_list):
#     avg_score = summ / num_lines
#     my_list[i] = avg_score 

# plt.plot(my_list)
# plt.xlabel('Base pair')
# plt.ylabel('Average quality score')
# plt.title(f"Quality score distribution for Index1")
# plt.show()
# plt.savefig(f"quality_score_distribution_Index1.png")



import argparse
import bioinfo
import matplotlib.pyplot as plt
import sys

def get_args():
    parser = argparse.ArgumentParser(description="Quality score per base distribution")
    parser.add_argument("-s", "--list_size", type=int)
    parser.add_argument("-n", "--name", help="Name of the file for naming the plot", type=str)
    return parser.parse_args()

args = get_args()

def init_list(size, value=0.0):
    return [value] * size

my_list = init_list(args.list_size)
num_lines = 0

for i, line in enumerate(sys.stdin):
    if i % 4 == 3:  # Only process the quality score lines
        num_lines += 1
        newline = line.strip()
        for j, ch in enumerate(newline):
            qual_score = bioinfo.convert_phred(ch)
            my_list[j] += qual_score

N = num_lines
for i, summ in enumerate(my_list):
    avg_score = summ / N
    my_list[i] = avg_score 

x_values = range(len(my_list))
plt.bar(x_values, my_list)
plt.xlabel('Base pair')
plt.ylabel('Average quality score')
plt.title(f"Quality score distribution for {args.name}")
plt.savefig(f"quality_score_distribution_{args.name}.png")
