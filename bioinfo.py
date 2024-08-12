#!/usr/bin/env python

# Author: Mahmoud Al Mahmoud

#     Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework. Started on 07/2024'''

__version__ = "0.4"    
                        
DNA_bases = ('ATCGNatcgn')
RNA_bases = ('AUCGNaucgn') 

def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score'''
    qual_score = ord(letter) - 33 
    return qual_score 

def qual_score(phred_score: str) -> float:
    '''Converts Phred scores from letters to their corresponding numbers'''
    score_sum = 0
    for val in phred_score: 
        score_sum += convert_phred(val)
        average = score_sum/len(phred_score)
    return(average) 

def validate_base_seq(base: str, RNAFLAG: bool =False):
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    bases = set(base)
    return bases.issubset(RNA_bases) if RNAFLAG else bases.issubset(DNA_bases)


def gc_content(seq: str):
    '''Returns GC content of a DNA or RNA sequence as a decimal between 0 and 1.'''
    assert validate_base_seq(seq)
    seq.upper()
    count = seq.count("G") + seq.count("C") 
    total = len(seq)
    return count/total     
  
def calc_median(lst: list)-> float:
    ''''Calculates and returns the median of a one-dimensional list that is already sorted.'''
    mid = (len(lst) - 1) // 2
    if len(lst) % 2 == 0:
        median = (lst[mid] + lst[mid+1]) / 2   
    else: 
        median = lst[mid]   
    return median

def oneline_fasta(file1, file2):
    '''Writes out each fasta sequence line as one line per record maintaining the fasta format'''
    with open(file1, "r") as fh1, open(file2, "w") as fh2: 
        current_seq = ""
    while True:
        line = fh1.readline() 
        if not line:
            fh2.write(current_seq)
            break
        elif line[0] == ">": 
            if len(current_seq) > 0: 
                fh2.write(current_seq + '\n')
            fh2.write(line)
            current_seq = ""
        else:
            current_seq = current_seq + line.strip() 
    return file2

if __name__ == "__main__":
    # write tests for functions above, Leslie has already populated some tests for convert_phred
    # These tests are run when you execute this file directly (instead of importing it)
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Your convert_phred function is working! Nice job")

    assert calc_median([1,2,100]) == 2, "calc_median function does not work for odd length list"
    assert calc_median([1,2]) == 1.5, "calc_median function does not work for even length list"
    print("Median successfully calculated")

    assert validate_base_seq("AATAGAT"), "Validate base seq does not work on DNA"
    assert validate_base_seq("AAUAGAU", True), "Validate base seq does not work on RNA"
    assert validate_base_seq("R is the best!")==False, "Not a DNA string"
    assert validate_base_seq("aatagat"), "Validate base seq does not work on lowercase DNA"
    assert validate_base_seq("aauagau", True), "Validate base seq does not work on lowercase RNA"
    assert validate_base_seq("TTTTtttttTTT")
    print(" validate_base_seq successfully calculated")

    assert gc_content("GCGCGC") == 1
    assert gc_content("AATTATA") == 0
    assert gc_content("GCATCGAT") == 0.5
   
    assert qual_score("A") == 32.0, "wrong average phred score for 'A'"
    assert qual_score("AC") == 33.0, "wrong average phred score for 'AC'"
    assert qual_score("@@##") == 16.5, "wrong average phred score for '@@##'"
    assert qual_score("EEEEAAA!") == 30.0, "wrong average phred score for 'EEEEAAA!'"
    assert qual_score("$") == 3.0, "wrong average phred score for '$'"