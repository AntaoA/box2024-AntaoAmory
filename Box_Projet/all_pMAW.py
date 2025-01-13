# Import TSV file : sequence name, k, number of MAW of length k, MAW of length k

import sys
import math
from bitarray import bitarray


# Reverse complement of a DNA sequence
def reverse_complement(seq):
    #Returns the reverse complement of a DNA sequence
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in seq[::-1]])

# Create a dictionary with the MAW of each sequence sorted by length and add a flag for MAW presence
def create_dict_per_length(dict_MAW, flag_seq):
    dict_length = {}
    for i,name in enumerate(dict_MAW):
        for k in dict_MAW[name]:
            if k not in dict_length:
                dict_length[k] = {}
            for maw in dict_MAW[name][k]:
                if maw == min(maw, reverse_complement(maw)):    # to avoid duplicates since maw and its reverse complement are both in the dictionary
                    if maw in dict_length[k]:
                        flag = dict_length[k][maw].copy()
                        flag[i] = 1
                        dict_length[k][maw] = flag
                    else:
                        flag = flag_seq.copy()
                        flag[i] = 1
                        dict_length[k][maw] = flag.copy()
    #sort the dictionary on both keys
    dict_length = {k: dict_length[k] for k in sorted(dict_length)}
    for k in dict_length:
        dict_length[k] = {maw: dict_length[k][maw] for maw in sorted(dict_length[k])}    
    return dict_length


# Precompute the number of sequences that contain a pattern and its reverse complement, for MAW
def pre_count_AW(dict_length):
    cache = {}
    for k in dict_length:
        for pattern in dict_length[k]:
            cache[pattern] = dict_length[k][pattern].copy()
    for k in dict_length:
        for pattern in dict_length[k]:
            for i in range(k):
                for j in range(i+1, k+1):
                    pat = min(pattern[i:j], reverse_complement(pattern[i:j]))
                    if pat in cache:
                        cache[pattern] |= cache[pat]
    return cache



# Count the number of sequences that contain a pattern and its reverse complement
def count_AW(dict_length, pattern, flag, cache):
    if pattern in cache:
        return cache[pattern].count()
    
    cache[pattern] = flag.copy()
   
    k = len(pattern)
    for i in range(min(dict_length), k+1):
        for j in range(k-i+1):
            pat = min(pattern[j:j+i], reverse_complement(pattern[j:j+i]))
            if pat in cache:
                cache[pattern] |= cache[pat]

    return cache[pattern].count()


# Add patterns to the potential pMAW
def add_pattern(dict_length, dict_pot, pattern, n, max_length, cache, nb_seq, nb_test, nb_test_valid, nb_new_test, nb_new_test_valid):
    k = len(pattern)
    if k+1 <= max_length:
        nb_test += 1
        if count_AW(dict_length, pattern, bitarray("0" * nb_seq), cache) < n:
            nb_test_valid += 1
            nb_new_test += 4
            if k+1 not in dict_pot:
                dict_pot[k+1] = set()
            pattern_A = min(pattern + "A", reverse_complement(pattern + "A"))
            pattern_C = min(pattern + "C", reverse_complement(pattern + "C")) 
            pattern_G = min(pattern + "G", reverse_complement(pattern + "G"))
            pattern_T = min(pattern + "T", reverse_complement(pattern + "T"))
            if pattern_A not in dict_pot[k+1]:
                nb_new_test_valid += 1
                dict_pot, nb_test, nb_test_valid, nb_new_test, nb_new_test_valid = add_pattern(dict_length, dict_pot, pattern_A, n, max_length, cache, nb_seq, nb_test, nb_test_valid, nb_new_test, nb_new_test_valid)
            if pattern_C not in dict_pot[k+1]:
                nb_new_test_valid += 1
                dict_pot,nb_test, nb_test_valid, nb_new_test, nb_new_test_valid = add_pattern(dict_length, dict_pot, pattern_C, n, max_length, cache, nb_seq, nb_test, nb_test_valid, nb_new_test, nb_new_test_valid)
            if pattern_G not in dict_pot[k+1]:
                nb_new_test_valid += 1
                dict_pot,nb_test, nb_test_valid, nb_new_test, nb_new_test_valid = add_pattern(dict_length, dict_pot, pattern_G, n, max_length, cache, nb_seq, nb_test, nb_test_valid, nb_new_test, nb_new_test_valid)
            if pattern_T not in dict_pot[k+1]:
                nb_new_test_valid += 1
                dict_pot,nb_test, nb_test_valid, nb_new_test, nb_new_test_valid = add_pattern(dict_length, dict_pot, pattern_T, n, max_length, cache, nb_seq, nb_test, nb_test_valid, nb_new_test, nb_new_test_valid)
        else:
            if k not in dict_pot:
                dict_pot[k] = set()
            dict_pot[k].add(pattern)

    return dict_pot, nb_test, nb_test_valid, nb_new_test, nb_new_test_valid

# Find all potential pMAW
def all_potential_pMAW(dict_length, n, max_length, cache, nb_seq, nb_test, nb_test_valid, nb_new_test, nb_new_test_valid):
    dict_pot_pMAW = {}
    for k in range(max_length+1):
        if k in dict_length: 
            for pattern in dict_length[k]:
                dict_pot_pMAW, nb_test, nb_test_valid, nb_new_test, nb_new_test_valid = add_pattern(dict_length, dict_pot_pMAW, pattern, n, max_length, cache, nb_seq, nb_test, nb_test_valid, nb_new_test, nb_new_test_valid)
                
    return dict_pot_pMAW, nb_test, nb_test_valid, nb_new_test, nb_new_test_valid

# Verify if a pattern is a pMAW
def verify_pMAW(dict_length, pattern, n, nb_seq, cache):
    left = min(pattern[:-1], reverse_complement(pattern[:-1]))
    right = min(pattern[1:], reverse_complement(pattern[1:]))
    count_patter = count_AW(dict_length, left, bitarray("0" * nb_seq), cache)
    count_attern = count_AW(dict_length, right, bitarray("0" * nb_seq), cache)
    if not count_patter >= n and not count_attern >= n:
        return True
    else:
        return False
    
# Find all pMAW
def all_pMAW(dict_length, dict_pot_pMAW, n, nb_seq, cache):
    dict_pMAW = {}
    for k in dict_pot_pMAW:
        for pattern in dict_pot_pMAW[k]:
            if verify_pMAW(dict_length, pattern, n, nb_seq, cache):
                if k not in dict_pMAW:
                    dict_pMAW[k] = set()
                dict_pMAW[k].add(pattern)
    return dict_pMAW


# TSV : k, number of pMAW of length k, pMAW of length k
def write_MAW_to_TSV(dict_pMAW, max_length):
    with open("pMAW.tsv", "w") as f:
        f.write("k\tnumber of pMAW of length k\tpMAW of length k\n")
        for i in range(max_length):
            if i in dict_pMAW:
                dict_pMAW[i] = sorted(dict_pMAW[i])
                f.write(str(i) + "\t" + str(len(dict_pMAW[i])) + "\t" + ",".join(dict_pMAW[i]) + "\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print('Usage: python all_pMAW.py <fasta_file> <p> <max_length>')
        sys
    file = sys.argv[1]
    p = float(sys.argv[2])
    if p <= 0 or p > 1:
        print("p must be greater than 0 and less than or equal to 1")
        sys.exit(1)
    max_length = int(sys.argv[3])
    dict_MAW = {}
    with open(file, "r") as f:
        lines = f.readlines()
        for i,line in enumerate(lines):
            if i == 0:
                continue
            name, k, _, list_MAW = line.split("\t")
            list_MAW = list_MAW[:-1]    # remove the last character which is a newline
            k = int(k)
            list_MAW = list_MAW.split(",")
            if name not in dict_MAW:
                dict_MAW[name] = {}
            (dict_MAW[name])[k] = list_MAW
    nb_seq = len(dict_MAW)
    n = math.ceil(p * nb_seq)
    flag = bitarray("0" * nb_seq)
    nb_test = 0
    nb_test_valid = 0
    nb_new_test = 0
    nb_new_test_valid = 0

    dict_length = create_dict_per_length(dict_MAW, flag)
    cache = pre_count_AW(dict_length)
    dict_pot_pMAW, nb_test, nb_test_valid, nb_new_test, nb_new_test_valid = all_potential_pMAW(dict_length, n, max_length, cache, nb_seq, nb_test, nb_test_valid, nb_new_test, nb_new_test_valid)
    dict_pMAW = all_pMAW(dict_length, dict_pot_pMAW, n, nb_seq, cache)
    print("gamma:", nb_test_valid/nb_test)
    print("delta:", nb_new_test_valid/nb_new_test)
    write_MAW_to_TSV(dict_pMAW, max_length)