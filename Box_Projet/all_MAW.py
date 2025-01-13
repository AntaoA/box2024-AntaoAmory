#Enumerate all minimal absent words for each DNA sequence shorter than a given length in a FASTA file

import pysais
import sys
from suffix_tree import Tree

# Read a FASTA file and return a list of sequences and a list of sequence names
def fasta_to_seq_list(fasta_file, nb_seq_max):
    #Reads a FASTA file and returns a list of sequences
    seq_list = []
    seq_name = []
    i = 0
    with open(fasta_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if i == nb_seq_max:
                break
            if line[0] == '>':
                seq_name.append(line[1:].strip())
                seq_list.append('')
            else:
                seq_list[-1] += line.strip()
                i += 1
    return seq_list, seq_name

                
# Reverse complement of a DNA sequence
def reverse_complement(seq):
    #Returns the reverse complement of a DNA sequence
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in seq[::-1]])


# Build the suffix array of sequence
def suffix_array(sequence):
    return pysais.sais(sequence)

# Build the suffix tree of sequence
def suffix_tree(sequence):
    t = Tree()
    t.add(1, sequence)
    return t

# Search a pattern in the suffix array of sequence
def search_pattern(sequence, pattern, suffix_tree):
    return suffix_tree.find(pattern)

# Find the next letter over the alphabet {A, C, G, T}
def next_letter(letter):
    if letter == 'A':
        return 'C'
    elif letter == 'C':
        return 'G'
    elif letter == 'G':
        return 'T'
    else:
        return 'A'
    

# Find all MA of a DNA sequence between two suffixes in the suffix array
def potential_MAW(init, goal, max_length):
    progress = init
    pot_MAW = []
    pos_letter = len(init) - 1
    if pos_letter > max_length:
        progress = progress[:max_length]
        pos_letter = max_length - 1
    flag_init = True # mark the first time we analyse init (in order to avoid infinite loop)
    while progress != goal:
        if pos_letter <= len(goal) and progress[:pos_letter+1] == goal[:pos_letter+1]:
            if pos_letter > max_length:
                break
            pos_letter += 1
            if goal[pos_letter] == 'A':
                progress += 'A'
            else:
                pot_MAW.append(''.join([progress,'A']))
                if goal[pos_letter] == 'C':
                    progress += 'C'
                else:
                    pot_MAW.append(''.join([progress + 'C']))
                    if goal[pos_letter] == 'G':
                        progress += 'G'
                    else:
                        pot_MAW.append(''.join([progress + 'G']))
                        progress += 'T'
                
        else:
            if init == progress and flag_init:
                flag_init = False
                progress += 'A'
                pot_MAW.append(progress)
                pos_letter += 1
            else:
                if progress[-1] == 'T':
                    progress = progress[:-1]
                    pos_letter -= 1  
                else:       
                    progress = progress[:-1] + next_letter(progress[-1])
                    if progress != goal[:pos_letter+1]: 
                        pot_MAW.append(progress)
    return pot_MAW

# Find all potential MAW of a DNA sequence
def all_potential_MAW(sequence, suffix_array, max_length):
    all_pot_MAW = []
    init = ""
    goal = sequence[suffix_array[0]:]
    all_pot_MAW.extend(potential_MAW(init, goal, max_length))
    goal = "error"
    for i in range(len(suffix_array)-1):
        init = sequence[suffix_array[i]:]
        goal = sequence[suffix_array[i+1]:]
        all_pot_MAW.extend(potential_MAW(init, goal, max_length))
    init = goal
    goal = ""
    all_pot_MAW.extend(potential_MAW(init, goal, max_length))
    all_pot_MAW = [maw for maw in all_pot_MAW if len(maw) <= max_length]
    return all_pot_MAW


# Verify if a pattern is a MAW of a DNA sequence
def verify_MAW(sequence, maw, max_length, st):
    if len(maw) > max_length:
        return False
    if len(maw) == 1:
        return False
    if search_pattern(sequence, reverse_complement(maw), st):
        return False

    pos_1_pattern_found = search_pattern(sequence, maw[:-1], st)
    pos_2_pattern_found = search_pattern(sequence, maw[1:], st)
    pos_1_rev_comp_found = search_pattern(sequence, reverse_complement(maw[:-1]), st)
    pos_2_rev_comp_found = search_pattern(sequence, reverse_complement(maw[1:]), st)
    
    if (pos_1_pattern_found or pos_1_rev_comp_found) and (pos_2_pattern_found or pos_2_rev_comp_found):
        return True
    else:
        return False

# TSV : sequence name, k, number of MAW of length k, MAW of length k
def write_MAW_to_TSV(fasta_file, max_length, nb_seq_max, output_file):
    seq_list, seq_name = fasta_to_seq_list(fasta_file, nb_seq_max)
    nb_seq = len(seq_list)
    print('Number of sequences:', nb_seq)
    print('Number of bases:', sum([len(seq) for seq in seq_list]))
    with open(output_file, 'w') as f:
        f.write('Sequence\tk\tnumber of MAW\tMAW\n')
        for i_seq in range(len(seq_list)):
            seq = seq_list[i_seq]
            name = seq_name[i_seq]
            print('Sequence', i_seq+1, 'out of', nb_seq, '(nb of bases:', len(seq), ')')
            if not set(seq) <= {'A', 'C', 'G', 'T'}:
                #print('Error: the sequence contains at least one base different from A, C, G, T : ', set(seq))
                continue
            sa = suffix_array(seq)
            st = suffix_tree(seq)
            all_pot_MAW = all_potential_MAW(seq, sa, max_length)
            all_MAW = set()
            for maw in all_pot_MAW:
                if verify_MAW(seq, maw, max_length, st):
                    all_MAW.add(maw)
                    all_MAW.add(reverse_complement(maw))
            sort_MAW_by_length = sorted(all_MAW, key = len)
            i = 0
            for k in range(2, max_length+1):
                length_k_MAW = []
                while(i < len(sort_MAW_by_length) and len(sort_MAW_by_length[i]) == k):
                    length_k_MAW.append(sort_MAW_by_length[i])
                    i += 1
                sorted_length_k_MAW = sorted(length_k_MAW)
                if sorted_length_k_MAW:
                    f.write(name + '\t' + str(k) + '\t' + str(len(sorted_length_k_MAW)) + '\t' + ','.join(sorted_length_k_MAW) + '\n')
    
if __name__ == '__main__':
    if len(sys.argv) != 4:
        print('Usage: python all_MAW.py <fasta_file> <max_length> <nb_seq_max>')
        sys
    fasta_file = sys.argv[1]
    max_length = int(sys.argv[2])
    nb_seq_max = int(sys.argv[3])
    output_file = 'MAW - ' + sys.argv[2] + ' ' + sys.argv[3] +'.tsv'
    write_MAW_to_TSV(fasta_file, max_length, nb_seq_max, output_file)