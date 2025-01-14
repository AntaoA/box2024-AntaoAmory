#! /usr/bin/env python3

import sys
import Levenshtein      # huge time saver


"""
    # Function calculating distance by memoization. Abandoned because Levenstein is MUCH more optimized.

def edDistRecursiveMemo1(a, b, memo, max_errors=float("inf"), current_score=0):

    if len(a) == 0:
        return len(b) + current_score
    if len(b) == 0:
        return len(a) + current_score

    if current_score > max_errors:
        return float("inf")

    key = (len(a), len(b))
    if key in memo:
        return memo[key] + current_score

    delta = 1 if a[-1] != b[-1] else 0
    substitution = edDistRecursiveMemo(a[:-1], b[:-1], memo, max_errors, current_score + delta)
    insertion = edDistRecursiveMemo(a[:-1], b, memo, max_errors, current_score + 1)
    deletion = edDistRecursiveMemo(a, b[:-1], memo, max_errors, current_score + 1)

    result = min(substitution, insertion, deletion)

    memo[key] = result - current_score
    return result
"""



def parse_fasta(file_path):
    """
    Reads a FASTA file and returns a dictionary {name: sequence}.
    """
    sequences = {}
    name = None
    sequence = []
    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):                        # Description line
                if name:
                    sequences[name] = "".join(sequence)
                name = line[1:]                             # Remove ">"
                sequence = []
            else:
                sequence.append(line)                       # Sequence line
        if name:                                            # Add the last sequence
            sequences[name] = "".join(sequence)
    return sequences


def parse_fastq(file_path):
    """
    Reads a FASTQ file and returns a list of tuples (name, sequence).
    """
    reads = []
    with open(file_path, "r") as file:
        while True:
            name = file.readline().strip()      # Line 1: Name like @ERR...
            if not name:
                break
            sequence = file.readline().strip()  # Line 2: Sequence
            file.readline()                     # Line 3: Ignored -- AAAAAAEEEEEE
            file.readline()                     # Line 4: Ignored -- +
            reads.append((name[1:], sequence))  # Remove "@"
    return reads


def reverse_complement(seq):
    """
    Returns the reverse complement of a DNA sequence, including ambiguous bases.
    """
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return "".join(complement.get(base, 'N') for base in reversed(seq))


def build_kmer_index(references, k):
    """
    Builds a dictionary of k-mers from multiple reference sequences.
    Returns a dictionary where each k-mer points to (seq_name, position).
    """
    index = {}
    for seq_name, sequence in references.items():
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i + k]
            if kmer not in index:
                index[kmer] = []
            index[kmer].append((seq_name, i))
    return index


def edit_distance(a, b, max_errors=float("inf")):
    """
    Calcul rapide de la distance d'édition avec la bibliothèque Levenshtein.
    """
    score = Levenshtein.distance(a, b)
    operations = Levenshtein.editops(a, b)
    return (score, operations) if score <= max_errors else (float("inf"), [])


def transcript_of_edits(operations, a, b):
    """
    Returns the transcript of edit operations from a to b.
    """
    transcript = []
    for operation in operations:
        if operation[0] == "replace":
            transcript.append(f"R{operation[2]}")
        elif operation[0] == "delete":
            transcript.append(f"D{operation[1]}")
        elif operation[0] == "insert":
            transcript.append(f"I{operation[2]}")
    return transcript


def align_reads(references, reads, k, max_errors, show_progress=False):
    """
    Aligns reads against the reference and returns the best alignments.
    """
    kmer_index = build_kmer_index(references, k)
    results = []

    total_reads = len(reads)
    for read_count, (read_name, read) in enumerate(reads):
        if show_progress and read_count % 1000 == 0:
            print(f"Alignment progress: {read_count / total_reads * 100:.2f}%")
        best_alignments = []
        rev_comp = reverse_complement(read)

        # Avoid re-aligning the same region
        processed_positions = {seq_name: set() for seq_name in references.keys()}

        # Align for both orientations (forward and reverse complement)
        for seq in [read, rev_comp]:
            for i in range(len(seq) - k + 1):
                kmer = seq[i:i + k]
                if kmer in kmer_index:
                    for ref_name, ref_pos in kmer_index[kmer]:
                        if ref_pos in processed_positions[ref_name]:
                            continue

                        reference = references[ref_name]

                        # Extend the reference region to match the read
                        start_pos = max(0, ref_pos - i)
                        end_pos = min(len(reference), start_pos + len(seq))
                        ref_seq = reference[start_pos:end_pos]

                        # Calculate edit distance
                        score, operations = edit_distance(seq, ref_seq, max_errors)
                        change = 0


                        # Extend the reference region to match the read
                        for operation in operations:
                            if operation[1] == 0:
                                if operation[0] == "delete":
                                    start_pos += 1
                                    change = 1
                                elif operation[0] == "insert":
                                    start_pos -= 1
                                    change = 1
                            elif operation[2] >= len(ref_seq) - 1:
                                if operation[0] == "delete":
                                    end_pos += 1
                                    change = 1
                                elif operation[0] == "insert":
                                    end_pos -= 1
                                    change = 1

                        if change:
                            ref_seq = reference[start_pos:end_pos]
                            score, operations = edit_distance(seq, ref_seq, max_errors)


                        # Update best alignments
                        if score <= max_errors:
                            transcription = transcript_of_edits(operations, seq, ref_seq)
                            strand = "F" if seq == read else "R"
                            best_alignments.append((
                                ref_name, read_name, strand, start_pos,
                                start_pos + len(seq) - 1, score, transcription
                            ))

                        # Mark processed positions
                        processed_positions[ref_name].update(range(start_pos, end_pos))


        # Keep only the best alignments
        best_alignments.sort(key=lambda x: x[5])  # Sort by score
        if best_alignments:
            best_score = best_alignments[0][5]
            best_alignments = [a for a in best_alignments if a[5] == best_score]
            results.extend(best_alignments)
    return results





def write_results(results, output_file):
    """
    Writes the alignment results to a tab-separated file.
    """
    with open(output_file, "w") as file:
        for res in results:
            file.write("\t".join(map(str, res)) + "\n")


def main():
    # Read command-line arguments
    args = sys.argv[1:]
    if len(args) < 3:
        print("Usage: python script.py <reference_file> <reads_file> <k> [max_errors=4] [--p]   # Optional: --p for showing progress")
        sys.exit(1)

    reference_file = args[0]
    reads_file = args[1]
    k = int(args[2])
    max_errors = 4  # Default value
    show_progress = False

    for arg in args[3:]:
        if arg.startswith("--p"):
            show_progress = True
        else:
            try:
                max_errors = int(arg)
            except ValueError:
                print(f"Invalid argument: {arg}")
                sys.exit(1)

    # Processing
    print("Reading files...")
    references = parse_fasta(reference_file)
    print(f"Reference: {len(references)} reference sequences")
    reads = parse_fastq(reads_file)
    print(f"Reads: {len(reads)} reads")

    print("Starting alignment...")
    results = align_reads(references, reads, k, max_errors, show_progress)
    print(f"{len(results)} alignments found")
    print(f"{len(results)/len(reads)*100:.2f}% of reads aligned")

    output_file = "results.tsv"
    write_results(results, output_file)
    print(f"Results written to {output_file}")


if __name__ == "__main__":
    main()
