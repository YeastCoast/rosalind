from itertools import product

def reverse_comp(nucleotide_sequence):
    NUCLEOTIDE_DICT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    output = [NUCLEOTIDE_DICT[i] for i in nucleotide_sequence]
    output.reverse()
    return ''.join(output)


def get_kmers(text, k):
    output = {}
    for i in range(len(text)-k):
        curr_kmer = text[i:i+k]
        if curr_kmer not in output:
            output[curr_kmer] = 0
        output[curr_kmer] += 1
    return output


def hamming_distance(str1, str2):
    return sum(1 for i, j in zip(str1, str2) if i != j)


def kmer_count_with_mismatches_and_reverse(kmers: dict, mismatches, k):
    all_kmers = product(["A", "C", "G", "T"], repeat=k)
    output_dict = {}
    for kmer in all_kmers:
        output_dict[''.join(kmer)] = sum(j for i, j in kmers.items() if (hamming_distance(i, kmer) <= mismatches)) + \
                            sum(j for i, j in kmers.items() if (hamming_distance(reverse_comp(i), kmer) <= mismatches))

    return output_dict


def get_most_frequent(kmers: dict):
    max_val = max(kmers.values())
    return ' '.join([i for i, j in kmers.items() if j == max_val])


if __name__ == '__main__':
    # test case
    genome = "ACGTTGCATGTCGCATGATGCATGAGAGCT"
    parameter = "4 1"

    k, d = [int(i) for i in parameter.split(' ')]
    kmers = get_kmers(genome, k)
    kmers_mismatch = kmer_count_with_mismatches_and_reverse(kmers, d, k)
    most_frequent_kmers = get_most_frequent(kmers_mismatch)
    print(most_frequent_kmers)

    # user case example
    genome = "ATAACACGTGCGCGACGATAACACGATTCCAGGTGCGCGACGATTCCAGGATAACACGATAACACGCTCGGGGATTCCAGGATAACACGCTCGGGGATTCCAGGCTCGGGGATTCCAGGGCTGAGAGGGCTGAGAGGCTCGGGGCTCGGGGCTCGGGGATTCCAGGGCTGAGAGGCTCGGGGCTCGGGGGCTGAGAGGATTCCAGGATTCCAGGATAACACGCTCGGGGCTCGGGGTGCGCGACGGCTGAGAGGGCTGAGAGGATAACACGCTCGGGGGCTGAGAGGATTCCAGGGCTGAGAGGGCTGAGAGGCTCGGGGATAACACGCTCGGGGATTCCAGGCTCGGGGATTCCAGGGCTGAGAGGATTCCAGGATTCCAGGATTCCAGGGCTGAGAGGATAACACGCTCGGGGGCTGAGAGGCTCGGGGATAACACGCTCGGGGATAACACGATAACACGATTCCAGGTGCGCGACGCTCGGGGATTCCAGGGCTGAGAGGATAACACGGCTGAGAGGGCTGAGAGGATTCCAGGCTCGGGGATAACACGTGCGCGACGATTCCAGGCTCGGGGCTCGGGGCTCGGGGCTCGGGGCTCGGGGATTCCAGGATTCCAGGGCTGAGAGGATTCCAGGATTCCAGGCTCGGGGTGCGCGACGATAACACGATAACACGGCTGAGAGGCTCGGGGTGCGCGACGGCTGAGAGGTGCGCGACGGCTGAGAGGCTCGGGGTGCGCGACGGCTGAGAGGATTCCAGGCTCGGGGCTCGGGGATTCCAGGGCTGAGAGGTGCGCGACG"
    parameter = "5 2"

    k, d = [int(i) for i in parameter.split(' ')]
    kmers = get_kmers(genome, k)
    kmers_mismatch = kmer_count_with_mismatches_and_reverse(kmers, d, k)
    most_frequent_kmers = get_most_frequent(kmers_mismatch)
    print(most_frequent_kmers)
