from math import inf

def get_prefix(text):
    return (text[:i] for i in range(len(text)+1))


def gc_skew(text, char1="G", char2="C"):
    char1_count = text.count(char1)
    char2_count = text.count(char2)
    return char1_count - char2_count

def minimize_skew(skews):
    min_skew = min(skews)
    return [i for i, j in enumerate(skews) if j == min_skew]


def skew_efficient(text):
    skew_dict = {'A': 0, 'C': -1, 'G': 1, 'T': 0}
    min_skew_inds = []
    min_skew = inf
    current_skew = 0
    for i, j in enumerate(text):
        current_skew += skew_dict[j]
        if current_skew < min_skew:
            min_skew_inds = []
            min_skew = current_skew
        if current_skew == min_skew:
            min_skew_inds.append(i)

    return min_skew_inds


if __name__ == '__main__':
    genome = "CCTATCGGTGGATTAGCATGTCCCTGTACGTTTCGCCGCGAACTAGTTCACACGGCTTGATGGCAAATGGTTTTTCCGGCGACCGTAATCGTCCACCGAG"
    prefixes = get_prefix(genome)
    skews = [gc_skew(i) for i in prefixes]
    min_skew = minimize_skew(skews)
    print(' '.join([str(i) for i in min_skew]))


    user_path = "data/BA1F.txt"
    with open(user_path, 'r') as f_r:
        genome = f_r.read().strip()

    print(' '.join([str(i) for i in skew_efficient(genome)]))
