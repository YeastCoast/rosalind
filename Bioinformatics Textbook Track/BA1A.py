
def pattern_count(text, pattern):
    """
    simple pattern count algorithm, find the number of occurrences of the pattern in the text;
    iterates over the text and checks at each position if the pattern is found by a second for loop
    inefficient with a time complexity of O(n*m)
    :param text:
    :param pattern:
    :return:
    """
    count = 0
    for i in range(len(text)-len(pattern)+1):
        current_count = 0
        for j in range(len(pattern)):
            if text[i+j] != pattern[j]:
                break
            else:
                current_count += 1
        if current_count == len(pattern):
            count += 1

    return count


def pythonic_pattern_count(text: str, pattern: str):
    return sum(1 if text[i: i+len(pattern)] == pattern else 0 for i in range(len(text)))


def kmp(text, pattern):
    lps_table = [0 for _ in pattern]
    pointer = 0
    counter = 0
    for i in range(1, len(pattern)):
        if pattern[i] != pattern[pointer]:
            counter = 0
            pointer = 0
        if pattern[i] == pattern[pointer]:
            counter += 1
            pointer += 1
            lps_table[i] = counter

    i = 0
    j = 0
    pattern_counter = 0
    while len(text)-i >= len(pattern) - j:
        if text[i] == pattern[j]:
            i += 1
            j += 1
        if j == len(pattern):
            pattern_counter += 1
            j = lps_table[j-1]
        elif i < len(text) and text[i] != pattern[j]:
            if j == 0:
                i += 1
            else:
                j = lps_table[j-1]

    return pattern_counter


if __name__ == '__main__':
    # test case
    text = "GCGCG"
    pattern = "GCG"
    print(pattern_count(text, pattern))
    print(pythonic_pattern_count(text, pattern))
    print(kmp(text, pattern))

    # user case example
    text = "CGGTGGCTCCGAGTGCTCATCCGAGTTCCGAGTATCCTTCCGAGTTCCGAGTTCTTCCGAGTTCCGAGTCTGCTCCGAGTTCCGAGTACTCCGAGTCGAGCTGCCTCATCCGAGTTTAGACTCCGAGTTCCGAGTTTCCGAGTGTGGCGTCCGAGTGGTTCCGAGTTATCTCCGAGTGCGTCCGAGTAGATCCGAGTTCCGAGTAAGTCCGAGTTTTCCGAGTACAGTTACACAAATCCGAGTCGTTCCGAGTACTCCGAGTGTCCGAGTTGTTCCGAGTACGTCCGAGTGCCTTTCTGTTCACGAATTCCGAGTCTAGTCCGAGTTCCGAGTTCCGAGTTTCCATGGATCCGAGTTCCGAGTTCCGAGTCGTCCGAGTTCCGAGTTCCGAGTTCCGAGTTCCGAGTCCTCCGAGTATCCGAGTGATCCGAGTTCCGAGTAAATCCGAGTTCCTCCGAGTTTCCGAGTCGAAGCGGTCCGAGTGATCCGAGTAGGGCACCATCCGAGTCTCCGAGTAATCCGAGTTCCGAGTCTTCTCCGAGTTCCGAGTGCCAAATTCCGAGTTCCGAGTTTGTCCGAGTAATTCAGATCCGAGTTTCCGAGTTCCGAGTGTCCGAGTCTTAAGCTTGCTCCGAGTTCCGAGTTCCGAGTTACTGTGCGGCTCGTCCGAGTATCCGAGTTCTCCGAGTTTCCGAGTGTCCGAGTGGTCCGAGTGTCCGAGTAATCATCCGAGTGTCCGAGTTGTCCGAGTACGGTCCGAGTGTCCGAGTCTCCGAGTACAAGACCACGTAGGCACTCCGAGTTCCGAGTCTTTCCGAGTTCCGAGTACAGAACTCCGAGTCTCCGAGTTCATCCGAGTTTTCCGAGTGCCCGATCCGAGTCCTCCGAGTCCATCCGAGTTCCGAGTGGTTCCGAGTATCCGAGTGTCCGAGT"
    pattern = "TCCGAGTTC"
    print(pattern_count(text, pattern))
    print(pythonic_pattern_count(text, pattern))
    print(kmp(text, pattern))
