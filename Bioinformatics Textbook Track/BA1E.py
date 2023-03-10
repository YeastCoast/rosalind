def kmers(text, k):
    return [text[i:i+k] for i in range(len(text)-k)]


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
    indices = []
    while len(text)-i >= len(pattern) - j:
        if text[i] == pattern[j]:
            i += 1
            j += 1
        if j == len(pattern):
            indices.append(i-len(pattern))
            j = lps_table[j-1]
        elif i < len(text) and text[i] != pattern[j]:
            if j == 0:
                i += 1
            else:
                j = lps_table[j-1]

    return indices


def get_clumps(indices, length, threshold):
    if len(indices) == threshold:
        return True
    else:
        for i in indices:
            current_clump = len([j for j in indices if i <= j <= i+length])
            if current_clump >= threshold:
                return True
        return False


if __name__ == '__main__':
    # test case
    genome = "CGGACTCGACAGATGTGAAGAAATGTGAAGACTGAGTGAAGAGAAGAGGAAACACGACACGACATTGCGACATAATGTACGAATGTAATGTGCCTATGGC"
    parameters = "5 75 4"
    k, L, t = [int(i) for i in parameters.split(' ')]
    kmers_lst = kmers(genome, k)
    indices = {i: kmp(genome, i) for i in kmers_lst}
    potential_clumps = {i: indices[i] for i in indices if len(indices[i]) >= t}
    clumps = {i: get_clumps(potential_clumps[i], L, t) for i in potential_clumps}
    clumps = ' '.join({i for i, j in clumps.items() if j is True})
    print(clumps)

    # test case
    genome = "CCAGTGAAAATTAGGTCACAGAACTCAATTGACTTTGGATCTTCCGGCCGGCGCACACGCATTGCTACTAAGACTAAAACTATCTTGTATTGCTTGTCGCCCTTACCTTGGCACATCGCGTCATGGGCCTACTGCAACGGAGAAAGTATCTATCATAAGTCGCCAAAGACGAAAGCCACACATTGTGGTTGCCATCACTAAACTTAGACAAGCCCCGGACAACATACAGTGGATTATTAGCTTGATGGCCCCCCAGGCCACCACCCACTACCTGTATCGGCGGGTGATTGGCCAGTTAAAATTCGTGCGGCGATAGATTTAAGACTACACCAACATAATGACGGGATTTAAGACTAGCTTCATATTTAAGACTAATTTAAGACTAGCTTGTGACCTAGATTTAAGACTAACCATGGTTACTTATTTAAGACTACCAACTTATAATCAGATTTAAGACTACTATCATTTAAGACTAAGATCCCCTCAATTTAAGATTTAAGACTATTTAAGACTAAGACTATGAAATTTAAGACTACGCGGCTATTTAATTTAAGAATTTAAGACTAAAGACTAATTTAAGACTATACAATTTAAGACTAAGACTATAAGACTAGCCAGCCGTAAGTCAATGTGTTGTCCGCTAACACCATGCCATTTAAGACTACTACGGGTCCAGGAAGAGTTCACAAATAGGATTAACGCTCTGCCCACAACCTAAGACCAATAGGTATTTAAGACTACCCGTATTTAAATTTAAGACTATGTATTTTGTTGTGAGCTGTCCCGCTAGATTCATGGGTTGAGTGTCGCGCCTTGTGTGAGCTGTCTGTGAGCTGTCTTATCATGTGAGCTGTCTGTTGTGAGCTGTCAGACTTGTGAGTGTGAGCTGTCATCACGGCAACATGTGAGCTGTCCCCGACCGTTTTGTGAGCTGTCTATTACGTCGCCCAACTATCTGTGAGCTGTGTGAGCTGTCATTGTGAGCTGTCGTGCGTCTGTGTGTGAGCTGTCAATAGTGGCACGCCTGTGAGTGTGAGCTGTCACTAGAATGTGAGCTGTCTCAGATAAATGGAGACCACTACTGTTGTGAGCTGTCCTCGGATGGTGTTAATAATTATAGTGTGAGTGTGAGCTGTCCCCCGTGTGAGCTGTCCTTGTGAGCTGTCTCTAGTGTTGTGAGCTGTCCTTGGTGTGAGCTGTTGTGAGCTGTCGTCTCAGTCCGAGTAGGTTGTCCGTGATGTTGTGAGCTGTCCTGTCGGATGTTAATCCCTCATGGTGTGAGCTGTCGGTACACGATTTGGTCCTAGGTGTGAGCTGTTGTGAGCTGTCGAACCCGATCTGAGCCGAATGAGCAGTTAATCTATAACAGCAGTGCTTTCTTGAGGCTCGCGAAAGAGGTGGTTCGGTGTATGCCTTAATTAGGCCCACCTTGTTTTCCCTGATATTCGACCTCACACGTTCCATGTAAGAGTATGGCCATTGCGCCTAGCCCGTACAACGTTTTATTAGTAATCCGTAACTAGGCCGGTCAGTACTGCAATGCGATTGATGCGTAGGTGCCGTGACCCCGTACGGTGCAAGGAATGGTCATCATGAGATTGGCACGGGCCGGCCTATCGCTCGTTGACGAGCTTTTGGCATCTCACTTGAACGCTAAGCACGCGTATCTGAAGTCGTTGTGCTAACAGGGAGGCCAGGATTGGCCCGCTCACGATGACCCACAAGGGTACAATACCCCGCATTAGTTGACTACTTAACTCTGATCTCTACAGTAGCTTAAGGTTACACGCTAGCAATGGGAGTGTTACGTTTTTTTTCTTTGTTACGCCCTTATCTTGTTAGTCAGGTGCAGCCTTTGCCAGGCCGCATAACATGTCCTCGCCCATTCTCTTAGACCAGGGTCATATTCGCGATTCGATGAGATATCTTTCGACTACTTTAGACACTTCTCGCAAACCCAGGGGTCTCATAGATCCTACTTGGTAACTCTGACTTCGACTGTCAATATTAAGGATGTCCGTCGGTTGAGCCCAACTTCACACGCGATAAATATTTACACAGCTGTTCCGGCGAATTCGGTTCCGGCGATTCCGGCGAATCATTCGGCAGGTGTTCCGGCGTTCCGGCGAATCCCTTCCGGCGAATGATGGTGAGCCTTATTCCCTTCGTTTCCGGCGAATATCGTACGGTATCTACGAGGGATTCCGGCGAATGAATTCCTCTCGAGCTTCTTCATTGTCGTTCCGGCGAATGTTCCGGCGAATGTATTTTCCGGCGAATCTTCCGGCGAATCGGCGAATGCCATTCTTCCGGCGTTCCGGCGAATGCGAATTCCTCGTTATTTCTTCCGGCGAATGTCGATGTGCGGCGTCGGAATTCCGGCGAATGCGACTATGTGCTAGGTTTTCCGGCGAATACGGTAATTCCGTTCCGGCGAATTCCGGCGAATGGCGAATAATAGAACGTTCCGGCGAATAAGTTCCGTTCCGGCGAATGGACCGACTTTCCGGCGAATAGTATAGGGAGTAGGACTGGGAGCCTTGACGTTGTGAGTCCTTTTTGGAAATTCCGGCGAATCCGTGGTTATTCCGGCGAATTCCGGCGAATCCAGCGGAAAGGACTGACGGTTCCGGCGAATCTCGGGCGTCCTTACGGTGCGGTGGATCTGCTCGGACGGTGCAGTTTTGTAGATGTACCTAACCCATTATTAGGAACCGCGAAGAAATCGCAATCGCCGGAAAACGCAGGAAAATACGTACTGTTATTGTCAGATGCGGTTACATGTCTCGGGAATCTTATCTTCTATGACGTTTGAACTCTTCCTGTTTCAATCCACCTAAGCACGGGCATACCGTGGTGCAACCTGATAACCTGTAATCCAAGAGCTGGCTGGGACAGGTTAGTAAGCTTCCTTCTTCTCTGGATCCTCAGAACGTCTCTGCGTCACTTTTTTGCCCCGACGGATGTTCCACTAAGGTCTCGCCCCTGCGGACTAAAATCAGCCTATGGCGTTCCTCGAGCTCGATCCCAGCAAGATCCCGGCTCTGGTTAACCCCCGGCAATTATCCGGGGGGCAAAAGTGACCGACAGAAATGTTTAGGCGGCCTGATGGCGGCGTTGGTGAGATGAGTGATCCTCAAGGTGTAATCACCTGCATTGAGCCTCGCTAACTGTCACGGTTAACCGTAGGGTCTAACGGTCAGTGCAAGCTGTCGGCGGTTTTCAGGAAGGATCTTTTCGGTCCTCAATAGTCCATCTCCCAAAGAATGCCGGGTAAACGATAAAAGGTAGACGAAGTTCGACACGCTGTTGGGAATTTGTATTACAATCGCACTAAAGCAACCAGTACTAGAACCGCTGAAACTCTTTTTCGTTGGGTCCCCCGCTGAAGCTGGTGTCAACTCTTCATACTCAACCGCCTTTACGATTACACCGTTAGGATGCCGGCATCATCTTCACGCGTTAGGATGCCTGCCGGTTAGGATGCCTAGATCACATGTATAGTTAGGATGCCGTGCTCATCAGGTTAGGATGCCTCGTTAGGATGCCCGTTAGGATGCCACCGTTGTTAGGATGCCAGCAGTTAGGATGCCCCCTGCGTGGCACTCGTTGCGGTTAGGATGCCGCCTGTTAGGATGCCTAGTAGGTTAGGTTAGGATGCCGGAGGTTAGGATGCCCTCGTTAGGATGCCGGAAAAGTAGCTGTTAGGATGCCTTTTGTCAAGTAATGTCTCAACAAGCCGGCAGGCACGGGTTAGGATGCCACTTCAACTATATACAAGAGGGACCTTATTGGTCAGGATCTAGTAGCCCTATTATCCCCCGACTAGTAGGTGTTAGGGTAACTAACTCTTAGTAACTAACTCCGTAACTAACTCGTAACTAGTAACTAACTCCCCGAGGGTAACTAACTCCTCCTCTGGTTAGGATGCCCCCTTAAACGTTAGGGTAACGTAACTAACTCCCCGATGGTGTTAGGATGCCGCCTACCGTAACTAACTCGCCGGTAACTAACTCAGTAACTAACTCAGAGCCGCCGAGGAGGGTCTATGATCACGTAACTAACTCAGAGAGTGTAACTAACTCAGTAACTAACTCGTTGGCGTAACTAGTAACTAACTCTGAAAGATCATGAAAGGTAACTAACTCTTATACGTAACTAACTCTGCGAGCCAGAAGGCAGTAACTAACTCCCTCCGAGCGTAACTGTAACTAACTCCTAACTCTAACTAACTCCTCCGTAGGGAGCTCGTAGTAACTAACTCCTCGTAACTAACTCCTCCCTAAGTAACTAACTCTATACCTGCTTTTTCGTAACTAACTCAGGATAGGCGAATGAGCGAGCTCGTAACTAACTCGCTCGTAGTAACTAACTCGCGCAAAGAGCTCGTAGGTAACTATACTAGTCCTTCATCAGGGTGAGCTCGTAGGGTCTCGAATGAGCGAGCTCGAGCTCGTAGGCGAATTTCGCGAGCTCGTAGGTAGGAGACCGTTGACTTTCCAGGAGCTCGTAGAGCTCGTAGGAGCTCGTAGGGAGACCTTTTGCTGACAAAGGAGAGCTCGTAGGGGAGCTCGTAGGACCCGGAGCTCGTAGGTCCTGAGCTCGTAGGTTGAGACAAAACCGACCAGGCTTGCGAAAGGTTAAGAGCGGCTTGCGAAAGCGCTTGCGAAAGTTGCGAAAGCGTAGGTTAATCTTAGCATGCGTTCATCCGAGCTCGCTTGCGAAAGGAGCTCGCTTGCGGCTTGCGAAAGGAAGCTTGCGAAAGTCGAGCCCGTGGTGCTTGCGAAAGTGGAAGTGCTACAGCTTGCGAAAGAGCCCATTGTCCCGTCGCTTGCGAAAGAGCTTGCGAAAGCTTGCGAAAGCGGAACGTGCTTGCGAAAGTAATGCAGAACCACGGCTTGCGAAAGGCGAGGCGCTGCTTGCGAAAGGGGATGGACGCTTGCGAAAGCTGCTTGCGAAAGGAAAGGACGCGATTGATGTTCTTATAGCGATGATCGTGCGTAAGAGGGCTTGCGAAAGTGGCTTGCGAAAGATATGCTTGCGAAAGACGGCTTGCGCTTGCGAAAGTTGCGAAAGGCTTGCGAAAGATCGTGGGCTTGCGAAAGCGTCAAATGCTTGCTTGCGAAAGAGGGGATCCCGCGCGGCTTGCGAAAGCGCTTGCGAAAGCGCTTGCGAAAGGGCTTGCGAAAGGGATATTGTTTATCGGGACACCCCGGTGACCGATGATCGAAATCCGCCCGTAAGTAGATCTAAATAATAAAGTGAGGATATCGAATATCTTAGCTTTGTCCCACATATTCACGACAGTAGGGAAGTGGGGATTCTAGACAGCATACGTCGGTTACGTTCTTTACATGTAAGAGGCTGTAATCCGACGAAACCTAACCCACAAAACGTGTCGAAGTTTAGAGCACTGCGGATGGAGACGTCATCACTAATTAACAGACCCCTGTATGCCGGTGCGGATACTTTTAACACAGACACTATGGCCCACAACAATCGAGATGGAGACGCTCCTTCCAGAGTCTGTGACAACGACGGAAGCTACGAGAGCCTGGGTCCGCTAGTTGATCGATGTGGTGGGTCGGAGATGCGCCGTTCAGATCACCAAGGAGATGGGGCTCGACATGAGTTGTGTAGGAGGGCGATAGATCACTTGCGTTGGGTCATGCGTATTTGGCCCCAATCAGAAAATCGAATCACGCTATCACATTTGCGTTGAAAAATTTCCAGTTCTTGTACTAAAACGGTTTCCGGAGTCCACTTATCGTTTAGTGGCTCACCGAAGAGCGCGCCGTATTGTCTATTATGTCAAATCCCCAACACCTACTCCTGGCGGAAAGGTAGTCGATTCCAACCCGCTTCTATATCGCCAGTTAGCTGGGGTGCTGTGACCACCTACTGAATTAGATGCTTTCGGCATCAACTGTCCCATTTCGAACCCTACTCAGCCGCGTTAGGTAGTGATATGTATAGGCAACAACTCCGTCTGTCTTCCACAATGTGATACCCCCACTAACCCTCGGGGGTTGTGGAGAGAAGGCGTGCCTCTCACTAGGTAATTGAGAAGAGGCTCTGCCCTATGCAAGCGCCCGGGTTTAGGGAGTCCTTGCCTACTACGAGACTCCCTTGCCTTATGTGTCTTTTGATCGCAGATTCTTGAACTATTAGGCCATTTAGCAGAATTACTCCTTTAGGTGACTCCGGCAACTGGTCGTTTGTGCAATTATACATTGATTCACTCATTTGATTCACTCTATTGATTCACTCTCTGATTCACTCTGCGACCACGCTTCGCTTGATTCATTGATTCACTCTTGATTCACTCTTCATTGATTCGCTCAAGTAAGTACCTTCTTTAGTACCTTCTTACCTTCTTCTCATGAGTGGAAGTACCTTCTTATTCACTCAGATTAACCTCTAAGAGTACCTTCTTTTGTGGGAGTAGCTTGTACCTTGGGAGTAGCTGGGAGTAGCTACGCTCAAGTACCTTCTTGAGTACTGGGAGTAGCTTAGCTAAAGTACCTAGTACCAGTGGGAGTAGCTACCTTCTTCTCATGGGAGTAGCTGCTTGGGAGTAGCTAGTGGGATGGGAGTAGCTCCTTCTTTTGGGAGTAGCTATGGGAGTAGCTAAGTACCATGGGAGTAGCTATGAGTGCTGGGAGTGGGAGTAGCTGGGAGTAGCTACAGTACCTTCTGGGAGTAGCTCCCGCTCATGAGTGTTCTGGGAGTAGCTGGGATGATGTGGGAGTAGCTTTCTTATGGGAGTAGCTCTTTTCATGAGTGCTCATGAGTGGGTGGGAGTAGCTTATCTCGACCTGGTGGGATTGGTGGGAGTAGCTGGTGACCTGGTATCACCTGGTATCTGTGGGAGTAGCTTGGTATCTCAGTACCTTCTTCTGGGAGTAGCTCTTCTTTGGGAGTAGCTCTGGGAGTTGGGAGTAGCTGGGACCTGGTATCCCCAATGGGAGTTGGGAGTAGCTATCCTGGTATCCCCCAAGACCTGGTATCAATGGCACCGGTCGACCTGGTATCCAACCCAACCACGACCTGGTATCCTGGTATCACCGGTCCCAACCACGAAAGGGACCTGGGACGACCTGGTATCTAATCCCCAACCCAACCACGAAAGGCCCAACCACGAAGCACCGGTCGACCTGACCTGGTATCGGTATCTAAAAAAGCGACCTGGTATCCATGACCTGGTATCACGCACCGGTCAGACCTGGTATCCTAATGCTGTTAAAAAGGCGACCTGGTATCCAAATTAAAAGCCATGCGACCTGGTATCCGAGACCTGGTATCAAAGCCATGCACTCATCGACCGACCTGGTATCAAAGCAAAAGCCATGCGCCAAAAAGCCATGCGACGACCTGGTATCTCTCTGACCGGAGAAAAGCCATGCGGGACATGCCATGAAGCGTAAAAGCCATGCCTACAAAAGCCATGCCCAAAAGCCATGCAAAAAAAGCCATGCAACCCCAAAAGCCATGCGAAAGAAAAGCCATGCATCGGTGTCCATCGCTTAAGTTAAAAGCCATGCAAAAAGCCATGCAAAAGCCATGCCCCATGCCTAAAAAGCCATGCTCTGAATAGTGAGTATGCAGGGAAATTTGTTGTCACTATGAATAGTGAGTGTGGACGGTGCTGAATAGTGAGGTGAGTGCGCCAAGGCACGGAGTGCTGAATATGAATAGTGAGTAGGGCTGACCTGAATAGTGAGTTCCGTTCAGCGTTTAGTGAATAGTGAGCAAAGTATGAATAGTGAGCTGCTTGAATAGTGAGTGGACAAGAGATATTGAATAGTGAGGCACAGTCCGACTGTCTGTTTGATAATGAATAGTGAGATGAATAGTGAGGAGGAATATCCCCTGGCGGTTCTGCGTGAATAGTGAGAATCGACATTCAACCCAATGGCTCCCGCGCTAGCCGTGAATAGTGAGTGAGCACTCATGAATGAATAGTGAATAGTGAGAAATACGATCCATCAATTGATGAATAGTGAGGAGTGGCTCAAATGAATAGTGAGTAGTGAGTGGATTATGTGAATAGTGAGTAGTGAGATGAATAGTGAGTCGACACGCTAAGGTGAATAGTGAGTATGAATGAATAGTGAGTGGCGGTGAAATGAATAGTGAGGTGAGACCGCGCATTGAATAGTGTGAATAGTGAGAGTGGGTCTTGCGGACTAATATTTTTGACGGGTTGCTAGCTTGTCATCTGACTGCTATGCTAGCTTGTTAATACTTGGATCCTACGATATTGCTAGCTTGTTGAAATGCTATGCTAGCTTGTGTGCTAGCTTGTTGCTAGCTTGTCTGTAGCTGCTAGCTTGTGCTAGCTTGTTGCTAGCTTGTATGCTAGTGCGACAGGTAAAAGGAATGTAAGGACAACAAGATGCACATGTTGTGCTAGCTTGTCTTTGCTAGCTTGTCTTGTACTAGAGAACGTGCTAGCTTGTTGCTAGCTTGTTGCTAGCTTGTGGGACGTCCTCATTGCTAGCTTGTTACATGCTATGCTAGCTTGTGTTTTATGTCAGGTAGGGCTGTTACCCCCCGTCCAATTCAAGTGCTAGCTTGTACCATAAAGCAATTTGCTAGCTTGTCGTTTAGCATGGGATGTGCGTGTGCTAGCTTGTCTAGCTTGTTGCTAGCTTGTTGCTAGCTTTGCTAGCTTGTCTAGCTTGTTGCTAGCTTGTTGCTAGCTTGTTGCTAGCTTGTTGCTAGCTTGTGAATCTGCACAGAATCTGCACAGAATCTGCACAGAATCTGCACAGAATCTGCACAGAATCTGCACAGAATCTGCACAGAATCTGCACAGAATCTGCACAGAATCTGCACAGAATCTGCACAGAATCTGCACAGAATCTGCACAGAATCTGCACAGAATCTGCACAGAATCTGCACAGAATCTGCACAGAATCTGCACAGAATCTGCACAGAATCTGCACAGAATCTGCACAGAATCTGCACAGAATCTGCACAGAATCTGCACAGAATCTGCACAGAATCTGCACAGAATCTGCACAGAATCTGCACAGAATCTGCACAGAATCTGCACAGAATCTGCACAGAATCTGCACAGAATCTGCACAGAATCTGCACAGAATCTGCACAGAATCTGCACAGAATCTGCACAGAATCTGCACAGAATCTGCACAGAATCTGCACA"
    parameters = "11 557 20"
    k, L, t = [int(i) for i in parameters.split(' ')]
    kmers_lst = kmers(genome, k)
    indices = {i: kmp(genome, i) for i in kmers_lst}
    potential_clumps = {i: indices[i] for i in indices if len(indices[i]) >= t}
    clumps = {i: get_clumps(potential_clumps[i], L, t) for i in potential_clumps}
    clumps = ' '.join({i for i, j in clumps.items() if j is True})
    print(clumps)
    exit()

