def hamming_distance(str1, str2):
    return sum(1 for i, j in zip(str1, str2) if i != j)


def approximate_pattern_search(text, pattern, mismatch):
    distances = [(hamming_distance(text[i:i+len(pattern)], pattern), i)
                 for i in range(len(text)-len(pattern)) if hamming_distance(text[i:i+len(pattern)], pattern)<=mismatch]
    return [i[1] for i in distances]


if __name__ == '__main__':
    genome = "AAGGCCACGTGAGATGCAGTTGACCAAGGCTGCTGGACTCATTCTCATAGACAGCTCGGGAATCTTTCGTACATGATGCAGCACATTTGACCACTCGGCCTCCTGACGCGGAGGTGCTACTCCAAACACAGAAAGGGGCAGATAGTGAATGAGCGTTTATCATGCTTGCCTGGCGACGATGCCATTGTACTAGAACATAGGCATGAAACATGAGAGGATTTGAGAGTTGGTTTAATCGATCACTGGACTACCATGCTGAATGAGCGCACAAGTCGGAAAGCGTTGTGATAGGCGTAGCATGTGTCATACATTCGTGATATACGGAGCTCCTGGTATTTGGCGGCGCTAAACGTCTACTCCCGCCCTACTGGCCAGATAGATGTTCGCTTACCTGAGGACTCGGGGAGTTTCAGGCGTGCAAGAAAGTCGAACATGCGGGCGTGGCCCGAACTTGTACCTGTTCCGTCATCGTTTACCCCAGTTAGAGTGGACCCTGTGCATATCGACCGAGTCTAAGTCGAGTTCCGGCAGGTGAAGGGGTCCACTAACAAAGACTACTCGAGAATACCATTTGTGTCCTGTAGTTACAGGTCGAGACGCCCACTACAATGGGACGGTACTTGGGCGGCTCAGGCGAAGTGATGTCTAAATGTTTTTCCATCAGTAGGAGTTTCCCTTGATTAAGTGACCTCTCGACGAGACAACCTACCCAAAATGTGCCCGGCCCAACAAGTGTTTCATGGCTTACGTATAACTGACTGTCTCTATGTCTGCGAACCCTGGGGAACTTTCAGGTCACTAATACAGCGAAGCTTTCCATGAACGCAGTCCCTTCTAAAATGTCGCAGTTCGAGGGACTATTGGGATTGCGGTCGACGTAGAATTAACAAGGAGAGCTGTACTCGAATGGTCCCCCCCCGGAGAGAGTCTAGTGACGCTTGCAGACCAGGCAGCTGTTTTAATTATAATCTCTCCGCCTCATAAGGCCAATTATGGGGCTCGGGTGGTCGTCAAGTGACAACAAACTAATTAATGGTCATAATCTAGCACATTCCCGTGTGGAAACAATGCAAGTCGTGTGTTGATGTGCGTATTGGGCAGCGAACCCAAGCAGGGCATCTACTGGATACAGAGGTGATCGCAACACATGCCTCGCGACTCATTATGTGCCTACACGCTAACCTTTGTACGGTGCATAAGGTCTAGCGAGCCATCTCACCTGCTAGTTACGCTACATCAACGCGGAACCAATTCGCAAGGGTATCGACCGAAGTGGCCACGCTGAAGAGGGTTTATCTCAGAGTCTTAGGACGGAGAAAATGGCACCACAAGCGATTCTATGCGCGGACACGTCAGCTGAACTGGGGACGGGGAGGCCCAGATCAGCCTCTCATCGCGTCAAGACCTGTGAACGCTCTGGTCACTATTAGTCAGTGTATGTAAATGCCGGGCCTCAGGCGCCGGATTCCTGGGGAGGTGTTCCTGTCAGTAGCAACCGGGTACAACCGCGGTAAGGCTGCTCTACAGTGAATTGAATCTCCGCTCCCGGCCTCAGCCCTATCACAGCTATGAGTAGCCGCCCCCCCCATGGGCAGCCTCTCCGGTTTAAGGTCCCAAGGCGTCAGGACAAAACCGGAGCTGTACACTACCAGCGAGGCATATAGCACAGGGTGTGCTCGACGCTGTAAGAGGAAATGACCTCGGGCAGAAAAAGTTTAGCGTTGTGCCGATTTGCCGGTTCTGGTACCAAGTATCCAAAAGTCTTCCCGGAGTGCCCTGAAAGGTGGTTGACGGAACCTCTTAAAAAGATACGTTTCAAGGAGTCGGGCGTTACTGGCCCAAATCCATCTCGACTTGGAGTCGTATCATCCTTCTATTCGTATAACCTACATGCGTAATTCGGACCAGCTTTCTCTAGGAGTGTTGCAAAGTTCTCTCATCTCGGACGCCCGTTTACGATAGCGACGCAGGCGCGCGACAGGGGGTGAGTTAACTAGCCTAATTGTATTGTCGTACCAAGGTCCCGCGTGGCGATGAGGCCTTTGACGTTCGGCTAAGACGGACAGTGGCAAGGTTGTCACTTTAGTCTTTCCCACGTATGGTTTCAACGTGGATAGCCGTACAACGCCCATCCGATATACGCCTGCGATTTATCTAGCTATAGTTAGCAGTGGGGATCCACTTTAGCATGAAAGGTATCGTCGAGACACTAGATAGGACTCGGGTCCTTTTGCGTCGATAAAACCCTCTGCATAATTTAGGCCAACAAGTATGCTTGATGGAACGCACTGGCAACACTCACTATGGACCTAACCATCTCCGTTTTCTGTCTGCATCCGTCTGCCCCTGTATGCAGTTGGCCCAGATACGGCTCCAGGTTAACTTTCTCACATTGGCCAAATCGGACCTAGAAGCACCCCCAGCCATGCAAATCTGGTTTAACGCCTAGCAAGACAACCTGTAAGTGACTGGTGGTGTTGAGCCTACCCACACGCCCTATAACGCTTAACTTGCAGTAGCGTAAGCTGTAGACTGTCTTAGCCTCTAGCATGCCGGCTTTTGTGCGGCGAGACATTGAAACAATCAAGTTGTCATGAGTACCCTCGACGCCCCGATGAGATGAGCTAGCCGCGATATCGGCCTTGACACCCTTATTATACTACATAGATTTCCGGTTGGGCCGGAAGTACACCTTTTGAACCGCTGCCCAGAGATGTCCACAGCGGGCACACAGCTTGGATATACTGGCTCTACTGAAACTGCCTGCATCAGGAGGGCAAACTGATCGTTAGTCAAGCCGCCAGGAGTCTTCGCGATACTCCCCAGCCGACCGGTATGGGCTTTAACTTTGCCACGTGTTTGTGTCGAGAGGGGAACCTTTTCTGCCGTCAACGGAGTCTGTCGTGGACCAACCCTAAGTATATGTTAGCCCATACCATTGCGGACTTCCGTTGACTGAAGAGACCGCTCAGTTAGGATTCTGTCACCAGCTCATGGAGCGCAACGAATCCTTAGGCGGAGACACCGTACGCTAGTCCGTGTAATGGAGGGCGTCTACCAAGATCCCGGACCACATAAAGTAACACACTAAGGATACAATATCAGTTCGGGCCAGGTGATTAGGGTTAAGTCTTTTATGTGATCGAACCCACATGTTACACCTTGATTGCCCATGTCGGGCGAATGGTGCACTGCATAACGACTCCTTTTACTAGCCCTATAACCCAAGGAGATTACTGCAGTGGCAAAAAAGAAAACTGCCATGCACTGGGCTCAATTTGCCTAATTATGCCTTATCAGCTTTTCCCATTTGGATTTGTCAGGTTCATCACGAAGGCTCCTAAACTGGGTGACGAGCGGAATCGGTTACGCTAGTTAGCCATTCGGAAATCGGGGCTGACATGATCCACGGGTGTGCCAATTGAGAATAGCAAGGGGAGGTGCTCCCCTGGTGTACGTGCCGGCAACATGACCGTAAGTCTTCTCGTACCTCTCGACAGATCACCGCCGTGGGAAACAGAAGTTCCAACTTCTCAACTGCGAAATTACAATCCCTCTATAGTACTCTTACATGTGGGGTGGTAACTATAGCTCGCGATTTAACTATACCAGGATTAAATTATTAATCCTCAAGATAACAGATAACTGTCCGTCCGGAGACGAGCGTAGTAAAGAGCCATGAGGATCAGGTTTGGCGCAACATCCCGGACAACCAAAGCATCTACAGGGTTTAACATTTCCGCGCTTCTATTTTCCACAAAGAAAAAAGACTGCACTCGACGATACAAAGGCGTTGCCTGTGAGATCCCCACGGCCACTTCTCCATGGATCCCCCGGTTCGGAATTATTGGGCCGCCATTCTGTGCGCGAGTTACCATACTGGAATGCCGCAAGCGCCACTCCTCTCGAACCCAACCAATAACTGAACGCGCATCTAATCTTGTAAAACCGTGATGCGGTAATCCGGATGAAATCCTACCAGGCTCTTTTTAAGGTGTGATTTGCCTAGATCTTTCTGGCCGCCTTCTGATTCCAGCTTCGCTCTTTTCCAGGGTCCCGGAACTTGCTAACCTTTTTACGTTGCAAGGCTGCGAATGGATTAGTTAGCATGAAAACTAGACTGGCGGGGTCTTCCCGCACCAATGCATGGGGGGTCGTGGTTTAGAGTCCGGGAGAACCAAGAACACCGGCCATGAAGATGCCCGTCCGGTTGGAAAGGAGAATGGTTCATTAAACTTTTATTCAAATAGAGCAAAATCCTGTGATGGCCACCTCGGCCAGTATACACCTGACTGCCACGGACGACCTACCAGACTTAACATCTCACGTTGGGCTCACAAGCCGTTGATCGTCATGGCGGGCACACTGGTGTGGGGCAACGTCTGAGTTTAGTAGTTCAAGCCAATAGCGGCACCTGCAATAAGTGAATGCTATCAATTCCCAGCAGTCCTCGGTGCAAGACAATTGGTCATTTCGCGGGTCTGCGCGGTTTAATAGCACTGATGTGAGGTGCCCCTATCAGGCATTCCGTCCTATCCAAAGCAGGCTCCATGACAGGAGTGCAGTCAAAACAAAGCTCCAGGGAGATAGAAACTACGGCGTGTGGCTCGACTATATCGGTCCCGGAACCAGTACCCTGCATCCTAAAACAGGACCCTTATTGACGAAGCTGTTGAGCTCTGTGGGGTTCTACAGTAGATTAGTCATCAACGATTGGCTATCCATACCATTGAGGGGGTGGAGCTCTACCAGATTTTGACGCTACCGTCGCAAAACCCTAAGGATAATGAACATCCAATCGGTTTCATGGAGCGAAAGTCCTCATCGGGGAATAGGTGTCCTGGCGCATCTATACCGGTGTTACGGTCATCTGCTGAAAGAGACGTTGCTTGTTCGCTAAGAGTCATCAGAAGAGCCGCAATGGGCGCTAGCAAATCTTGGTTTTCAGCCCTCTATCATTTTTGCCGCGGGCACGTCACGAAAAACCGATATCTACAGGTTATACCTGGTGAGTCAAACTGGCATATGTGACCTTATCTAGCAGTCTGAGGCAAACTCTACACTCCCCGGTAGTTGGGTAGAAGAATTGCAATTCAGGGTGCCGAATGTGGGTAGGTGACCGACCGTGGTGCAATCTAATCTGGCACTTATTTGTACTTACTTAGGCGAACAATCGGGCGTTGGACTTTTCAGGTTCTGGCTTGCCTATACGTGCATTCCGGTTTATGCCTAAACTCGTAAGCGATCTAGCACGCCCGGGTCCCCTGACAATCTCCGTGTTGGGGAACCCCACCAATAATACTGAGATTCGGGCGAACCCATCTAGAGCTAGCACCTTATGTTGAAAATATACAGGGGGCTGCGCGCCGGACGGGAATTTGCCCTCCCCCTTATTGCTACTAGGTCCTCAGCTGAATGAGGCGAAAGCGATCACTCCATCAACCCCATCTACATGGCGCGCTAGAACCTCCAGCATGGGCGTCAATGGGACCGTTTCAAGCGCAATGTACATGAGCCACCAAAGTCAGCCCGGACCGGTACCCTAGCCGATCAAGCAATAGCTTTTAGGGTTACAATGGGTTACATTCAAATGCGGTGATTACCAGCTCCCACAGGGGGATAACCGTCGCTCCGGGCATAACGATCGGGATTCAGAGCGCATCGGTGTTTGTATGTGTGGTGGAGCGTCTAAACCCCAAAACTTGACAAAAACCGACGACTTTTTTTGTTCCCCACCTCATGTAAATTTACGCTTCCGCTACAAATGTGGTCTGTAGTCGTGATAATTCGTCTATAGTGACGATATTACGAACTTGCCTAGCCTTTATTATGCTCTCGCGAACTCATCTGAACCATAGTATTTTGGGGCCCGATTTGCTCTAGATACGACAGATCGATAAGCTCTCGTGTTTATTGAGTTAAGAATCGTTATTCGCAACCGGAATAACCCCTGCTGTTAGTGAGAATTGTACAAATAGCGGCCCGCATGGTAGACAACCGTCGAGGCTCTGCACGAGAATCTACGGCTCGGGCTTCCTCGTTGGCACATCTGCGCTTTATAAAGTGTCGTATTCGAACTACCGTAAGCTCCGATCCCACATACCCAATCGCGGAGAAGACAACTAACGTGGTGCTACCGTGCTGAATTCGTTTGGTTATCAACGTGCTCTGCCAACGTAGGATCCATACGAAGTGGGGTTTAATGAGGGGGGTGGTCTCAAAACCTGGGGCATAGCCAAGTTAATGATTCGTACCCCAGCCCAGAGTATAAACGTAGTAAATTTTCGATCTTTCACTGTAGAGGGTGATGTATCACTTAACGTCACTATAGAGTTGCGGGTGAAGAACAGGCGTACCTGGACAGATAGGCTCAGCTCACCCCCGACATGCGCTTATGACGATGCCAGCCAGCATCGAGCCAACATGAGGCACACATGAGTTCTCCGTATCCGTGACTGGATTCATCGCGGGTCACGGCAATAACCTTTGAAGACAGTTAAGCAATGGGATTAGAACTAAACGTAACTATCACGATTTATTGCGTCGGAATCTGTGAGTTAGAAAAAGGCCCTCCGCGATTGTCCTTTTGGTCATCTCCGATTGCGCCATGTGACTGTTTTAGACAGCCTGGGGCATAAATTATGAGCGACGCTACCACAGCTCACGATCCCGGCGTAAGGGTTTCAGAAATTGTGGGAGACCACGGTTGGCAGTGGGCCGGGACGATGCTAATCCATGTCAATTAAGTTTGGATCGGCATGGCAAGTTTATGTTTAGCAACATCTGGACTCCAGGCACGGAGCGGAGTGGGACACCTATGGGCAGCGGCCTAGGCCCCGCCCGCTCTGGAACCAAGAGTGGCGGAAAACCAGATGCGAACTTTACGTAAACTTTGCACGCTTTACGATGCAGTTCTGGAGGTTAGAACCTGCGATACAGGCACCTACGGTTGTTTGGTAATCTCAACTACACTGAGGTCGTGGGCCTTTCCATATTTTCCGGGTGCTCTGATAGACTTTCCTACCTTGAGACCGGACGCCCGAAGCGAGACACCATGGCAACCTTGGTTACCGTGGACCACCCACCGTAGCGGGCAAGTCGAGAGTAACGTCGGGGGCATATATGCGGGGGCCACCCTCCGCCCAAGCTAGCTCCTTCTCGAGGTCGGGACTTCAGTTACCAAGATATCATTAGTATTACTTATTCTCGCAGCAATACGAAAACTCATAGTGATGTCTTGCATACACCTTGTTGGACAGATTAATGGAGTGAACTGGTGATTATCAATGGAGCGGGATGCCGCGCTCCATACTTGAGGTGCGTCTCTCAATTAAGCACGCTATGTCCCATCGCAGATGGTACCCGCTAGACATAAGATTGGTTAGGCTAAGGTGGAGCCTGAGCCATATGGCAAGAGTCATCGTGTTAATGTAATGTTCATGCTGACTTATGAGGAACCTCGCTCCGATCACCGGTATGAATGACCCTTTGCATATAATGCTAGAGTGGCGGCACCCCTTATAAATTATGACTAATCAGGATCGTCGGACCAGACATTACCTCTTTAGTAGTGTTGCACGATGCGAACTGTACGACATTGGAGCATATATAATAGCGTGACACGTCGATCTGATATCCCGAGTCCAGAATACTCCTCCAGAGATCCGTCCCTCTTGGACCTTTCGCGCGCAACAACTAAGAGTGGGATGGAGCACAACTGATCGCCTCGATGATCTATTTTTGATCCAACTCGTACGATGGGAGTAGCGTTTGTGAGACCATTCCCACCACGAATTATGGGTGTGAGCTTGGTGATCGCTACGAAGAGCGACATTGGCAGTATCCGGTGAATCTACATGTGTATGAGTGCGATACCCATGCAGTAGTTCGACCGTGCCGGACGCCCGGTCCTAATCCGTTCAGCGACATAGTTTTGGTATGATGGCCTGAGGTCAGCCGCGAACACACGCAAGAGCTACGCCCTGATAGGCAAGGGATGTATTGAAATTCCATCCTCCGGATCTTGTCCCGTGGAGAATTTTGGGAGAACATTACTTCTCCTCGTTCGCCATAAATCCAACTTCCGCCAGCGCTGCTGCCTTATGAAGATCTAGTGAGGGGTTAGGTTGCTGTGTTGCAAGACAGGCACCTTACATCTACAAGTCACCGGGAAATGCAATAAGGTAAAATGAGCGGACAGAAAAACTGCGGCGTAAATCGCAAAGTGAAATCATTACTAAAACTTATCGGCAGTTGGTGCGAATGATTGAACAGAGGAGTTTCCTGTCTGGGTGAGCACGCAGATCGCTGGTAGGGTGGAAAGTCAATGATCTCCGAATGCTGCATCGGTTTGCGACCGAGTCGTTGTTGGCCACATAGCGAGTTAACGTTTTTGCAAATGGTGCTGCGTAGGTTGATCGATGTTACGGACAAGACTAGTCGATAAAACGACGACGGCGAGGTTCATTTTCATGCTTACCACGTAGGAATGGAAGTGAAGGCACAGCTTCTCGTGACTAATGTTTAAGTAGACTTAGGGTAAATGACTAACCCTTGTGAATTTACTTAATCACGAGATAAGATCGTTCAAATGGGTATTATGGCCCGCAATCTCTTGCATTACGCGACGCCGCAACTCTACAAAACAGCAAACACGCCCTCCTTTGTGCTCCTGAAGCGCGCAGCTGGTGAGACCTCTATATAAAGCCAACGCTCACTGCGTAGCCTACCGGGCTATATCCAAGCTAGTCGTTCAGATTGGCCGTCCCTTTTTAGGAATAGAAACGGGACCATTCCGTGCGCAAGTAGACTAAATTCATTACGTTACCCGCTGGGCATCGTTTATAGACACATTAAGCGCAATATCCGTCTTCAAGCGCACCGGGGCCTTCCTCCAATGAATCACGCCTGTCGGTACATTCTTGGATTATGAGTCATCGGAACGCCACAAGAGCAGTGAGCTAAATATTGGCGAGCTATTATGGTACAGTCGCTTGCCTTCCGTTCTTTAACGCGTGGTGGACTTGAACGAAAGAGGGCAGCTCTGCGAGTTGGAGAGGGGAGTCTTTTCACAATGTAGTTAATGCTTCTGAAGGCCAAGTGCTCTGCGCTGGGCGATTACTAACCGGTTAGGAGACACTACTTTACGTACTAGGTCGTGTAGGGAAAAATGGTGCGCAACTCGAGCACATAGGACGCGGGATTTGACGTTCAAACATCACCGACGCCTCAATACGTTGATCACGGGCCTCCCAGCCTGTTTTATGCACCATACGGCGGTGTTCGAACTGGTGTCTTGGCGGCTCTATCTATAGACTGACGACTGTGCATGGTTTTAATAAAACGTCTCCATCGCCATCATAGTGCACAGGTCAACAAATTGCCTCACCTAGGAGCTGCCGGAAGGGTCTTATGCTGAACTCAGAACTGATATAGTCGACGTCGTAGTTCGGCCTCACCGCCTAAGCATCAGGCTTCGATTTAAGCTCAAGTTCGTAATCCACGACGGTGCCTCCCATTGGTGTGTTTCGTTCGTGTTCGGCCGCAAGCCACAGGTGCAACCTGAAGTGAGGTAAGCTAGGTGCTCTTTTCTTGCTTATTTCGAGCTAATCCACTCTAATGCAGGGGTACACTTTAATGTTGATTCAGAAAAATGTCTAGATTCCTCTCCAGTTGCATAGCTCGTAATGGATGACGCGCGTCGCTTTGCTGCGCACAGTCATGCTGACTTCGTGAAGGGCAACACGAGTGGCAGTCGATCACATGCATCGGCGGGCATCATTGACGCCCCCCTTAGCCCGAGGATGTTATCAGTACAATTACAAGTGGAGGGCGACCAAGCAGGCAAAAGGCTTCGATGGATCACACCCCGCTGCGTGTATGTTACTTCCACATTTAGATGGGCCTGGAACTCCCCTCGTCGACATGACATTATGTTGCGCTCCATCTCCCGCTACCAACGGTAGGTTTGCGATTGACTTCAACTGAAGTTACAAGCCTAGATGGTCGCTAGAAATTGGCATTTGCTAAACACTCTCGGCCCGTTATGCCTGCGCTTAGGGAAGGGTCAAGGCGATTGAACGTGCTGATTATTGACGAAGGTTGTTCAGGAGATGATATTCAAATTGACTACATGGGAACGACTCCTGGAGGTCACTTTTTCCCAAGTGGTTAGCTCATACTGAGGAATCAGCTCGAAAGCTTAGAAACATGCGAACGGGAGGCCTCTGAACTCGAGAGCTCGGAGCCTAACAGACACCAGACCCCTAAGTTGCCTTATACGGCCGCGATGAGAGACCCGACAGTAAGGCAGCAAAAAGGTTTAGTTCGATCGCACAGGGCTTAAGTCCGCCATTCTAGGGTTGGTCATATTAGTGTTCCATACCACGTCTAAGGGGGGGCGATTAGCCTCAATCTCAAAGGCTAATACCGGATAGCCAGCCGTACACACTCAACTAAGGTCTAGTGGCACTATGCTCCTAGTGATCCGTAGATAGGTCTTAGAAGCTTAAGTCAGATATCCGAATCAGGTCAACCCAACCCGCCACGGACCGAACCTAGTGATAATAAAGTCCTAGCAGGTTGGGGAACGGTCATGACGGTCGTGCCTTAAGATTCCGTACAGTCGCGGCATATACTTCCGTATAAGGGGGGTGGTTATAATCATAAAAGGAGATGAATGCTCCGCGTAGTCTAACCCTTCTGTAACCGGTACATTTGACTGTTAAGATCTGGCCATTAATGAACGGTCTCATCTCCCTCACCAGGAGGTCCGGGGGCACCCCTTATCCCCTCCTACTCCATACCCGCGAGTACTCACCTAGTGTCCCATCGTCCTGATGTTAGCGTCGAGCGGAGCTGTACAAGAGTCGACGTTGGACCGCTGAACATATGTGAGAGAGCGCTTGGAAGTGGTTTGCGTAGCCTTTTGGTCCGGACCTCATGGTTAGTGGCTTGGGCAGAGATTCCCGCGGATCTATTATCGTCTCGGTCGCTAGACTAAAGCTACCTATATGCACGTGGAACTTTGCAAAAGAACCGCGTTGAGAGCAATCGTTCGCCTGACCACTGACTCTCTAGGGAACCCGGCTTGGGTGAAGTAGTACTGCCCACTAAATAGATCCCACTGGCACCTGCTCTTCTCAAGTCCAAACTCCAGCCCATACCAAAAGGTGCCGGAGTTTAGTAACACACATGTTAATTCCTGGTGACGGCTCCGGTTTAAATGGCTCAACAGCTCCTCACTAAGGCCAGACCTAATTGACGTATGGAAGAATTCTTTAAGTAGCTGGGGCCACTGAAACACGCTTCTCCCAAGATGAGTCGAGCGGACCCGTCATTTAAATCCGCTGTTGGCCAGTCCAGAATAGGGTTACATCTTGTACAGAATGCCACCGACTGTCTCGGCCCCTAGAACTTCCGGTGGCTCGCTCACACATAGAGGCACCTCACCGGAATGAGGTAGATCGTGGACGGGACCCAACTGACGGAGGCCGAGTCATGACGCCGGAAGTGTACAAACAACATAACACGGGCGTCTTGCGGTCCGTAAGTTGGGTCCTCGACAATGAGTGGAGTAACAACAGAAAGTACGCCGGTGTGCATAATCAAGCTAAAAAAGGAACTGGAAATATTTGATCAAATCAAAAAGGTGCTTAAATTGAGCGATTATAAACACATATTGTCTGCAGGATTCGGGTCTCACCGATCTAACGCCGAGAACTGTTGGCAACCCCCCCGAATACAGTGAGGCCCTGATCCAGAATACGTGGGATTTTCTGTCTAAGCCAGTTACGACACTTGCTAACAAAACTGCGGACCGGGACCTGTGGCCGCAGTTGGGTCTTCAGTCAAGCACACGGTGAGCTTGTCGGCCGTATGGAGACGTACTAATGTGACTCTTCATAGGCGCTGGAGATTGAGTCAATGAGTTCCTACCCGCGCCCTGTCCCGGAAAAAGGAATTATTAAGTCTGAAAATATATTTTATTCTGCACAATAACCTGTTCGTCCCGGGCCTAGACCAGGTAGGGAACGCCCGGACGAAGAACTATTCTAGACGCCAGTACGTAGCAGATTTGCATACAGAGTATCTCTAGGGCTTAAAACTTGTCGTTGAAGCATGTGCTCACATTTCAAATCGATCTTATGATGAACCGTAGGTACAACCTGACGATATAGACGGGAGTGGAACCCTAGTAAGTCGGAATACTCTGTAACCACTAAATCAAAATTGCCCGTCTGGGCACTCAACCACTAATGCGTGCTTAACTTGCGTATTTGAATTACAGGCTTCCCCGTGAGTCTTAACTCGAGAGATCGTCATTTTGGCTAGAACGATCGGTCTTATAATGTCAAGCAGAGAAGCGACCTTGCAAAGCGTGACATCGCCAACCTGTGACTGCGTCCCATAAGCACACGGGTCGCTGGTATTGAGAGACCAGGATCCCCCGGCCGAGTAAAGACGGCACTATTACAGACTAGTATAGTTATCGATGCAAAACGTGTGGTCAGAGGTTCTTTGTTGCATCAGTGGTCAGGCCTTCAAGGGCTTAAAATAGGGGAAAGTCTGTTTTATGGTTCTTCTTATGACGAAGCCTCTAAGGGTATTGTTAAGAGCCAGATGGAGCTCACGAACCTCAGACGAGGATAAGTAACTGACCCTATTCACACACAGGCTTGATTCTGTGAGAATGTCCGTACGAGTGATACGCAGGGTTTCGGCACGCACGGGACGATTACAGGTAATCGAGGAGTCACACAAGAGCTTCCATTAGCTCTTAGTACCTCCGATCTATCAAGCTGGGTTGGTTCGACTTGTAGCGAGGTATCCCCAACTCTGCAAACCGTGAATGAGATCTACCCGCTTTAGGCGTTCCGTCATGTATAATTCTGGACAAGTTAGACCACTTCGCCCGACTATAATCATCGCACCTGTAAGACCCAAATGTTAACATCTGTACAAGAGATGTTGGATACTGGTCACTTGTCGACATTAAGATTCTTAGCACAACAAGTAGTAGGCAGAGCAGCTGTCTATAAGCTTCATTTTTGCTCACGTATAAACAGAAATGTTGACTGGAGCTTGTCACGCTGAAGGGGTCTGTGGAATCGGGGCATCCCTTCAACCTCGTTGGAGTGGTTGTACCATACCGTGTGTAACAAACTGATGAGTATGGATTCATGTGCGACTTAGCAAATTATTACTACCTTCCACGCCGCAGTGGAGAACTCAGGTACCCTATAAATTGAAGGAGGCTGCGAGGGGCTGGGTCCTCAACCGACACACAATGATAAGCTCAACGTTTGCTTGGCCCAATCCCCTATCGTTACCTGAGAGAGTTTCGCTACCTCTTGGAAGCCACGAAAACATAAACGTAGGAAGGATTATACAATCCACTCGTGTGGCGCTTAGGAGTCCCAAAGCTACAAACCGGCAAGTATAGTGATATTTACACTAAGAAACTACATCATCGTAAACCACCGATGCGTGAGCGATAGGGTACCTTCTACGCGTTGATAGCAGGCCTAATCCGACTTGTCGACGAGAGGAATCGTTGAAGATCCAAAATTAGGGCGGTGGTCGCTTCAACCATGGCTCAGGTTCCACCAACCCGTGCTTTCTTCATAGCTTTACATGGATAGGGTGTCACCCGACTGTTAGCACCTCCTTCCGCCGTAAAGTCTCAGCTCTCTTATACGCCGATTCCACCGTGTCGGCCGGTCCTGAGAATTCTCTACGTGGAGTGAGCCGGCTTCCCTCGCCCGCAGTAAAGCTAGTCAGTGATATGGTGTCTCTTCACCCTACTACTCTTAAATAGGCCATTTAATAGTCACGTAGCCTGTTGGTCTTAGGTGAAGAGTACGCTCGGGGGAGCACCAGCAACTTTACAAACTGCCTTCGAACTTAGCTGATGGCTGTCTGGTGCCAAGAGTACTATCCGCTAACTGGGCCTCCCTGGTCAGTATGGCAATCCATCTGAGTTGGCAATACGCTTTAAATCGAAGCGTGAATCTAACGTCTCTGGGATGATAGTGCTCTGGGGGCAATCTCCAATACAAGAAGGTAGCGTGACACTCATATTGGAATGTTTGGCATGTTAACGCCATATTGGCGAAGCTGGTTAAAACAATCACGAATGGGCCCCCCTTTCCACTCGACCAGGCGACTCACGTGCCGGGACACGAGGACCCTACAGCAGATCGCATCCAGCTGTCGGCGAACAGGATCTGGGATTAAGGTAGCGAGCAAACCCCGTAAGCCGTGGTTTTGTGCTGAATGCTACGAATGCTTTAACATGCTCTCGATATACCCTTCAAGCCTTTGACGTCACCGTTCGACTGTCAATTCGGGATATGACTGGCTTCTTGAGATTGCCTTACTGTTCAAGCAAACGTTTAAGAGGCCCGCATAGCCAATAGCGGCCAGCCTCTGTTTGTAATGGGTCCCGCGTTCGAGGAAACCCCCAAAGGTGGTACAGGGTGCATGCTATAGGCATGAGGGTTAGGCACTTAAGCGACTGACCGAGCTTGATTGAAGCTCGTCCGGCACGATGGATATGTAAAGATTCTACGAGTTATTCGTCCCGCATAATGCCGGCCACCAGTTGTCGGTATTGAGCCCAAGGTGCAGCACAGGCTTGATTAAGCGACATGCAATTGGAGTGTTGACCATTAGGAGGGGGCTCTGCAAGGGTTCGTGTCTCGCATAAGGTGATCAGCCTTTCCTCTAGCACGGCGTTAAAGTTGCACAGAGGTTAAATTGGGCCTTGCAATCTCACCACCCTAGACTCCGAGTAAGGGCTGACAAACTCGGTTTTGGACACATCAGGTTCGTCGGTCCCGTTTATTTCAGCTGAGTCCACCGAGTTGTACCGGCGGTTAACGTGAGGAAATGGCGCCCATTTAAGCAGACCCCTAGTTACCAAGTATTGAATTGTTGTATGAGGGCCGATGGCCACCTATGATCTTTCGGACTTCTTCATTTGCCAAATCTCGATCACCTGAGAAAATTACCTGTCACGGCTTACTCATGGGTTACCTTTCCCATACCGGGGCGTTCCGAGTTCCTGGTGCATCAGAGACAAATACAGTTGAGTTATTTCGTACCTTGTTTTACCTCAGCGTGCTCGTCTGTGGACCGATCACGACATCGCCGGTTCTCCTGTACACTCGAGGCTATGGCCCTGACAAAAGTGCCCCCCCTCGAATGGCACTCAGTGAATCTCTTTTAAGTCCGTAATAGCAGCCCCATTGCTAAACGGCCAGGTACGTGGCGTTTCGCGTTTGCGGTCAGCATAGACATCCCCGCAAGCCGGGACCTCTCGACGCCTCGCGTCGATCGTGTGGTATACATTAATACGGGTCACCTGAGCTAAGACCTTCCAGGGCCTAACCCGCGACTTCCCCCTTTATGGACTTCCTTTTGCTATACTATAAAGTACTGAAGCCAAATGATTATACGCACTTGCCGCTTGGCCGTACACCCTCGGTCCCTAACATTTCCATTAGGACAGTTTGCAGTCCCCAAAGGTTCGTTCGCGAGTCGGATTATCCAGACTCCAGTCAAGCAATTTGGAAGCTACTTCTCCTAACATGATGAGAAGGTTAAGATCGGTCATTCAGGCTCTCCCCAATACACAAAGATTGCGGGCATACATAGGTGTATAAATCGCCGGGATGTGCGTTGAAATGAGTCTTTTATGTGAATCAGCTCCTAGATGGCACTTGCCGTTGTTTTAGCCGAGCGCTGCCGCAGACGCAACACTTAGACATTTATTTGGGTCTAATAGAAACTCCGAGCCGTGAGATTTCCCGCGATTGCGGGAATAGATTGTAAGGAACACAGATGCTTTGACAGGTGAGACTGTACGTCATCGTGTTGCAAACTGAACGATTTGGCAATATTAAACACTAACCCTACGGGAAATTTTGCTGCACGTGGCCCAGCCGGTGCTACCAAACCCTAGCTCCCATGACGACTCGTTTCACTAGTAGGAGATAATAACTCAGCGGAATCATGCGGTCCGTATTAGCAGTCTCATATATACGGGGGAAGTTCTCTTTGGTAGTACCTGCAGTGCCGGAGAAGGGCCCCGTACGACGCTGAATGGTTAGTCGAGGCTGAACTAACAGAGACAGCACACGTGCCAGTGTGTGGCTGAATCAACACCCTCCAGATATTCTACTGTCTATTTGAGCGATCACCTCGTTGTGCGAAGATGCACGGTCAGCCCTTGGTGTCTGCTAGTTTGCGCGGATTCGGGCGCAAGAATCACTGCTAAAAATGAGGTCCTTAGAATTTACACCCTTTCGCGCATAAGACAGGCAGGTGCGGACAAACTTTCGACACTTCGTTAGCGCTTACAGGCTCCGTTCGAAGGTCCCCGCATACGTGGTCGCACACATTAACGCGCCTGTTCGGTCAGGGCTTTATAGATCTCCCCCAAGTATGGTACGTCCCGTGCGCTATGAGATTCAAGTGACTGATGAGTAGGATATTGCGTACCAAACCTTTAAAAGCCTTCACAACAAGACTTATTGCGAGGAAATCCCAGGTA"
    pattern = "ATCTGCCACA"
    mismatches = 4
    indices = approximate_pattern_search(genome, pattern, mismatches)
    print(' '.join([str(i) for i in indices]))