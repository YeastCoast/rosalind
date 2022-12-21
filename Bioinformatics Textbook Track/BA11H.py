import numpy as np

massTableDictionary = {"A": 71,  "C": 103, "D": 115, "E": 129,
                       "F": 147, "G": 57,  "H": 137, "I": 113,
                       "K": 128, "L": 113, "M": 131, "N": 114,
                       "P": 97,  "Q": 128, "R": 156, "S": 87,
                       "T": 101, "V": 99,  "W": 186, "Y": 163}
value_aa_dict = {massTableDictionary[i]: i for i in massTableDictionary}
AA_DICT = {"X": 4, "Z": 5}
VALUE_DICT = {4: "X", 5: "Z"}

def spectrum_string_to_tuple(spectr_str: str) -> tuple:
    return tuple(int(i) for i in spectr_str.split(' '))


def get_size_of_spectral_dict(spectrum: tuple, threshold: int, max_score: int, amino_acid_dict: dict):
    size_dict = {i: {j: 0 for j in range(max_score+1)} for i in range(len(spectrum)+1)}
    size_dict[0][0] = 1

    for i in range(1, len(spectrum)+1):
        for t in range(max_score+1):
            size_dict[i][t] = sum([size_dict[i-amino_mass][t-spectrum[i-1]] for amino_mass in amino_acid_dict.values() if
                              (i - amino_mass) >= 0 and (max_score >= t - spectrum[i - 1] >= 0)])

    return sum(j for i, j in size_dict[len(spectrum)].items() if i >= threshold)


if __name__ == '__main__':
    # test case 1
    spectrum = "4 -3 -2 3 3 -4 5 -3 -1 -1 3 4 1 3"
    spectrum_tup = spectrum_string_to_tuple(spectrum)
    threshold = 1
    max_score = 8

    size = get_size_of_spectral_dict(spectrum_tup, threshold, max_score, AA_DICT)
    print(size)

    # test case 2
    spectrum = "14 -4 -3 -3 5 9 0 14 2 1 -4 6 -1 13 2 -5 13 -8 -8 3 0 -10 14 4 14 14 8 -8 1 3 -10 -2 2 -9 3 6 13 -10 6 -8 12 2 8 -1 -5 -6 -6 10 3 -3 12 -4 14 3 11 14 15 12 -7 -5 -2 11 13 -9 15 -8 -10 5 -8 5 6 -9 -2 7 -6 -1 -2 12 12 -3 9 0 3 0 5 -6 3 3 -7 6 0 -6 5 8 -7 5 3 13 13 -2 -9 0 2 13 13 12 7 -2 -10 -5 -7 7 13 11 14 -4 -9 15 -10 5 -7 -6 -7 -6 11 5 9 8 -4 7 1 -9 12 2 8 12 -6 0 2 -5 10 11 14 15 -1 3 -3 3 -3 12 15 4 -2 14 13 8 -10 2 -3 0 -6 8 3 10 0 9 10 13 15 6 9 -10 -9 1 -3 -10 8 1 -10 2 -1 14 -3 15 -1 0 1 6 -7 5 12 6 -9 2 1 -2 14 -5 1 -8 -6 11 -5 2 -3 -8 7 -6 -10 8 6 13 -8 -5 -10 12 -5 -9 8 9 0 10 15 -1 4 2 -8 9 1 -9 -6 -8 -1 -1 5 10 -4 7 3 11 4 12 6 6 13 -3 12 -3 1 7 11 6 13 8 3 -6 5 11 4 -1 15 10 -8 -7 0 4 7 5 -4 8 -3 -4 -8 9 -2 -3 13 1 12 4 -1 13 -1 -5 -5 7 7 -7 -5 6 6 -2 -5 7 10 14 11 12 -9 6 -3 4 15 -8 11 -3 -7 5 -4 7 9 15 -9 8 13 6 -2 -3 9 6 5 14 10 -7 -9 -8 10 2 -3 -1 2 3 12 13 6 -2 8 -5 5 -3 -8 10 3 0 12 -7 10 6 15 8 7 -2 8 14 -2 13 -1 8 15 -7 -7 -7 7 -3 -2 5 -4 -3 15 11 -4 9 11 13 15 8 4 -6 7 12 14 6 -10 -5 -9 4 -9 13 -3 0 12 3 12 -5 11 1 15 -8 5 3 -5 7 15 -2 -9 0 0 1 1 -1 -4 -1 5 12 12 -5 8 5 14 12 5 -9 2 -10 -9 4 -2 6 5 -3 -7 7 5 -8 -10 8 -1 7 3 6 -6 14 -8 6 -5 -8 -10 14 -2 12 4 5 -2 9 -4 1 -5 -3 -6 -8 -9 -10 10 4 9 11 -7 6 -4 5 13 -8 -7 -3 10 14 4 10 6 4 0 13 -3 11 -9 2 -8 6 -8 4 -1"
    spectrum_tup = spectrum_string_to_tuple(spectrum)
    threshold = 37
    max_score = 200

    size = get_size_of_spectral_dict(spectrum_tup, threshold, max_score, massTableDictionary)
    print(size)

    # user case
    spectrum = ""
    spectrum_tup = spectrum_string_to_tuple(spectrum)
    threshold = 0
    max_score = 0

    size = get_size_of_spectral_dict(spectrum_tup, threshold, max_score, massTableDictionary)
    print(size)
