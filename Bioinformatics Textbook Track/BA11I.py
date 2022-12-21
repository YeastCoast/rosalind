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
            size_dict[i][t] = sum([1/len(amino_acid_dict)*size_dict[i-amino_mass][t-spectrum[i-1]] for amino_mass in amino_acid_dict.values() if
                              (i - amino_mass) >= 0 and (max_score >= t - spectrum[i - 1] >= 0)])

    return round(sum(j for i, j in size_dict[len(spectrum)].items() if i >= threshold), 15)


if __name__ == '__main__':
    # test case 1
    spectrum = "4 -3 -2 3 3 -4 5 -3 -1 -1 3 4 1 3"
    spectrum_tup = spectrum_string_to_tuple(spectrum)
    threshold = 1
    max_score = 8

    size = get_size_of_spectral_dict(spectrum_tup, threshold, max_score, AA_DICT)
    print(size)

    # test case 2
    spectrum = "-10 11 3 10 11 12 -6 -5 4 4 -2 9 6 -8 9 -6 -1 10 -6 14 4 13 1 -6 5 -7 13 0 -1 12 -2 11 7 -10 9 13 14 -7 7 -9 -6 4 14 2 -9 1 12 13 15 6 15 13 -6 -10 -10 -8 -8 -7 -10 -7 -6 -4 6 9 -6 7 11 -1 -8 1 9 -5 6 7 -3 -10 -9 -1 4 7 7 -6 14 -6 12 15 7 8 11 -5 8 -8 12 -3 -1 -7 -6 9 13 12 -3 7 7 6 3 1 2 4 10 11 -10 -3 14 9 6 8 -9 1 5 -6 -8 5 -7 6 -6 -7 4 1 -3 7 5 10 11 12 0 -10 12 13 11 3 9 8 -10 9 -8 0 15 4 1 1 -4 12 2 4 0 15 -10 4 -10 -10 6 -5 -5 0 10 -5 8 1 14 6 -3 12 9 -7 -4 -9 -9 7 2 6 4 -10 -9 8 -4 -5 0 7 -4 -3 5 12 -10 3 -6 -10 6 10 -6 3 -5 15 4 14 -1 10 -9 13 11 -7 -5 -3 14 15 6 -3 -8 -5 0 12 0 12 2 8 -1 6 2 4 -6 3 11 -4 -10 1 -5 0 14 -5 -6 -1 15 13 12 -10 6 4 0 14 -1 5 15 13 4 -6 13 12 7 14 6 15 10 -9 1 -8 10 9 6 6 2 9 -2 5 11 -4 -6 -10 -7 10 9 8 -6 1 -8 2 -1 -1 -4 -2 0 9 11 -6 9 11 5 5 14 7 -10 14 -4 7 4 14 14 14 8 2 5 14 -4 13 7 10 14 -7 -6 11 -7 -2 -6 -3 1 -7 7 10 15 -6 -2 0 14 1 9 -7 5 -3 -5 5 -5 0 -4 1 3 11 9 -4 -3 -4 0 1 -4 15 -8 -3 0 0 11 -9 11 5 -9 1 -1 -7 -3 8 -9 11 5 4 4 -7 11 -1 -4 -5 7 -7 7 3 6 13 -1 11 -3 13 11 4 3 2 3 0 12 -6 3 12 -10 -8 -9 12 -2 12 5 -3 5 11 5 1 -2 3 5 1 11 6 -6 -2 0 -7 15 14 15 -10 0 6 13 9 10 -2 10 2 8 6 -6 5 -2 1 13 8 14 1 -4 11 11 -8 0 8 5 5 9 -1 -7 3 15 -7 -8 -3 11 9 0 10 2 1 13 4 0 -6 15 15 -1 10 3 1 2"
    spectrum_tup = spectrum_string_to_tuple(spectrum)
    threshold = 30
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
