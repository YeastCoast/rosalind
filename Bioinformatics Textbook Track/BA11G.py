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


def protein_to_mass(peptide: str, amino_acid_dict: dict) -> tuple:
    """
    function to turn a peptide sequence into a tuple of the masses of each amino acid
    :param peptide: string representing the peptide
    :param amino_acid_dict: dictionary with one-letter amino acid code as keys and mass rounded to the nearest integer
                            as value
    :return: tuple containing the masses of each amino acid in the peptide
    """
    return tuple(amino_acid_dict[i] for i in peptide)


def masses_to_vector(masses: tuple) -> dict:
    """
    convert an object accessible by index containing amino acid masses to a binary vector
    :param masses: for example a tuple containing the masses of each amino acids of a peptide
    :return: dictionary with indices as keys and binary vector as value
    """
    return {i: [0 for _ in range(masses[i] - 1)] + [1] for i in range(len(masses))}


def score_peptide(vector: tuple, spectrum: tuple) -> int:
    """
    get the score of a protein vector against a spectrum, have to be equal length
    :param vector: a binary vector representing a peptide e.g. (0, 0, 0, 1, 0, 0, 0, 0, 1) -> XZ
    :param spectrum: tuple representing a protein spectrum
    :return: an integer representing the score of the vector against the spectrum
    """
    if len(vector) != len(spectrum):
        raise ValueError
    return sum(int(i)*int(j) for i, j in zip(vector, spectrum))


def scan_proteome(proteome: str, spectra: list, amino_acid_dict: dict, threshold: int) -> list:
    matches = ['' for _ in spectra]
    scores = [0 for _ in spectra]
    spectra_len = [len(i) for i in spectra]
    masses = protein_to_mass(proteome, amino_acid_dict)
    vector = masses_to_vector(masses)
    for i in range(len(proteome)):
        start_vector = vector[i]
        for j in range(i+1, len(proteome)):
            start_vector += vector[j]
            check_lens_eq = [len(start_vector) == spectrum_len for spectrum_len in spectra_len]
            if all([len(start_vector) > spectrum_len for spectrum_len in spectra_len]):
                break
            elif any(check_lens_eq):
                current_index = check_lens_eq.index(True)
                current_score = score_peptide(start_vector, tuple(spectra[current_index]))
                if current_score >= threshold and current_score > scores[current_index]:
                    best_match = proteome[i:j+1]
                    score = current_score
                    matches[current_index] = best_match
                    scores[current_index] = score
            else:
                pass

    return matches, scores


def parse_rosalind_input(file_path: str):
    with open(file_path, 'r') as f_r:
        all_lines = [i.strip() for i in f_r.readlines()]
    proteome = all_lines[-2]
    threshold = int(all_lines[-1])
    spectra = [[int(j) for j in i.split(' ')] for i in all_lines[:-2]]
    return spectra, proteome, threshold


def format_output(matches):
    output_str = '\n'.join([i for i in matches if i != ''])
    print(output_str)


if __name__ == '__main__':
    # test case on rosalind.info
    spectrum1 = "-1 5 -4 5 3 -1 -4 5 -1 0 0 4 -1 0 1 4 4 4"
    spectrum2 = "-4 2 -2 -4 4 -5 -1 4 -1 2 5 -3 -1 3 2 -3"
    spectra = [spectrum1, spectrum2]
    proteome = "XXXZXZXXZXZXXXZXXZX"
    threshold = 5
    spectr_tups = [spectrum_string_to_tuple(i) for i in spectra]
    matches, scores = scan_proteome(proteome, spectr_tups, AA_DICT, threshold)
    format_output(matches)
    print('\n')

    # extra test case on rosalind; you can paste your input into this file
    spectra, proteome, threshold = parse_rosalind_input("data/BA11G.txt")
    matches, scores = scan_proteome(proteome, spectra, massTableDictionary, threshold)
    format_output(matches)
    print('\n')

