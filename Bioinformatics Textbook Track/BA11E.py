from math import inf
from numpy import diff


massTableDictionary = {"A": 71,  "C": 103, "D": 115, "E": 129,
                       "F": 147, "G": 57,  "H": 137, "I": 113,
                       "K": 128, "L": 113, "M": 131, "N": 114,
                       "P": 97,  "Q": 128, "R": 156, "S": 87,
                       "T": 101, "V": 99,  "W": 186, "Y": 163}

value_aa_dict = {massTableDictionary[i]: i for i in massTableDictionary}
VALUE_DICT = {4: "X", 5: "Z"}


def create_dag(vector: list, masses: dict) -> dict:
    """
    function to create the directed acyclic graph from a protein spectrum vector
    :param vector: List of integers representing the protein spectrum.
    :param masses: Dictionary of with mass of amino acid as key and one-letter code as value.
    :return: Dictionary representing the constructed dag. Keys represent the index of the vector. The values are list
    of tuples representing edges in the graph. Each tuple consists of the Node the edge points to and the value of node
    it points to. E.g. {0: [(4, 0), (5, 0)], 1: [(5, 0), (6, 1)], ....}
    """

    vector = [0] + vector  # add sink to vector

    # create dag with dictionary comprehension, each entry of the vector is linked to all entries that are within the
    # distance of the mass of an amino acid, e.g. {4: X, 5: Z} 0 -> 4, 5; 1 -> 5, 6
    dag = {i: [(i+k, vector[i+k]) for k in masses if i+k < len(vector)] for i, j in enumerate(vector)}
    return {i: dag[i] for i in dag if len(dag[i]) != 0}


def bellman_ford(graph: dict, input_lst: list, start: int = 0) -> dict:
    """
    implementation of the bellman ford algorithm to find the maximum weighted path of the dag, which represents the
    optimal peptide of the given spectrum
    :param graph: dictionary representing a dag with following structure {node: [(next_node, value of next_node),..],..}
                  the nodes are a sequence of numbers from 0 to the length of the original protein spectrum vector
    :param start: start node
    :param input_lst: list representing the protein spectrum
    :return: dictionary representing the prefixes of the optimal peptide
    """
    weight_graph = [-inf for _ in [0]+input_lst]  # initialize
    weight_graph[start] = 0
    prefix = {start: []}
    for key, values in graph.items():  # iterate over the graph and fill the weighted graph
        for value in values:
            curr_weight = weight_graph[key] + value[1]
            if curr_weight > weight_graph[value[0]]:  # check if current edge represent a better path
                weight_graph[value[0]] = curr_weight
                prefix[value[0]] = prefix[key] + [value[0]]

    return prefix


def prefix_to_peptide(prefix: list, masses: dict) -> str:
    """
    convert a prefix to a peptide sequence
    :param prefix: List of integers representing the longest prefix of a protein sequence.
    :param masses: Dictionary of with mass of amino acid as key and one-letter code as value.
    :return: string representing the protein sequence corresponding to the prefix
    """
    return ''.join([masses[i] for i in diff(prefix)])


if __name__ == '__main__':
    # input the content from the file into input_str
    input_str = "-4 -8 27 16 17 -8 -20 -6 15 7 -11 1 11 21 -9 -9 -3 20 27 -14 -20 4 -15 2 18 -19 -6 -7 6 18 24 15 5 26 0 -1 -5 -11 20 16 20 -17 10 21 -11 27 -10 -2 26 17 0 -10 -7 -4 7 -5 -8 16 -15 -2 25 29 -15 -5 -18 5 0 13 18 4 25 -5 -3 -16 -5 2 -14 15 14 -9 8 7 18 12 -19 5 -9 -15 30 -10 24 23 17 30 13 26 3 12 -14 14 17 21 24 2 -11 12 -1 7 4 2 8 18 -17 -18 24 -13 8 -4 -13 30 -8 16 21 30 17 8 4 -11 -8 27 12 10 -15 -16 -20 -19 -10 -18 3 17 5 -12 -6 -15 -11 2 21 17 0 -13 15 6 4 5 -10 -4 -13 -2 17 -18 19 30 6 -3 18 14 5 20 20 -9 26 25 16 -6 -17 18 15 -10 -9 18 11 16 -11 -19 -11 -7 -1 -16 0 -20 3 19 -4 -11 12 29 1 -19 -18 -12 -8 -16 9 16 10 20 30 24 30 -5 9 -13 25 30 4 23 -5 11 12 27 10 -5 14 18 -17 -16 17 23 24 29 -3 -1 3 -12 17 30 -19 -1 20 22 17 8 28 -10 -1 10 17 8 16 -9 -2 8 29 -15 -5 4 18 23 26 20 -7 28 4 22 10 -1 20 -7 -2 -4 -13 -13 -9 -5 16 4 19 3 2 -7 13 -19 1 24 27 -4 14 -6 -14 8 15 21 -9 18 29 7 29 5 25 -3 28 0 1 20 10 4 19 -15 -16 -16 -7 19 15 26 30 11 3 30 14 21 29 -2 -15 9 0 -16 -20 2 -3 -16 2 26 15 15 9 -8 30 -14 30 5 -8 15 -16 20 28 -19 -18 4 20 -15 -2 12 -20 2 17 -20 18 -1 0 0 25 -4 14 2 7 -7 10 -2 -11 22 28 26 22 24 -12 -8 -7 -17 -15 16 -6 -16 2 28 0 25 6 7 27 -17 5 -2 17 -4 -4 4 -17 30 21 -11 -14 -10 15 -1 -7 -12 -5 15 13 19 19 19 9 -11 -19 16 -3 15 5 9 14 19 21 -19 -5 -10 16 24 21 -10 2 -17 22 3 13 -18 -4 26 30 19 30 -13 -3 -9 -8 30 29 -16 5 -7 23 -17 8 2 25 14 25 8 8 -16 -18 16 -4 -12 -4 22 6 0 -11 -9 -19 24 20 -9 7 -17 2 1 11 21 -18 -11 -2 -6 20 -9 -15 18 29 -13 20 5 17 21 15 7 12 14 -12 6 -1 0 18 -3 13 -16 26 17 -1 -2 20 0 21 -14 30 13 14 4 -18 25 26 1 -15 -11 -11 -4 1 -15 20 14 14 -17 -6 12 0 5 29 21 -14 -19 -5 4 -6 -6 30 2 -17 6 -12 11 3 -3 -4 8 12 26 7 -10 15 -2 24 27 -9 -2 -12 19 8 11 1 -12 24 -20 -16 -15 -13 8 21 -17 14 -4 21 16 -20 -3 -19 3 20 25 28 -9 3 16 7 4 -4 19 -14 14 4 -13 2 -9 12 5 18 25 -7 -4 -4 -15 -15 28 30 27 7 -10 -17 21 20 -19 -18 5 1 24 -20 -14 8 -18 17 9 -17 -14 -4 5 4 -20 2 8 16 -19 0 -9 1 6 20 -2 15 -4 0 30 -1 6 17 21 -9 -14 23 16 -1 -13 -1 25 4 -1 -11 -9 20 -12 -2 19 21 -7 -2 9 7 17 -12 12 -11 -17 7 22 -11 -20 -13 29 20 -1 13 18 -16 6 1 25 12 11 29 3 21 15 6 16 28 4 7 11 23 27 13 -20 28 3 10 -2 22 -8 21 -9 3 -6 -20 -3 8 -19 24 -16 -20 -17 17 27 16 30 -5 3 -3 -1 -17 -4 29 4 12 -17 4 -1 10 3 -14 20 -14 -13 -8 20 -19 -14 10 -18 9 29 18 -10 26 -8 18 7 26 -1 14 -19 15 -2 -13 -16 -2 5 5 -4 -7 12 30 26 -6 -7 1 30 -17 13 23 25 -13 20 7 1 -9 -4 0 8 0 29 -1 -7 25 -15 22 -8 -16 14 6 -3 -8 6 28 1 1 -20 11 -6 3 12 22 29 -19 20 -3 -11 20 -10 -16 -3 26 -3 -16 27 27 13 -19 14 26 4 23 20 26 6 -15 28 19 -7 -18 11 -5 5 -13 -8 30 6 24 -17 -4 13 12 -3"
    input_lst = [int(i) for i in input_str.split(' ')]

    # for the test case on the rosalind webpage use VALUE_DICT; for exercise use value_aa_dict

    dag = create_dag(input_lst, value_aa_dict)
    prefixes = bellman_ford(dag, input_lst, 0)
    protein_string = prefix_to_peptide([0] + prefixes[len(input_lst)], value_aa_dict)
    print(protein_string)
