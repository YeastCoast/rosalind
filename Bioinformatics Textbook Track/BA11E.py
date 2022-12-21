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
    input_str = "0 0 0 4 -2 -3 -1 -7 6 5 3 2 1 9 3 -8 0 3 1 2 1 8"
    input_lst = [int(i) for i in input_str.split(' ')]

    # for the test case on the rosalind webpage use VALUE_DICT; for exercise use value_aa_dict

    dag = create_dag(input_lst, VALUE_DICT)
    prefixes = bellman_ford(dag, input_lst, 0)
    protein_string = prefix_to_peptide([0] + prefixes[len(input_lst)], VALUE_DICT)
    print(protein_string)
