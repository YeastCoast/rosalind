def pattern_to_number(text):
    value_dict = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    return sum(len(value_dict)**(i-1) * value_dict[text[len(text)-i]] for i in range(len(text), 0, -1))


if __name__ == '__main__':
    text = "CTTAATAGCGGATCAACTCCAG"
    print(pattern_to_number(text))
