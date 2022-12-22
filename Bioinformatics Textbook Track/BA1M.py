def number_to_pattern(number, k):
    value_dict = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
    alphabet_len = len(value_dict)
    output_str = ""
    for i in range(k, 0, -1):
        key = number % alphabet_len
        number = int(number/alphabet_len)
        output_str = value_dict[key] + output_str
    return output_str


if __name__ == '__main__':
    integer = 6978
    k = 9
    print(number_to_pattern(integer, k))
