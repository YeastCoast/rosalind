def generate_d_neighborhood(text, d):
    alphabet = {'A', 'C', 'G', 'T'}
    combinations = []
    text_len = len(text)

    def helper(start, array, comb=d, text_len=text_len):
        if len(array) == comb:
            combinations.append(array)
            return 0
        for i in range(start, text_len):
            new_array = array + [i]
            helper(i+1, new_array)

    def helper2(index, comb, output_str=""):
        if index == text_len:
            neighboords.append(output_str)
            return 0
        if index in comb:
            for nuc in alphabet:
                helper2(index+1, comb, output_str+nuc)
        else:
            helper2(index+1, comb, output_str+text[index])

    helper(0, [])

    neighboords = []
    for comb in combinations:
        comb = tuple(comb)
        helper2(0, comb)

    return set(neighboords)


if __name__ == '__main__':
    # test case
    text = "ACG"
    d = 1
    neighboord_list = '\n'.join(generate_d_neighborhood(text, d))
    print(neighboord_list)
