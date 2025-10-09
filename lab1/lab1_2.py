def letter_percentages(seq: str):
    total = len(seq)
    result = {}
    for ch in seq:
        result[ch] = result.get(ch, 0) + 1
    for ch in result:
        result[ch] = int((result[ch] / total) * 100)
    return result

print(letter_percentages("ACGGGCATATGCGC"))