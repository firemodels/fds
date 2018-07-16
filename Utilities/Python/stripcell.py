def stripcell(CELL_ARRAY):
    C = []
    for elem in CELL_ARRAY:
        if isinstance(elem, str):
            C.append(elem)
    return C
