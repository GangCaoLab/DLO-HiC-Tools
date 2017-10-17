def get_base_map():
    base_map = ['\0' for i in range(256)]
    base_map[ ord('A') ] = 'T'
    base_map[ ord('T') ] = 'A'
    base_map[ ord('C') ] = 'G'
    base_map[ ord('G') ] = 'C'
    base_map[ ord('a') ] = 't'
    base_map[ ord('t') ] = 'a'
    base_map[ ord('c') ] = 'g'
    base_map[ ord('g') ] = 'c'
    base_map[ ord('N') ] = 'N'
    base_map[ ord('n') ] = 'n'
    base_map = bytes(''.join(base_map))
    return base_map

base_map = get_base_map()

def reverse_complement(seq):
    res = seq[::-1]
    res = res.translate(base_map)
    return res
