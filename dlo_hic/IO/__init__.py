from dlo_hic.utils.parse_text import is_comment

def load_chr_len(file_chr_len):
    """ 
    load a list of (chromesome, length) pair.

    :file_chr_len: a Tab split file record,
    the order and length(in basepair) of chrosome, file format like this:

    chr1    248956422
    chr2    242193529
    ...
    chrM    16569
    
    """
    chr_len = list()
    with open(file_chr_len) as f:
        for line in f:
            if is_comment(line):
                continue
            chr_, len_ = line.strip().split("\t") 
            len_ = int(len_)
            chr_len.append((chr_, len_))

    return chr_len
