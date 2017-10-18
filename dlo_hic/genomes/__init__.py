import re

from human import hg18, hg19, hg38
from mouse import mm8, mm9, mm10

def sort_chr_len(chr_len):
    num_chrlens = []
    other_chrlens = []
    for chr_, len_ in chr_len:
        if re.match("chr[0-9]+$", chr_):
            num_chrlens.append((chr_, len_))
        else:
            other_chrlens.append((chr_, len_))
    num_chrlens.sort(key=lambda (chr_, len_):int(chr_.replace("chr", "")))
    other_chrlens.sort()
    return num_chrlens + other_chrlens

supported_genomes = {
    'hg38': sort_chr_len(hg38),
    'hg19': sort_chr_len(hg19),
    'hg18': sort_chr_len(hg18),
    'mm10': sort_chr_len(mm10),
    'mm9' : sort_chr_len(mm9 ),
    'mm8' : sort_chr_len(mm8 ),
}
