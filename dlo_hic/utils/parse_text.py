import re

def is_comment(line):
    if re.match("^[ \t]*#", line):
        return True
    else:
        return False


def add_chr_prefix(chr_a, chr_b):
    if not chr_a.startswith("chr"):
        chr_a_ = 'chr' + chr_a
    else:
        chr_a_ = chr_a

    if not chr_b.startswith("chr"):
        chr_b_ = 'chr' + chr_b
    else:
        chr_b_ = chr_b

    return chr_a_, chr_b_


def parse_line_bedpe(line, to_short=False):
    """ parse bedpe line 
    if to_short:
        return  chr_a, pos_a, chr_b, pos_b, score
    else:
        return chr_a, a_s, a_e, chr_b, b_s, b_e, name, score, strand_a, strand_b
    """
    items = line.strip().split()
    chr_a, a_s, a_e, chr_b, b_s, b_e, name, score, strand_a, strand_b = items[:10]
    score = int(score)
    a_s, a_e, b_s, b_e = [int(i) for i in (a_s, a_e, b_s, b_e)]
    if to_short:
        pos_a, pos_b = (a_s + a_e)//2, (b_s + b_e)//2
        chr_a, chr_b = add_chr_prefix(chr_a, chr_b)
        return chr_a, pos_a, chr_b, pos_b, score
    else:
        return chr_a, a_s, a_e, chr_b, b_s, b_e, name, score, strand_a, strand_b


def parse_line_bed6(line):
    """ parse bed6 line
    return chr, start, end, name, score, strand
    """
    items = line.strip().split()
    chr_, start, end, name, score, strand = items[:6]
    start, end, score = [int(i) for i in (start, end, score)]
    return chr_, start, end, name, score, strand


def parse_line_short(line):
    """ parse 'short format interaction' line 
    return  chr_a, pos_a, chr_b, pos_b, score
    """
    items = line.strip().split()
    chr_a, pos_a, chr_b, pos_b, score = items[:5]
    pos_a, pos_b, score = int(pos_a), int(pos_b), int(score)
    chr_a, chr_b = add_chr_prefix(chr_a, chr_b)
    return chr_a, pos_a, chr_b, pos_b, score
