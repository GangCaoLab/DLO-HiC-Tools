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


def parse_line_bedpe(line):
    """ parse bedpe line 
    return  chr_a, pos_a, chr_b, pos_b, val
    """
    items = line.strip().split()
    chr_a, a_s, a_e, chr_b, b_s, b_e, val = items[:7]
    a_s, a_e, b_s, b_e, val = [int(i) for i in (a_s, a_e, b_s, b_e, val)]
    pos_a, pos_b = (a_s + a_e)//2, (b_s + b_e)//2
    chr_a, chr_b = add_chr_prefix(chr_a, chr_b)
    return chr_a, pos_a, chr_b, pos_b, val


def parse_line_short(line):
    """ parse 'short format interaction' line 
    return  chr_a, pos_a, chr_b, pos_b, val
    """
    items = line.strip().split()
    chr_a, pos_a, chr_b, pos_b, val = items[:5]
    pos_a, pos_b, val = int(pos_a), int(pos_b), int(val)
    chr_a, chr_b = add_chr_prefix(chr_a, chr_b)
    return chr_a, pos_a, chr_b, pos_b, val
