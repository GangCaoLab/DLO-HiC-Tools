import re

def is_comment(line):
    if re.match("^[ \t]*#", line):
        return True
    else:
        return False


def parse_line_bedpe(line):
    """ parse bedpe line 
    return  chr_a, pos_a, chr_b, pos_b, val
    """
    items = line.strip().split()
    chr_a, a_s, a_e, chr_b, b_s, b_e, val = items[:7]
    a_s, a_e, b_s, b_e, val = [int(i) for i in (a_s, a_e, b_s, b_e, val)]
    pos_a, pos_b = (a_s + a_e)//2, (b_s + b_e)//2
    return chr_a, pos_a, chr_b, pos_b, val


def parse_line_short(line):
    """ parse 'short format interaction' line 
    return  chr_a, pos_a, chr_b, pos_b, val
    """
    items = line.strip().split()
    chr_a, pos_a, chr_b, pos_b, val = item[:5]
    return chr_a, pos_a, chr_b, pos_b, val
