import os
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


class Bedpe(object):
    """ The abstract of bedpe record. """
    def __init__(self, line):
        items = parse_line_bedpe(line)
        self.items = items
        self.chr1 = items[0]
        self.chr2 = items[3]
        self.name, self.score = items[6], items[7]
        self.start1, self.start2 = items[1], items[4]
        self.end1, self.end2 = items[2], items[5]
        self.strand1, self.strand2 = items[8], items[9]

        self.center1 = (self.start1 + self.end1) // 2
        self.center2 = (self.start2 + self.end2) // 2

        self.fields = ("chr1", "start1", "end1", "chr2", "start2", "end2", "name", "score", "strand1", "strand2")

    def is_rep_with(self, another, dis=10):
        """ Judge another bedpe record is replection of self or not. """
        if (self.chr1 != another.chr1) or (self.chr2 != another.chr2):
            return False
        else:
            if dis == 0:
                return (self.start1 == another.start1) and \
                       (self.start2 == another.start1) and \
                       (self.end1 == another.end1) and \
                       (self.end2 == another.end2) and \
                       (self.strand1 == another.strand1) and \
                       (self.strand2 == another.strand2)
            if abs(self.start1 - another.start1) > dis:
                return False
            if abs(self.end1 - another.end1) > dis:
                return False
            if abs(self.start2 - another.start2) > dis:
                return False
            if abs(self.end2 - another.end2) > dis:
                return False
            if self.strand1 != another.strand1:
                return False
            if self.strand2 != another.strand2:
                return False
        return True

    def exchange_pos(self):
        """
        exchange position of chr1 and chr2
        """
        self.chr1,    self.chr2    = self.chr2, self.chr1
        self.start1,  self.start2  = self.start2, self.start1
        self.end1,    self.end2    = self.end2, self.end1
        self.strand1, self.strand2 = self.strand2, self.strand1

    def to_upper_trangle(self):
        """
        rearrange item to upper trangle form.abs
        upper trangle:
                chr1    chr2    chr3 ...
        chr1     x       x       x   ...
        chr2             x       x   ...
        chr3                     x   ...
        ...                          ...
        """
        if self.chr1 > self.chr2:
            self.exchange_pos()
        elif self.chr1 == self.chr2:
            if self.start1 > self.start2:
                self.exchange_pos()                

    def __str__(self):
        line = "\t".join([self.chr1, str(self.start1), str(self.end1),
                          self.chr2, str(self.start2), str(self.end2),
                          self.name, str(self.score), self.strand1, self.strand2])
        return line

    @property
    def pos1(self):
        return self.end1

    @property
    def pos2(self):
        return self.start2

    def to_pairs_line(self, pos1='start', pos2='start'):
        """
        convert to pairs format line.
        about pairs format:
        https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md
        """
        self.to_upper_trangle()
        pos1 = self.start1 if pos1 == 'start' else self.end1
        pos2 = self.start2 if pos2 == 'start' else self.end2
        line = "\t".join([self.name, 
                          self.chr1, str(pos1),
                          self.chr2, str(pos2),
                          self.strand1, self.strand2])
        return line


def parse_line_pairs(line):
    """ parse pairs file format file line
    return read_id, chr1, pos1, chr2, pos2, strand1, strand2 """
    items = line.strip().split()
    read_id, chr1, pos1, chr2, pos2, strand1, strand2 = items
    pos1, pos2 = int(pos1), int(pos2)
    return read_id, chr1, pos1, chr2, pos2, strand1, strand2


class Pairs(object):
    """ The abstract of bedpe record. """
    def __init__(self, line):
        items = parse_line_pairs(line)
        self.items = items
        self.name = items[0]
        self.chr1, self.chr2 = items[1], items[3]
        self.pos1, self.pos2 = items[2], items[4]
        self.strand1, self.strand2 = items[5], items[6]

        self.fields = ("name", "chr1", "pos1", "chr2", "pos2", "strand1", "strand2")

    def __str__(self):
        return "\t".join([self.name,
                          self.chr1, str(self.pos1),
                          self.chr2, str(self.pos2),
                          self.strand1, self.strand2])

    def to_upper_trangle(self):
        """
        rearrange item to upper trangle form.abs
        upper trangle:
                chr1    chr2    chr3 ...
        chr1     x       x       x   ...
        chr2             x       x   ...
        chr3                     x   ...
        ...                          ...
        """
        if self.chr1 > self.chr2:
            self.chr1, self.pos1, self.strand1 = self.chr2, self.pos2, self.strand2
        if self.chr1 == self.chr2:
            if self.pos1 > self.pos2:
                self.pos1, self.strand1 = self.pos2, self.strand2

    def is_rep_with(self, another, dis=10):
        """ Judge another bedpe record is replection of self or not. """
        if (self.chr1 != another.chr1) or (self.chr2 != another.chr2):
            return False
        else:
            if dis == 0:
                return (self.pos1 == another.pos1) and \
                       (self.pos2 == another.pos2) and \
                       (self.strand1 == another.strand1) and \
                       (self.strand2 == another.strand2)
            span1 = abs(self.pos1 - another.pos1)
            span2 = abs(self.pos2 - another.pos2)
            if span1 > dis:
                return False
            if span2 > dis:
                return False
            if self.strand1 != another.strand1:
                return False
            if self.strand2 != another.strand2:
                return False
        return True


def infer_interaction_file_type(path):
    """
    Inference a interaction file in Pairs or Bedpe format.
    """
    from .filetools import open_file
    f_st = os.stat(path)
    if f_st.st_size == 0:
        raise IOError("Empty file")
    with open_file(path) as f:
        while True:  # skip header lines
            line = f.readline()
            if not is_comment(line):
                break
        else:
            raise IOError("{} not have any contents.")
    try:  # try Bedpe
        Bedpe(line)
        return Bedpe
    except:
        pass
    try:  # try Pairs
        Pairs(line)
        return Pairs
    except:
        pass
    raise NotImplementedError("Only support pairs and bedpe file format.")

            