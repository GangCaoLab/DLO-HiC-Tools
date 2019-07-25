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


def parse_line_bedpe(line, to_short=False, extends=False):
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

    if extends:
        return (chr_a, a_s, a_e, chr_b, b_s, b_e, name, score, strand_a, strand_b), items[10:]

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

    fields = ("chr1", "start1", "end1", "chr2", "start2", "end2", "name", "score", "strand1", "strand2")

    def __init__(self, chr1, start1, end1, chr2, start2, end2, name, score, strand1, strand2, extends=None):
        self.chr1 = chr1
        self.chr2 = chr2
        self.name, self.score = name, score
        self.start1, self.start2 = start1, start2
        self.end1, self.end2 = end1, end2
        self.strand1, self.strand2 = strand1, strand2

        self.items = [getattr(self, f) for f in self.fields]

        self.center1 = (self.start1 + self.end1) // 2
        self.center2 = (self.start2 + self.end2) // 2
        self.extends = extends
        if extends:
            self.etag1 = extends[0]
            self.etag2 = extends[1]

    @classmethod
    def from_line(cls, line):
        items, extends = parse_line_bedpe(line, extends=True)
        return cls(*items, extends=extends)

    def is_rep_with(self, another, dis=10, by_etag=False):
        """ Judge another bedpe record is replection of self or not. """
        if (self.chr1 != another.chr1) or (self.chr2 != another.chr2):
            return False
        if by_etag:
            return (self.etag1 == another.etag1)     and (self.etag2 == another.etag2) and \
                   (self.strand1 == another.strand1) and (self.strand2 == another.strand2)
        else:
            if dis == 0:
                return (self.start1 == another.start1) and \
                       (self.start2 == another.start2) and \
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
        fields = [self.chr1, str(self.start1), str(self.end1),
                  self.chr2, str(self.start2), str(self.end2),
                  self.name, str(self.score), self.strand1, self.strand2]
        if self.extends:
            line = "\t".join(fields + self.extends)
        else:
            line = "\t".join(fields)
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

    fields = ("name", "chr1", "pos1", "chr2", "pos2", "strand1", "strand2")

    def __init__(self, name, chr1, pos1, chr2, pos2, strand1, strand2):
        self.name = name
        self.chr1, self.chr2 = chr1, chr2
        self.pos1, self.pos2 = pos1, pos2
        self.strand1, self.strand2 = strand1, strand2

        self.items = [getattr(self, f) for f in self.fields]

    @classmethod
    def from_line(cls, line):
        items = parse_line_pairs(line)
        return Pairs(*items)

    def __str__(self):
        return "\t".join([self.name,
                          self.chr1, str(self.pos1),
                          self.chr2, str(self.pos2),
                          self.strand1, self.strand2])

    def exchange_pos(self):
        """
        exchange position of chr1 and chr2
        """
        self.chr1,    self.chr2    = self.chr2, self.chr1
        self.pos1,    self.pos2    = self.pos2, self.pos1
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
            if self.pos1 > self.pos2:
                self.exchange_pos()                

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
        Bedpe.from_line(line)
        return Bedpe
    except:
        pass
    try:  # try Pairs
        Pairs.from_line(line)
        return Pairs
    except:
        pass
    raise NotImplementedError("Only support pairs and bedpe file format.")


def parse_rest(rest):
    import re
    pattern = "[*^]"
    m = re.search(pattern, rest)
    if not m:
        raise ValueError("There are no cutting site in rest-seq")
    cutting = m.span()[0]
    seq = re.sub(pattern, '', rest)
    return cutting, seq
