from dlo_hic.utils.filetools import open_file


class FastaRec(object):
    def __init__(self, name=None, seq=None):
        self.name = name
        self.seq = seq

    def __str__(self):
        return ">"+self.name + "\n" + self.seq + "\n"


def read_fasta(file):
    if isinstance(file, str):
        file = open_file(file)

    rec = FastaRec()
    for line in file:
        line = line.strip()
        if line.startswith(">"):
            if rec.name is not None:
                yield rec
            rec = FastaRec(line[1:], "")
        else:
            rec.seq += line
    if rec.name is not None:
        yield rec

    file.close()

