
cdef class Fastq:

    cdef readonly str seqid, seq, quality

    def __init__(self, str seqid, str seq, str quality):
        self.seqid = seqid
        self.seq = seq
        self.quality = quality


def read_fastq(str path):
    """
    Read fastq file, yield fastq recored. Support Gziped file.
    """
    from .filetools import open_file
    cdef str seqid, seq, quality
    cdef size_t idx
    cdef str line
    with open_file(path) as f:
        for idx, line in enumerate(f):
            if idx % 4 == 0:
                seqid = line.strip().split()[0][1:]
            elif idx % 4 == 1:
                seq = line.strip()
            elif idx % 4 == 2:
                continue
            else:
                quality = line.strip()
                rec = Fastq(seqid, seq, quality)
                yield rec


def write_fastq(rec, fh):
    """
    Write Fastq record to file.
    """
    fh.write("@"+rec.seqid+"\n")
    fh.write(rec.seq+"\n")
    fh.write("+\n")
    fh.write(rec.quality+"\n")

