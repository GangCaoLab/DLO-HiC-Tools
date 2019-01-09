from collections import namedtuple

Fastq_ = namedtuple("Fastq", ["seqid", "seq", "quality"])

class Fastq(Fastq_):
    
    def __getitem__(self, key):
        if isinstance(key, slice):
            sub_rec = Fastq(
                self.seqid,
                self.seq.__getitem__(key),
                self.quality.__getitem__(key),
            )
            return sub_rec
        else:
            return super().__getitem__(key)


def read_fastq(path):
    """
    Read fastq file, yield fastq recored. Support Gziped file.
    """
    from .filetools import open_file
    with open_file(path) as f:
        for idx, line in enumerate(f):
            if idx % 4 == 0:
                seqid = line.strip()[1:]
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

