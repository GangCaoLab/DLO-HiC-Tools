import os
import gzip
import io
import pytest

from dlo_hic.utils.fastqio import Fastq, read_fastq, write_fastq


fq1_seqid   = "test_reads1"
fq1_seq     = "AAAGGGCCCTTT"
fq1_quality = "#AAA<J<A7FFF"
fq2_seqid   = "test_reads2"
fq2_seq     = "TTTTAAAGGGCCCTTT"
fq2_quality = "AAAA#AAA<J<A7FFF"


def create_example_fq(fq="/tmp/example.fq", gz=False):
    if gz:
        fq = fq + ".gz"
        f = io.TextIOWrapper(gzip.open(fq, "wb"))
    else:
        f = open(fq, "w")

    f.write("@" + fq1_seqid + "\n")
    f.write(fq1_seq + "\n")
    f.write("+\n")
    f.write(fq1_quality + "\n")
    f.write("@" + fq2_seqid + "\n")
    f.write(fq2_seq + "\n")
    f.write("+\n")
    f.write(fq2_quality + "\n")

    f.close()
    return fq


def test_Fastq():
    fq = Fastq("test", "GAAAC", "AAAAA")
    assert fq.seqid == "test"
    assert fq.seq == "GAAAC"
    assert fq.quality == "AAAAA"
    with pytest.raises(AttributeError):
        fq.seqid = "other"
    with pytest.raises(AttributeError):
        fq.seq = "other"
    with pytest.raises(AttributeError):
        fq.quality = "other"


def test_read_fq():
    def test_(fq_file):
        fq_iter = read_fastq(fq_file)

        fq1 = next(fq_iter)
        assert fq1.seqid == fq1_seqid 
        assert fq1.seq == fq1_seq
        assert fq1.quality == fq1_quality
        fq2 = next(fq_iter)
        assert fq2.seqid == fq2_seqid
        assert fq2.seq == fq2_seq
        assert fq2.quality == fq2_quality

        with pytest.raises(StopIteration):
            next(fq_iter)

    fq_file = create_example_fq(gz=False)
    test_(fq_file)
    os.remove(fq_file)
    fq_file = create_example_fq(gz=True)
    test_(fq_file)
    os.remove(fq_file)


def test_write_fq():
    fq1_rec = Fastq(fq1_seqid, fq1_seq, fq1_quality)
    fq2_rec = Fastq(fq2_seqid, fq2_seq, fq2_quality)

    ex_fq_file = "/tmp/example.fq"

    with open(ex_fq_file, 'w') as fh:
        write_fastq(fq1_rec, fh)
        write_fastq(fq2_rec, fh)

    exp_fq_file = "/tmp/example.fq.exp"
    create_example_fq(exp_fq_file, gz=False)
    with open(ex_fq_file) as f1, open(exp_fq_file) as f2:
        assert f1.read() == f2.read()

    with io.TextIOWrapper(gzip.open(ex_fq_file, 'wb')) as fh:
        write_fastq(fq1_rec, fh)
        write_fastq(fq2_rec, fh)

    exp_fq_file = "/tmp/example.fq.exp"
    exp_fq_file = create_example_fq(exp_fq_file, gz=True)
    with gzip.open(ex_fq_file) as f1, gzip.open(exp_fq_file) as f2:
        assert f1.read() == f2.read()