import numpy as np

from dlo_hic.utils.fastqio import Fastq
from dlo_hic.utils.align import Aligner


cdef class LinkerTrimer:

    cdef tuple linkers
    cdef bint  _single_linker
    cdef str adapter
    cdef bint _trim_adapter
    cdef str rest_left, rest_mid, rest_right
    cdef str _rest_left_side, _rest_right_side
    cdef int   mismatch
    cdef float _max_err_rate
    cdef int   _min_overlap
    cdef int   mismatch_adapter
    cdef float _max_err_rate_adapter
    cdef int   _min_overlap_adapter
    cdef int pet_len_max
    cdef int pet_len_min
    cdef int pet_cut_len

    def __init__(self, tuple linkers, str adapter, tuple rest,
                  int mismatch, int mismatch_adapter,
                  tuple pet_len_range, int pet_cut_len):
        self.linkers = linkers
        self._single_linker = linkers[1] is None
        self.adapter = adapter
        self._trim_adapter = not (adapter == "")
        self.rest_left = rest[0]
        self.rest_mid = rest[1]
        self.rest_right = rest[2]
        self._rest_left_side  = rest[0] + rest[1]
        self._rest_right_side = rest[1] + rest[2]
        self.mismatch = mismatch
        self._max_err_rate = float(mismatch) / len(linkers[0])
        self._min_overlap = len(linkers[0]) - mismatch
        self.mismatch_adapter = mismatch_adapter
        self._max_err_rate_adapter = float(mismatch_adapter) / len(adapter)
        self._min_overlap_adapter = len(adapter) - mismatch_adapter
        self.pet_len_min = pet_len_range[0]
        self.pet_len_max = pet_len_range[1]
        self.pet_cut_len = pet_cut_len

    def trim(self, object fq_rec):
        """
        Perform adapter cut and linker trim.

        Return
        ------
        tuple, fields:

            fq_rec : Fastq

            flag : int

                Bit    Description
                ------------------
                1      linker unmatchable
                2      intra-molecular linker (AA, BB)
                4      adapter unmatchable
                8      added base to left PET
                16     added base to right PET
                32     left PET len less than threshold
                64     left PET len large than threshold
                128    right PET len less than threshold
                256    right PET len large than threshold

            PET1 : Fastq
            PET2 : Fastq

        """
        cdef int flag = 0
        cdef object PET1, PET2
        cdef str pet1_seq, pet2_seq, pet1_quality, pet2_quality
        cdef str seq, quality
        cdef int m_start, m_end, a_start, a_end, s_, e_, m_
        cdef int pet2_end = -1
        cdef int pet1_start = 0
        cdef int pet1_end, pet2_start
        cdef int cost, a_cost
        cdef float min_err_rate
        cdef tuple alignment, alignment_adapter
        cdef int i

        seq = fq_rec.seq
        aligner = Aligner(seq, self._max_err_rate)
        aligner.min_overlap = self._min_overlap

        # align linker with sequence
        if self._single_linker:
            alignment = aligner.locate(self.linkers[0])
            if not alignment:
                # unmatch
                flag |= 1
                return (fq_rec, flag, None, None)
            else:
                # AA matched
                flag |= 2
        else:
            for idx in range(4):
                alignment = aligner.locate(self.linkers[idx])
                if alignment:
                    if idx < 2:  # AA/BB matched
                        flag |= 2
                    break
            else:  # all linker unmatch
                flag |= 1
                return (fq_rec, flag, None, None)

        m_start, m_end, s_, e_, m_, cost = alignment
        pet1_end = m_start
        pet2_start = m_end

        # align adapter with PET2 sequence
        if self._trim_adapter:
            pet2_seq = seq[pet2_start:]
            aligner = Aligner(pet2_seq, self._max_err_rate_adapter)
            aligner.min_overlap = self._min_overlap_adapter
            alignment_adapter = aligner.locate(self.adapter)
            if not alignment_adapter:
                flag |= 4
            else:
                a_start, a_end, s_, e_, m_, a_cost = alignment_adapter
                pet2_end = pet2_start + a_start

        # extract pet seq
        pet1_seq = seq[:pet1_end]
        if pet2_end > 0:
            pet2_seq = seq[pet2_start:pet2_end]
        else:
            pet2_seq = seq[pet2_start:]
            pet2_end = len(seq)

        # cut pet
        if len(pet1_seq) > self.pet_len_max:
            flag |= 64
            pet1_seq = pet1_seq[-self.pet_cut_len:]
            pet1_start = pet1_end - self.pet_cut_len
        elif len(pet1_seq) < self.pet_len_min:
            flag |= 32

        if len(pet2_seq) > self.pet_len_max:
            flag |= 256
            pet2_seq = pet2_seq[:self.pet_cut_len]
            pet2_end = pet2_start + self.pet_cut_len
        elif len(pet2_seq) < self.pet_len_min:
            flag |= 128

        # add base
        if pet1_seq[-len(self._rest_left_side):] == self._rest_left_side:
            flag |= 8
            pet1_seq = pet1_seq + self.rest_right
        if pet2_seq[:len(self._rest_right_side)] == self._rest_right_side:
            flag |= 16
            pet2_seq = self.rest_left + pet2_seq

        # construct pet
        if flag & 8 != 0:
            pet1_quality = fq_rec.quality[pet1_start:pet1_end] + len(self.rest_right) * fq_rec.quality[-1:]
        else:
            pet1_quality = fq_rec.quality[pet1_start:pet1_end]
        if flag & 16 != 0:
            pet2_quality = len(self.rest_left) * fq_rec.quality[0] + fq_rec.quality[pet2_start:pet2_end]
        else:
            pet2_quality = fq_rec.quality[pet2_start:pet2_end]

        PET1 = Fastq(fq_rec.seqid, pet1_seq, pet1_quality)
        PET2 = Fastq(fq_rec.seqid, pet2_seq, pet2_quality)

        return (fq_rec, flag, PET1, PET2)


COUNT_ITEM_NAMES = [
    "linker unmatchable",
    "intra-molecular linker",
    "inter-molecular linker",
    "adapter unmatchable",
    "added base to left PET",
    "added base to right PET",
    "left PET length less than threshold",
    "left PET length large than threshold",
    "right PET length less than threshold",
    "right PET length large than threshold",
    "valid reads",
    "all",
]


cpdef process_chunk(list chunk, tuple args):
    """
    Process a chunk of fastq records.
    """

    cdef int flag
    cdef object PET1, PET2

    linker_trimer = LinkerTrimer(*args)
    counts = np.zeros(len(COUNT_ITEM_NAMES), dtype=np.int64)

    out_chunk = []
    for fq_rec in chunk:
        fq_rec, flag, PET1, PET2 = linker_trimer.trim(fq_rec)
        if flag & 1 == 0:
            if flag & 2 != 0:
                counts[1] += 1
            else:
                counts[2] += 1

            if flag & 4 != 0:
                counts[3] += 1
            if flag & 8 != 0:
                counts[4] += 1
            if flag & 16 != 0:
                counts[5] += 1
            if flag & 32 != 0:
                counts[6] += 1
            if flag & 64 != 0:
                counts[7] += 1
            if flag & 128 != 0:
                counts[8] += 1
            if flag & 256 != 0:
                counts[9] += 1

            if (flag & 32 == 0) and (flag & 128 == 0):
                counts[10] += 1
        else:
            counts[0] += 1

        out_chunk.append( (fq_rec, flag, PET1, PET2) )
        counts[11] += 1
    return out_chunk, counts
