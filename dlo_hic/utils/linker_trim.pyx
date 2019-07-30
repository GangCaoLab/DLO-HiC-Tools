from collections import Counter

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
        cutting_idx, rest_seq = rest
        self.rest_left  = rest_seq[:cutting_idx]
        self.rest_mid   = rest_seq[cutting_idx:-cutting_idx]
        self.rest_right = rest_seq[-cutting_idx:]
        self._rest_left_side  = self.rest_left + self.rest_mid
        self._rest_right_side = self.rest_mid  + self.rest_right
        self.mismatch = mismatch
        self._max_err_rate = float(mismatch) / len(linkers[0])  # assume all linker have same length
        self._min_overlap = len(linkers[0]) - mismatch
        self.mismatch_adapter = mismatch_adapter
        if self._trim_adapter:
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
        cdef int exact_start, len_linker
        cdef int idx
        cdef float error_rate
        cdef bint unmatch
        cdef int min_cost, min_idx

        alignment_adapter = None

        seq = fq_rec.seq
        aligner = Aligner(seq, 1)
#        aligner.min_overlap = self._min_overlap

        # align linker with sequence
        if self._single_linker:

            # try exact match
            exact_start = seq.find(self.linkers[0])
            if exact_start >= 0:  # exact matched
                len_linker = len(self.linkers[0])
                alignment = (exact_start, exact_start+len_linker, 0, len_linker, len_linker, 0)
                flag |= 2
            else:
                # do alignment
                alignment = aligner.locate(self.linkers[0])
                error_rate = alignment[-1] / len(self.linkers[0])
                if error_rate > self._max_err_rate:
                    # unmatch
                    flag |= 1
                    return (fq_rec, flag, None, None, alignment, alignment_adapter)
                else:
                    # AA matched
                    flag |= 2

        else:  # multiple linkers
            unmatch = 1
            min_idx = -1
            min_cost = 10000
            for idx in range(4):
                # try exact match
                exact_start = seq.find(self.linkers[idx])
                if exact_start >= 0:  # exact matched
                    unmatch = 0
                    min_idx = idx
                    len_linker = len(self.linkers[idx])
                    alignment = (exact_start, exact_start+len_linker, 0, len_linker, len_linker, 0)
                    break
                alignment = aligner.locate(self.linkers[idx])  # do alignment
                error_rate = alignment[-1] / len(self.linkers[0])
                if error_rate <= self._max_err_rate:  # matched
                    unmatch = 0
                    if alignment[-1] < min_cost:
                        min_idx = idx
                        min_cost = alignment[-1]
            if unmatch == 1:
                flag |= 1
                return (fq_rec, flag, None, None, alignment, alignment_adapter)
            elif min_idx < 2:  # AA/BB matched
                flag |= 2

        m_start, m_end, s_, e_, m_, cost = alignment
        pet1_end = m_start
        pet2_start = m_end

        # align adapter with PET2 sequence
        if self._trim_adapter:
            pet2_seq = seq[pet2_start:]
            aligner = Aligner(pet2_seq, 1)
#            aligner.min_overlap = self._min_overlap_adapter
            alignment_adapter = aligner.locate(self.adapter)
            if not alignment_adapter or alignment_adapter[-1]/len(self.adapter) >= self._max_err_rate_adapter:
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

        return (fq_rec, flag, PET1, PET2, alignment, alignment_adapter)


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
    "total",
]


def init_counts():
    return {
        'flag': np.zeros(len(COUNT_ITEM_NAMES), dtype=np.int64),
        'pet1_len': Counter(),
        'pet2_len': Counter(),
        'linker_match_score': Counter(),
        'adapter_match_score': Counter(),
    }


cpdef process_chunk(list chunk, tuple args):
    """
    Process a chunk of fastq records.
    """

    cdef int flag
    cdef object PET1, PET2

    linker_trimer = LinkerTrimer(*args)

    counts = init_counts()

    out_chunk = []
    for fq_rec in chunk:
        fq_rec, flag, PET1, PET2, align, align_ada = linker_trimer.trim(fq_rec)

        # count flags
        if flag & 1 == 0:
            if flag & 2 != 0:
                counts['flag'][1] += 1
            else:
                counts['flag'][2] += 1

            if flag & 4 != 0:
                counts['flag'][3] += 1
            if flag & 8 != 0:
                counts['flag'][4] += 1
            if flag & 16 != 0:
                counts['flag'][5] += 1
            if flag & 32 != 0:
                counts['flag'][6] += 1
            if flag & 64 != 0:
                counts['flag'][7] += 1
            if flag & 128 != 0:
                counts['flag'][8] += 1
            if flag & 256 != 0:
                counts['flag'][9] += 1

            if (flag & 1 == 0) and (flag & 32 == 0) and (flag & 128 == 0):
                counts['flag'][10] += 1
        else:
            counts['flag'][0] += 1

        counts['flag'][11] += 1

        # count pet length
        if PET1 or PET2:
            if (flag & 32 == 0) and (flag & 128 == 0):
                counts['pet1_len'].update({str(len(PET1.seq)): 1})
                counts['pet2_len'].update({str(len(PET2.seq)): 1})

        # count linker match score
        if align:
            l_score = align[-2] - align[-1]
            counts['linker_match_score'].update({str(l_score): 1})
        
        # count adapter match score
        if align_ada:
            a_score = align_ada[-2] - align_ada[-1]
            counts['adapter_match_score'].update({str(a_score): 1})

        out_chunk.append( (fq_rec, flag, PET1, PET2, align, align_ada) )

    return out_chunk, counts
