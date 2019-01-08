"""
Guess (phred)format of FASTQ file.

The code is from:
https://github.com/brentp/bio-playground/blob/master/reads-utils/guess-encoding.py
"""

import operator
import optparse
import sys
import gzip
import io

from collections import Counter

import logging

log = logging.getLogger(__name__)

#  Note that the theoretical maximum for all encodings is 126.
#  The upper limits below are for "typical" data only.
RANGES = {
    'Sanger': (33, 73),
    'Illumina-1.8': (33, 74),
    'Solexa': (59, 104),
    'Illumina-1.3': (64, 104),
    'Illumina-1.5': (66, 105)
}

# The threshold to decide between Illumina-1.3 and Illumina-1.5
# based upon how common "B" is. The threshold insists it is
# within the Nth most common quality scores.
# N.B. needs to be conservative, as this is applied per input line.
N_MOST_COMMON_THRESH = 4


def get_qual_range(qual_str):
    """
    >>> get_qual_range("DLXYXXRXWYYTPMLUUQWTXTRSXSWMDMTRNDNSMJFJFFRMV")
    (68, 89...)
    """

    qual_val_counts = Counter(ord(qual_char) for qual_char in qual_str)

    min_base_qual = min(qual_val_counts.keys())
    max_base_qual = max(qual_val_counts.keys())

    return (min_base_qual, max_base_qual, qual_val_counts)


def get_encodings_in_range(rmin, rmax, ranges=RANGES):
    valid_encodings = []
    for encoding, (emin, emax) in ranges.items():
        if rmin >= emin and rmax <= emax:
            valid_encodings.append(encoding)
    return valid_encodings


def heuristic_filter(valid, qual_val_counts):
    """Apply heuristics to particular ASCII value scores
       to try to narrow-down the encoding, beyond min/max.
    """

    if 'Illumina-1.5' in valid:
        # 64–65: Phread+64 quality scores of 0–1 ('@'–'A')
        #        unused in Illumina 1.5+
        if qual_val_counts[64] > 0 or qual_val_counts[65] > 0:
            valid.remove('Illumina-1.5')

        # 66: Phread+64 quality score of 2 'B'
        #     used by Illumina 1.5+ as QC indicator
        elif 66 in map(operator.itemgetter(0),
                       qual_val_counts.most_common(N_MOST_COMMON_THRESH)):
            print("# A large number of 'B' quality scores (value 2, ASCII 66) "
                  "were detected, which makes it likely that this encoding is "
                  "Illumina-1.5, which has been returned as the only option.",
                  file=sys.stderr)
            valid = ['Illumina-1.5']

    return valid


def open_(path):
    if path.endswith(".gz"):
        fh = gzip.open(path)
        fh = io.TextIOWrapper(fh)
    else:
        fh = open(path)
    return fh


def guess_fq_format(input_fq, n=-1):
    gmin = 99
    gmax = 0
    valid = []

    with open_(input_fq) as input_file:
        for i, line in enumerate(input_file):
            if i % 4 != 3:
                continue
            lmin, lmax, qual_val_counts = get_qual_range(line.rstrip())

            if lmin < gmin or lmax > gmax:
                gmin, gmax = min(lmin, gmin), max(lmax, gmax)
                valid = get_encodings_in_range(gmin, gmax)

                valid = heuristic_filter(valid, qual_val_counts)

                if len(valid) == 0:
                    msg = "no encodings for range: " "{}".format((gmin, gmax))
                    raise IOError(msg)

                if len(valid) == 1 and n == -1:
                    # parsed entire file and found unique guess
                    break

            if n > 0 and i > n:
                # parsed up to specified portion; return current guess(es)
                break

    return valid


def guess_fq_phred(input_fq, default=33):
    try:
        valid = guess_fq_format(input_fq, -1)
    except IOError:
        log.warning("Can not inference FASTQ phred. use phred 33")
        return default
    if len(valid) != 1:
        log.warning("Can not inference FASTQ phred. use phred 33")
        return default
    else:
        fmt = valid[0]
        if fmt in {'Sanger', 'Illumina-1.8'}:
            return 33
        else:
            return 64