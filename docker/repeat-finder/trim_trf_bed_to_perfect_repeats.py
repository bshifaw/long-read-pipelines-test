#!/usr/bin/env python3

"""
This script attempts to mimic the GangSTR paper's steps for generating a set of reference repeats, where, after they
run TandemRepeatsFinder, they attempt to trim the output intervals to remove any non-perfect repeats near the edges, and
leave a shorter sequence of perfect repeats.  The paper's methods section has this description:

To avoid errors in the local realignment step of GangSTR,
-- all repeating regions are trimmed until they no longer contain any imperfections in their first and last four
copies of the motif.
-- Next we require that the trimmed repeating region is a perfect repetition of the motif. This step ensures there
are no errors in longer STRs that may pass the trimming step.
-- Finally, we set a threshold of at least four surviving copies for motifs of length 2-8bp and at least three copies
for motifs of length greater than 8bp.
"""


import argparse
from collections import defaultdict
import gzip
import logging
import os
import pyfaidx
import tqdm
import unittest

logging.basicConfig(level=logging.INFO)

def init_args():
    p = argparse.ArgumentParser()
    p.add_argument("reference_fasta_path", help="genome reference fasta path")
    p.add_argument("input_bed_path", help="input .bed file generated by running TandemRepeatFinder, followed by convert_dat_output_to_bed.py")
    p.add_argument("output_bed_path", help="output .bed file")

    return p


def parse_args(argparser):
    args = argparser.parse_args()

    if not os.path.isfile(args.reference_fasta_path):
        argparser.error("File not found: {}".format(args.reference_fasta_path))

    if not os.path.isfile(args.input_bed_path):
        argparser.error("File not found: {}".format(args.bed_path))

    return args


def find_perfect_repeat(repeat_unit, nuc_sequence):
    """Attempts to find a sub-sequence within nuc_sequence that corresponds to 3 or 4 exact repeats of the repeat_unit.
    If found, the offset start and end positions of this sub-sequence (expanded by as many repeat_units as possible)
    are returned as a tuple of 2 integers.  If not found, a ValueError is raised.

    For example:
        find_perfect_repeat("GT", "ACTGTGTGTGTGTGTGTACTGTGTT") would return (3, 17)

        find_perfect_repeat("GT", "ACTGTGTGTA") would raise a ValueError because there's no sequence of 4 "GT" repeats


    Args:
        repeat_unit (str): a short sequence of nucleotides
        nuc_sequence (str): a long sequence of nucleotides

    Return:
        2-tuple: (start, end) positions of the first sequence of 3 or 4 exact repeats of the repeat_unit within nuc_sequence
    """

    if len(repeat_unit) <= 8:
        min_tandem_repeats = 4
    else:
        min_tandem_repeats = 3

    seed = repeat_unit * min_tandem_repeats

    match_start = nuc_sequence.find(seed)
    if match_start == -1:
        raise ValueError("seed (= {} * repeat_unit) not found in sequence".format(min_tandem_repeats))

    match_end = match_start + len(seed)

    # try to extend seed to the right
    while nuc_sequence[match_end: match_end + len(repeat_unit)] == repeat_unit:
        match_end += len(repeat_unit)

    return match_start, match_end


def main():
    args = parse_args(init_args())

    # read bed file
    fopen = gzip.open if args.input_bed_path.endswith(".gz") else open
    input_bed_file = fopen(args.input_bed_path)
    output_bed_file = fopen(args.output_bed_path, "w")

    counter = defaultdict(int)
    fasta = pyfaidx.Fasta(args.reference_fasta_path, one_based_attributes=False)
    for line in tqdm.tqdm(input_bed_file, unit=" lines"):
        if line.startswith("#"):
            continue

        counter["total input rows"] += 1

        # parse .bed line
        fields = line.strip().split()

        chrom = fields[0].replace("chr", "")
        start = int(fields[1])
        end = int(fields[2])
        repeat_unit = fields[3]

        if "_" in chrom or len(chrom) >= 5:
            continue

        counter["input rows in chromosomes"] += 1
        full_interval_sequence = str(fasta[chrom][start:end])

        try:
            start_offset, end_offset = find_perfect_repeat(repeat_unit, full_interval_sequence)
        except ValueError as e:
            counter[str(e)] += 1
            continue

        fields[1] = start + start_offset
        fields[2] = start + end_offset

        output_bed_file.write("\t".join(map(str, fields)) + "\n")
        counter["output rows"] += 1

    logging.info("Parsed {} lines from {}".format(counter["total input rows"], args.input_bed_path))
    for key, value in sorted(counter.items(), key=lambda x: x[1]):
        logging.info("   ==> {} {}".format(value, key))

    output_bed_file.close()
    logging.info("Wrote {} lines to {}".format(counter["output rows"], args.output_bed_path))

if __name__ == "__main__":
    main()


class Tests(unittest.TestCase):

    def test_compute_repeat_unit(self):

        self.assertRaises(ValueError,
            lambda nuc_sequence=3*"A" + 3*"C" + 5*"T": find_perfect_repeat("C", nuc_sequence))
        self.assertEqual(
            find_perfect_repeat("C", 3*"A" + 10*"C" + 5*"T"),
                (3, 3 + 10))
        self.assertEqual(
            find_perfect_repeat("CGT", "CGTCGTCGTCGTCGT"),
                (0, 5*3))
        self.assertEqual(
            find_perfect_repeat("CGT", "ACGTACGTCGTCGTCGTCGTTTTT"),
                (5, 5+5*3))
        self.assertEqual(
            find_perfect_repeat("CGTACAGAT", "CGTACAGATCGTACAGATCGTACAGAT"),
                (0, 9*3))


