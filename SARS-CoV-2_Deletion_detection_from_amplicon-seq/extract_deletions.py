# standard library imports
import argparse

# third-party imports
import pysam

# local imports
import dvg
import fasta

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("sam_file", help="SAM/BAM file containing alignments")    
    parser.add_argument("--primer-bed", help="BED file containing primer coordinates")
    parser.add_argument("--min-deletion-length", dest="min_deletion_length", metavar="LENGTH", type=int, action="store", default=0, help="Minimum length of deletion to be considered")
    parser.add_argument("--ignore-secondary", dest="ignore_secondary", action="store_true", default=False, help="Ignore secondary alignments")
    parser.add_argument("--virema", dest="virema", action="store_true", default=False, help="Input is Virema output with paired end reads ending in _1 and _2")
    parser.add_argument("--output-reads-without-deletions", dest="output_reads_without_deletions", action="store_true", default=False, 
                        help="Output a line for each read that does not have a deletion")
    args = parser.parse_args()

    primers = dvg.parse_primer_bed(args.primer_bed) if args.primer_bed else None

    input_samfile = pysam.AlignmentFile(args.sam_file, "r")
    deletion_table(input_samfile, primers, args)

table_fields = [
    "read_name", 
    "read_num",
    "is_supplementary",
    "chrom",
    "read_start",
    "start",
    "end",    
    "strand",
    "left_overhang",
    "right_overhang",
    "left_overhang_primer_frac",
    "right_overhang_primer_frac",
    "cigar",
    "aligned_length"]

def deletion_table(input_samfile, primers, args):
    primer_overlapper = dvg.PrimerOverlapper(primers) if primers else None

    print(*table_fields, sep="\t") 
    for read in input_samfile:
        if args.ignore_secondary and read.is_secondary: continue
        skips = dvg.read_skip_intervals(read, args.min_deletion_length)
        is_aligned = read.reference_name is not None and read.reference_start is not None
        read_values = {
            "read_name": dvg.read_id(read, args.virema),
            "read_num": dvg.read_num(read, args.virema),
            "is_supplementary": "Y" if read.is_supplementary else "N",
            "chrom": read.reference_name if is_aligned else "",
            "read_start": read.reference_start if is_aligned else "",                   
            "strand": ('+' if read.is_forward else '-') if is_aligned else "",
            "cigar": read.cigarstring if read.cigarstring is not None else "",
            "aligned_length": read.query_alignment_length if is_aligned else ""}

        for start, end, left_overhang, right_overhang in skips:
            skip_values = read_values.copy()
            skip_values.update({
                "start": start, 
                "end": end,
                "left_overhang": left_overhang,
                "right_overhang": right_overhang})

            if primers:
                left_primer_overlap = primer_overlapper.overlap_length(start - left_overhang, start)
                right_primer_overlap = primer_overlapper.overlap_length(end, end + right_overhang)
                assert(left_primer_overlap <= left_overhang)
                assert(right_primer_overlap <= right_overhang)
                skip_values.update({
                    "left_overhang_primer_frac": left_primer_overlap / left_overhang,
                    "right_overhang_primer_frac": right_primer_overlap / right_overhang})

            print_row(skip_values, table_fields)

        if not skips and args.output_reads_without_deletions:
            print_row(read_values, table_fields)

def print_row(values, fields, missing_string=""):
    print(*[values.get(field, missing_string) for field in fields], sep="\t")

if __name__ == "__main__":
    main()
