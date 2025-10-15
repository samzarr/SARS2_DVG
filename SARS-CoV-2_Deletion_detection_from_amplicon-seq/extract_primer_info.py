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
    parser.add_argument("--virema", dest="virema", action="store_true", default=False, help="Input is Virema output with paired end reads ending in _1 and _2")
    args = parser.parse_args()

    input_samfile = pysam.AlignmentFile(args.sam_file, "r")
    primer_info_table(input_samfile, 
                      args.virema)

table_fields = [
    "read_name", 
    "read_num",
    "is_supplementary",
    "chrom",
    "start",
    "end",
    "strand",
    "cigar",
    "closest_primer_name",
    "closest_primer_pool",
    "closest_primer_dist",
    "primer_overlap",
    "primer_overlap_frac"]

def primer_info_table(input_samfile, virema):
    print(*table_fields, sep="\t") 
    for read in input_samfile:
        is_aligned = read.reference_name is not None and read.reference_start is not None
        read_values = {
            "read_name": read.query_name if not virema else read.query_name[:-2],
            "read_num": (1 if read.is_read1 else 2) if not virema else int(read.query_name[-1]),
            "is_supplementary": "Y" if read.is_supplementary else "N",
            "chrom": read.reference_name if is_aligned else "",
            "start": read.reference_start if is_aligned else "",
            "end": read.reference_end if is_aligned else "",                
            "strand": ('+' if read.is_forward else '-') if is_aligned else "",
            "cigar": read.cigarstring if read.cigarstring is not None else "",
            "primer_overlap": read.get_tag("ZO") if read.has_tag("ZO") else "",
            "primer_overlap_frac": read.get_tag("ZF") if read.has_tag("ZF") else "",
            "closest_primer_name": read.get_tag("ZN") if read.has_tag("ZN") else "",
            "closest_primer_pool": read.get_tag("ZP") if read.has_tag("ZP") else "",
            "closest_primer_dist": read.get_tag("ZD") if read.has_tag("ZD") else ""}

        print_row(read_values, table_fields)

def print_row(values, fields, missing_string=""):
    print(*[values.get(field, missing_string) for field in fields], sep="\t")

if __name__ == "__main__":
    main()
