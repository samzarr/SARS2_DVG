import bisect
import sys
import argparse

# third-party imports
import pysam

# local imports
import dvg

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("primer_bed", help="BED file containing primer coordinates")
    parser.add_argument("sam_file", help="SAM file containing alignments")
    args = parser.parse_args()

    primers = dvg.parse_primer_bed(args.primer_bed)
  
    input_samfile = pysam.AlignmentFile(args.sam_file, "r")
    with pysam.AlignmentFile("-", "w", template=input_samfile) as output_samfile:
        annotate_reads_with_primer_info(input_samfile, output_samfile, primers)

def annotate_reads_with_primer_info(input_samfile, output_samfile, primers):
    primer_overlapper = dvg.PrimerOverlapper(primers)
    primer_finder = dvg.PrimerFinder(primers)
    for read in input_samfile:
        if read.is_mapped:
            annotate_read_with_primer_overlap(read, primer_overlapper)
            annotate_read_with_closest_primer(read, primer_finder)
        output_samfile.write(read)

def annotate_read_with_primer_overlap(read, primer_overlapper):
    # blocks are 0-based coordinates
    blocks = read.get_blocks()
    overlap_length = sum(primer_overlapper.overlap_length(start, end) for start, end in blocks)
    query_length = read.infer_query_length()
    read.set_tag("ZO", overlap_length)
    read.set_tag("ZF", overlap_length / query_length)

def annotate_read_with_closest_primer(read, primer_finder):
    if not read.is_reverse:
        closest_primer = primer_finder.find_closest_upstream_primer(read.reference_start)
    else:        
        closest_primer = primer_finder.find_closest_downstream_primer(read.reference_end)

    if closest_primer is None:
        print("No closest primer found for read", read.query_name, read.reference_start, file=sys.stderr)
        return
    
    if not read.is_reverse:
        assert(closest_primer.start <= read.reference_start)
        read.set_tag("ZD", read.reference_start - closest_primer.start)
    else:
        assert(closest_primer.end >= read.reference_end)
        read.set_tag("ZD", closest_primer.end - read.reference_end)
    read.set_tag("ZN", closest_primer.name)
    read.set_tag("ZP", closest_primer.pool)

if __name__ == "__main__":
    main()
