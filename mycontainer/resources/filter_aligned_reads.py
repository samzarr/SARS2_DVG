# standard library imports
import argparse
import sys
import collections

# third-party imports
import pysam

# local imports
import dvg

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("primer_bed", 
                        help="BED file containing primer coordinates")
    parser.add_argument("sam_file", 
                        help="SAM/BAM file containing alignments")
    parser.add_argument("--sam", 
                        dest="output_sam", 
                        action="store_true", default=False, 
                        help="Output SAM format instead of BAM")    
    parser.add_argument("--min-deletion-length", 
                        dest="min_deletion_length", 
                        metavar="LENGTH", 
                        type=int, action="store", default=0, 
                        help="Minimum length of deletion to be considered")
    parser.add_argument("--max-overhang-primer-frac", 
                        dest="max_overhang_primer_frac", 
                        metavar="FRACTION", 
                        type=float, action="store", default=1, 
                        help="Reads with deletions that have overhangs with primer overlap fraction greater than this value will be excluded")
    parser.add_argument("--min-aligned-length", 
                        dest="min_aligned_length", 
                        metavar="LENGTH", 
                        type=int, action="store", default=0, 
                        help="Reads with aligned length less than this value will be excluded")
    parser.add_argument("--primer-pool-matching",
                        dest="primer_pool_matching",
                        action="store_true", default=False,
                        help="Require that pairs of reads come from the same primer pool")
    parser.add_argument("--max-primer-dist", default=1000, type=int,
                        help="Maximum distance of read start to the closest primer start to be considered originating from that primer")
    parser.add_argument("--virema", 
                        dest="virema", 
                        action="store_true", default=False, 
                        help="Input is Virema output with paired end reads ending in _1 and _2")

    args = parser.parse_args()

    primers = dvg.parse_primer_bed(args.primer_bed)

    input_samfile = pysam.AlignmentFile(args.sam_file, "r")
    read_ids_to_exclude = find_reads_to_exclude(input_samfile, primers, args)
    input_samfile.close()

    all_read_ids_to_exclude = set()
    for reason, read_ids in read_ids_to_exclude.items():
        print(f"Filtering out {len(read_ids)} reads for reason: {reason}", file=sys.stderr)
        all_read_ids_to_exclude.update(read_ids)
    print("Filtering out %d total reads" % len(all_read_ids_to_exclude), file=sys.stderr)

    input_samfile = pysam.AlignmentFile(args.sam_file, "r")
    output_samfile = pysam.AlignmentFile("-", "w" if args.output_sam else "wb", template=input_samfile)
    filter_reads(input_samfile, output_samfile, all_read_ids_to_exclude, args.virema)

def find_reads_to_exclude(input_samfile, primers, args):
    read_ids_to_exclude = collections.defaultdict(set)
    primer_overlapper = dvg.PrimerOverlapper(primers)
    check_for_primer_associated_deletion = args.max_overhang_primer_frac < 1
    primer_pools = {read_num: dict() for read_num in (1,2)}
    reason_counts = collections.defaultdict(int)
    for read in input_samfile:
        if read.query_alignment_length < args.min_aligned_length:
            read_ids_to_exclude["alignment length to short"].add(dvg.read_id(read, args.virema))
        elif (check_for_primer_associated_deletion and
              has_primer_associated_deletion(read, primer_overlapper, args)):
            read_ids_to_exclude["has primer associated deletion"].add(dvg.read_id(read, args.virema))
        
        if args.primer_pool_matching:
            update_primer_pools(read, primer_pools, args)

    if args.primer_pool_matching:
        primer_pool_read1 = set(primer_pools[1].keys())
        primer_pool_read2 = set(primer_pools[2].keys())
        read_ids_to_exclude["mismatching primer pool"].update(primer_pool_read1 ^ primer_pool_read2)
        for read_id in (primer_pool_read1 & primer_pool_read2):
            if primer_pools[1][read_id] != primer_pools[2][read_id]:
                read_ids_to_exclude["mismatching primer pool"].add(read_id)

    return read_ids_to_exclude

def has_primer_associated_deletion(read, primer_overlapper, args):
    skips = dvg.read_skip_intervals(read, args.min_deletion_length)
    for start, end, left_overhang, right_overhang in skips:
        left_primer_overlap = primer_overlapper.overlap_length(start - left_overhang, start)
        right_primer_overlap = primer_overlapper.overlap_length(end, end + right_overhang)
        if (left_primer_overlap / left_overhang > args.max_overhang_primer_frac or
            right_primer_overlap / right_overhang > args.max_overhang_primer_frac):
            return True
    return False

def update_primer_pools(read, primer_pools, args):
    read_id = dvg.read_id(read, args.virema)
    read_num = dvg.read_num(read, args.virema)
    pool = read_primer_pool(read, args)
    if read_id in primer_pools[read_num]:
        prior_pool = primer_pools[read_num][read_id]
        print(f"Warning: read {read_id} end {read_num} has multiple alignments with different primer pools", file=sys.stderr)
    primer_pools[read_num][read_id] = pool

def read_primer_pool(read, args):
    if (read.has_tag("ZD") and
        read.get_tag("ZD") <= args.max_primer_dist and
        read.has_tag("ZP")):
        return read.get_tag("ZP")
    else:
        return None

def filter_reads(input_samfile, output_samfile, read_ids_to_exclude, virema):
    for read in input_samfile:
        if dvg.read_id(read, virema) not in read_ids_to_exclude:
            output_samfile.write(read)

if __name__ == "__main__":
    main()
