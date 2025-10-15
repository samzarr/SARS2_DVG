# standard library imports
import argparse
import sys

# third-party imports
import pysam

# local imports
import dvg
import fasta

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("genome_fasta", help="FASTA file containing the reference genome")
    parser.add_argument("sam_file", 
                        help="SAM/BAM file containing alignments")
    parser.add_argument("--sam", 
                        dest="output_sam", 
                        action="store_true", default=False, 
                        help="Output SAM format instead of BAM")    

    args = parser.parse_args()

    genome = fasta_to_dict(args.genome_fasta)

    input_samfile = pysam.AlignmentFile(args.sam_file, "r")
    output_samfile = pysam.AlignmentFile("-", "w" if args.output_sam else "wb", template=input_samfile)
    standardize_alignments(input_samfile, output_samfile, genome)

def fasta_to_dict(fasta_filename):
    fasta_dict = {}
    with open(fasta_filename, "r") as f:
        for record in fasta.Reader(f):
            fasta_dict[record.title.split()[0]] = record.sequence
    return fasta_dict

def standardize_alignments(input_samfile, output_samfile, genome):
    for read in input_samfile:
        standardize_alignment(read, genome)
        output_samfile.write(read)

def standardized_shift(start, end, seq):
    shift = 0
    while (end + shift < len(seq)) and seq[start + shift] == seq[end + shift]:
        shift += 1
    return shift

def has_deletion(read):
    return 'D' in read.cigarstring or 'N' in read.cigarstring

def standardize_alignment(read, genome):
    # only reads that are mapped and have deletions require standardization
    if read.is_unmapped or not has_deletion(read):
        return read

    new_cigartuples = standardize_cigartuples(read.cigartuples,
                                              read.reference_start, 
                                              genome[read.reference_name])

    old_cigar = read.cigarstring
    read.cigartuples = new_cigartuples
    new_cigar = read.cigarstring

    if old_cigar != new_cigar:
        print("Shifted alignment", 
              "read =", read.query_name, 
              "old =", old_cigar, 
              "new =", new_cigar, 
              file=sys.stderr)

def standardize_cigartuples_stack(cigartuples, ref_start, ref_seq):
    in_stack = []
    out_stack = []
    for op, length in reversed(cigartuples):
        push_tuple(in_stack, op, length)
    while in_stack:
        op, length = in_stack.pop()
        if op in (pysam.CDEL, pysam.CREF_SKIP):
            shift = standardized_shift(ref_start, ref_start + length, ref_seq)
            if shift > 0:
                temp_stack = []
                while shift > 0 and in_stack:
                    next_op, next_length = in_stack.pop()
                    next_ref_length = cigar_ref_length(next_op, next_length)
                    if next_ref_length > shift:
                        push_tuple(temp_stack, next_op, shift)
                        push_tuple(in_stack, next_op, next_length - shift)
                        shift = 0
                    else:
                        push_tuple(temp_stack, next_op, next_length)
                        shift -= next_ref_length
                if in_stack:
                    push_tuple(in_stack, op, length)
                while temp_stack:
                    push_tuple(in_stack, *temp_stack.pop())
                continue
        push_tuple(out_stack, op, length)
        ref_start += cigar_ref_length(op, length)
    return out_stack

def push_tuple(stack, op, length):
    if stack and stack[-1][0] == op:
        stack[-1] = (op, stack[-1][1] + length)
    else:
        stack.append((op, length))

def standardize_cigartuples(cigartuples, ref_start, ref_seq):
    if not cigartuples:
        return []
    op, length = cigartuples[0]
    if op in (pysam.CDEL, pysam.CREF_SKIP):
        shift = standardized_shift(ref_start, ref_start + length, ref_seq)
        if shift > 0:
            prefix, suffix = split_cigartuples(cigartuples[1:], shift)
            if suffix:
                reordered = compact_cigartuples(prefix + [(op, length)] + suffix)
            else: # no need for deletion
                reordered = prefix
            return standardize_cigartuples(reordered, ref_start, ref_seq)

    rest = standardize_cigartuples(cigartuples[1:],
                                   ref_start + cigar_ref_length(op, length),
                                   ref_seq)
    return compact_cigartuples([(op, length)] + rest)

def split_cigartuples(cigartuples, ref_length):
    if ref_length == 0 or not cigartuples:
        return [], cigartuples
    op, length = cigartuples[0]
    op_ref_length = cigar_ref_length(op, length)
    if op_ref_length > ref_length:
        return [(op, ref_length)], [(op, op_ref_length - ref_length)] + cigartuples[1:]
    else:
        prefix, suffix = split_cigartuples(cigartuples[1:], ref_length - op_ref_length)
        return [(op, length)] + prefix, suffix

def cigar_ref_length(op, length):
    return length if op in dvg.BAM_REF_CONSUMING else 0

def compact_cigartuples(cigartuples):
    """Combine adjacent operations of the same type in a CIGAR string"""
    compacted = []
    i = 0
    while i < len(cigartuples):
        op, length = cigartuples[i]
        j = i + 1
        while j < len(cigartuples) and cigartuples[j][0] == op:
            length += cigartuples[j][1]
            j += 1
        compacted.append((op, length))
        i = j
    return compacted

if __name__ == "__main__":
    main()
