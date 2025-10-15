import csv
import bisect
import pysam
from collections import namedtuple

Primer = namedtuple("Primer", ["start", "end", "chrom", "name", "pool", "strand"])

def parse_primer_bed(filename):
    primers = []
    with open(filename, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            # these look 0-based
            primer = Primer(int(row[1]), int(row[2]), row[0], row[3], row[4], row[5])
            primers.append(primer)
    return primers

def read_id(read, virema):
    return read.query_name if not virema else read.query_name[:-2]

def read_num(read, virema):
    return (1 if read.is_read1 else 2) if not virema else int(read.query_name[-1])

BAM_REF_CONSUMING = (pysam.CREF_SKIP, 
                     pysam.CMATCH, 
                     pysam.CEQUAL, 
                     pysam.CDIFF,
                     pysam.CDEL)

BAM_MATCH = (pysam.CMATCH, 
             pysam.CEQUAL, 
             pysam.CDIFF)

class PrimerOverlapper:
    def __init__(self, primers):
        primer_intervals = [(primer.start, primer.end) for primer in primers]
        self.nonredundant_primers = nonredundant_intervals(primer_intervals)
        self.primer_starts, self.primer_ends = zip(*self.nonredundant_primers)
    def overlap_length(self, start, end):
        return primer_overlap(start, end, self.primer_starts, self.primer_ends)

class PrimerFinder:
    def __init__(self, primers):
        self.forward_primers = sorted(filter(lambda p: p.strand == '+', primers), key=lambda p: p.start)
        self.forward_primer_starts = [p.start for p in self.forward_primers]
        self.reverse_primers = sorted(filter(lambda p: p.strand == '-', primers), key=lambda p: p.end)
        self.reverse_primer_ends = [p.end for p in self.reverse_primers]
    def find_closest_upstream_primer(self, position):
        i = bisect.bisect_right(self.forward_primer_starts, position)
        if i == 0:
            return None
        return self.forward_primers[i - 1]
    def find_closest_downstream_primer(self, position):
        i = bisect.bisect_left(self.reverse_primer_ends, position)
        if i == len(self.reverse_primers):
            return None
        return self.reverse_primers[i]

def nonredundant_intervals(intervals):
    if not intervals:
        return []
    intervals.sort()
    nonredundant = [intervals[0]]
    for start, end in intervals[1:]:
        if start > nonredundant[-1][1]:
            nonredundant.append([start, end])
        else:
            nonredundant[-1][1] = max(nonredundant[-1][1], end)
    return nonredundant

def primer_overlap(start, end, primer_starts, primer_ends):
    i  = bisect.bisect_right(primer_ends, start)
    overlap = 0
    for j in range(i, len(primer_ends)):
        if end > primer_starts[j]:
            overlap += (min(end, primer_ends[j]) - 
                        max(start, primer_starts[j]))
        else:
            break
    return overlap

def read_skip_intervals(read, min_deletion_length):
    cigar_tuples = read.cigartuples
    if cigar_tuples is None:
        return []
    skip_intervals = []
    reference_pos = read.reference_start
    overhang = 0
    for operation, length in cigar_tuples:
        if (operation in (pysam.CREF_SKIP, pysam.CDEL) and 
            length >= min_deletion_length):
            if skip_intervals:
                skip_intervals[-1][-1] = overhang
            skip_intervals.append([reference_pos, reference_pos + length, overhang, None])
            overhang = 0
        elif operation in BAM_REF_CONSUMING:
            overhang += length
        
        if operation in BAM_REF_CONSUMING:
            reference_pos += length
    
    if skip_intervals:
        skip_intervals[-1][-1] = overhang

    return skip_intervals
