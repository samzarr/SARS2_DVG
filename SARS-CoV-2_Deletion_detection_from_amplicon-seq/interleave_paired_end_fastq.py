import sys
import argparse
import fastq

def main():
    # Create argument parser
    parser = argparse.ArgumentParser(description='Interleave two fastq files.')
    
    # Add arguments
    parser.add_argument('read1_fastq', type=str, help='Read1 input fastq file')
    parser.add_argument('read2_fastq', type=str, help='Read2 input fastq file')
    
    # Parse arguments
    args = parser.parse_args()

    with open(args.read1_fastq, 'r') as f1, open(args.read2_fastq, 'r') as f2:
        interleave_fastq(f1, f2, sys.stdout)

def interleave_fastq(file1, file2, output_file):
    reader1 = fastq.Reader(file1)
    reader2 = fastq.Reader(file2)
    writer = fastq.Writer(output_file)

    for rec1, rec2 in zip(reader1, reader2):
        assert(rec1.name == rec2.name)
        rec1.name = rec1.name + "_1"
        rec2.name = rec2.name + "_2"
        writer.write(rec1)
        writer.write(rec2)

if __name__ == '__main__':
    main()