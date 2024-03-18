import pysam
import sys

def main():
    bam_file = sys.argv[1]
    out_file = sys.argv[2]
    out_QC_file = open(sys.argv[3],"w")
    out_skipped_reads = open(sys.argv[4],"w")

    out_QC_file.write("read_name"+"\t"+"clip_5"+"\t"+"clip_3"+"\n")

    with pysam.AlignmentFile(bam_file, 'rb') as inbam, pysam.AlignmentFile(out_file, 'wb', header=inbam.header) as outbam:
        for i, read in enumerate(inbam):
            print(f"Processing read {i}: {read.query_name}")
            if read.is_unmapped or read.query_sequence is None or read.query_qualities is None:
                print(f"Skipping read {i}: {read.query_name} due to missing information")
                out_skipped_reads.write(str(read.query_name)+"\n")
                continue

            newcigar = []
            clip_5 = 0
            clip_3 = 0

            changed = False
            inseq = False
            for op, length in read.cigar:
                if op == 5:  # H
                    changed = True
                elif op == 4:  # S
                    changed = True
                    if not inseq:
                        clip_5 = length
                    else:
                        clip_3 = length
                else:
                    inseq = True
                    newcigar.append((op, length))

            read.cigar = newcigar
            orig_length = len(read.seq)

            s = read.seq
            q = read.qual

            if clip_3:
                read.seq = s[clip_5:-clip_3]
                if q:
                    read.qual = q[clip_5:-clip_3]
            else:
                read.seq = s[clip_5:]
                if q:
                    read.qual = q[clip_5:]

            newtags = []
            if clip_5:
                newtags.append(('ZA', clip_5))
            if clip_3:
                newtags.append(('ZB', clip_3))

            newtags.append(('ZC', float(clip_5 + clip_3) / orig_length))

            read.tags = read.tags + newtags

            outbam.write(read)

            if clip_5 or clip_3:
                out_QC_file.write(str(read.query_name)+"\t"+str(clip_5)+"\t"+str(clip_3)+"\n")


if __name__ == '__main__':
    main()
