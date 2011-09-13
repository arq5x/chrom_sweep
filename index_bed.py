#!/usr/bin/env python
import sys
from pybedtools import BedTool, IntervalFile
        
def index(ivls):
    """

    """
    chroms = []
    chrom_start_offsets = {}
    chrom_end_offsets = {}
    chrom_num_lines = {}
    chrom_max_interval = {}

    byte_offset = 0
    prev_chrom  = None    
    for line in open(ivls):
        line_bytes = len(line)
        fields = line.strip().split("\t")
        curr_chrom = fields[0]
        # we've had a chrom change!
        # record the start byte for the new chrom and
        # the end byte for the previous chrom
        if curr_chrom != prev_chrom:
            chroms.append(curr_chrom)
            if prev_chrom is not None:
                chrom_end_offsets[prev_chrom] = byte_offset
            chrom_start_offsets[curr_chrom] = byte_offset
            # initialize the current chrom's info
            chrom_num_lines[curr_chrom] = 1
            chrom_max_interval[curr_chrom] = int(fields[2]) - int(fields[1])
            # update our chroms
            prev_chrom = curr_chrom
        # same chromosome block
        else:
            # increment the total lines for this chrom
            # and update the longest interval seen si necessaire.
            chrom_num_lines[prev_chrom] += 1
            curr_ivl_length = int(fields[2]) - int(fields[1])
            if curr_ivl_length > chrom_max_interval[prev_chrom]:
                chrom_max_interval[prev_chrom] = curr_ivl_length
        # we've consumed more bytes
        byte_offset += line_bytes
    # handle last chrom
    chrom_end_offsets[curr_chrom] = byte_offset
    
    idx_filename = ivls + ".idx"
    idx          = open(idx_filename, 'w')
    for chrom in chroms:
        idx.write("\t".join([chrom, str(chrom_start_offsets[chrom]), \
                                    str(chrom_end_offsets[chrom]), \
                                    str(chrom_num_lines[chrom]), \
                                    str(chrom_max_interval[chrom])]))
        idx.write("\n")
    idx.close()

if __name__ == "__main__":
    bed_file = sys.argv[1]
    index(bed_file)

