#!/usr/bin/env python
import sys
import os
import index_bed
from pybedtools import IntervalFile

def overlaps(a, b):
    """
Return the amount of overlap, in bp
between a and b.
If >0, the number of bp of overlap
If 0, they are book-ended.
If <0, the distance in bp between them
"""
    return min(a.end, b.end) - max(a.start, b.start)
    
def after(a, b):
    """
Is a after (to the right of) b?
If so, never need to be looked at again.
"""
    return a.start >= b.end
    
def scan_cache(a, cache, hits):
    """
Scan the cache of "in-play" intervals from B
for overlaps. If a is "after" (right of)
a cached intervals, the cached interval can be
removed from further consideration. We can trust
this because the files are position sorted.
"""
    if a is None:
        return cache
        
    temp_cache = []
    for b in cache:
        if not after(a,b):
            temp_cache.append(b)
            if overlaps(a,b) > 0:
                hits.append(b)
    return temp_cache

def report_hits(a, hits):
    """
Report the number of overlaps b/w a and B
"""
    print str(a) + "\t" + str(len(hits))
    
def get_next(ivls):
    """
Get the next interval in a file (ivls).
Return None if at EOF
"""
    try:
        return ivls.next()
    except StopIteration:
        return None

def sweep(A, B, A_card, B_card):
    """
Sweep through A and B (interval files) in one pass
and detect overlaps on the fly.
A is treated as the "query"
B is treated as the "database"
TODO: How to handle chromosome changes?
"""
    hits = []
    B_cache = []

    # grab the first interval from each file
    a = get_next(A)
    b = get_next(B)
    a_cnt = b_cnt = 1
    # loop through each ivl in A (query) in search of overlaps
    while a_cnt <= A_card:
        # scan B's cache for overlaps
        B_cache = scan_cache(a, B_cache, hits)
        # keep consuming B (the database) until no more overlaps w/ a
        while ((b_cnt <= B_card) and (overlaps(a, b) > 0)):
            hits.append(b)
            B_cache.append(b)
            b = get_next(B)
            b_cnt += 1
        # force consumption of next b if was not advanced b/c of no overlaps
        if len(hits) == 0 and b is not None:
            B_cache.append(b)
            b = get_next(B)
            b_cnt += 1
        # report a's overlaps and move on to the next a
        report_hits(a, hits)
        hits = []
        a = get_next(A)
        a_cnt += 1
    

if __name__ == "__main__":
    """
1. If indices don't yet exist for the BED files, create them
2. Use the index to jump to each chrom in the query (A) file
3. Jump to the same chrom in the database (B) file.
4. Sweep through the chrom in the query and database and report hits for A

INDEX entry structure:
0 = chrom,
1 = start_byte,
2 = end_byte,
3 = num_records,
4 = max_interval size (future, for binary search)
"""

    # what are the BED files
    A_file = sys.argv[1]
    B_file = sys.argv[2]

    # expected index file name
    A_idx_file = A_file + ".idx"
    B_idx_file = B_file + ".idx"

    # open up the BED files.
    A = IntervalFile(A_file) # The Query File
    B = IntervalFile(B_file) # The Database File
    
    # create index files if they don't yet exist
    if not os.path.exists(A_idx_file):
        index_bed.index(A_file)
    if not os.path.exists(B_idx_file):
        index_bed.index(B_file)
    
    # load the indices for A and B
    A_map = [] # list of chrom/offset tuples
    for line in open(A_idx_file):
        fields = line.strip().split("\t")
        A_map.append((fields[0], int(fields[1]), int(fields[2]), int(fields[3]), int(fields[4])))
    B_map = [] # list of chrom/offset tuples
    for line in open(B_idx_file):
        fields = line.strip().split("\t")
        B_map.append((fields[0], int(fields[1]), int(fields[2]), int(fields[3]), int(fields[4])))

    # sweep for each chrom in A, the "query file"
    for entry in A_map:
        a_chrom = entry[0]
        a_offset = entry[1]
        b_offset = entry[1]

        # jump to the offset for this chrom
        A.seek(a_offset)
        B.seek(b_offset)
        
        # sweep this chrom in A and B
        A_cardinality = entry[3]
        B_cardinality = entry[3]
        sweep(A, B, A_cardinality, B_cardinality)

