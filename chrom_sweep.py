#!/usr/bin/env python
import sys
from pybedtools import IntervalFile

def overlaps(a, b):
    """
    Return the amount of overlap, in bp
    between a and b.
    If >0, the number of bp of overlap
    If 0,  they are book-ended.
    If <0, the distance in bp between them
    """
    if b is not None:
        return min(a.end, b.end) - max(a.start, b.start)
    else:
        return -1


def after(a, b):
    """
    Is a after (to the right of) b?
    If so, never need to be looked at again.
    """
    return a.start >= b.end


def scan_cache(qy_i, db_cache, hits):
    """
    Scan the cache of "in-play" intervals from the DATABASE
    for overlaps.  If qy_i is "after" (right of)
    a cached intervals, the cached interval can be
    removed from further consideration. We can trust
    this because the files are position sorted.
    """
    if qy_i is None:
        return db_cache
        
    temp_cache = []
    for db_i in db_cache:
        if (qy_i.chrom == db_i.chrom) and not after(qy_i, db_i):
            temp_cache.append(db_i)
            if overlaps(qy_i, db_i) > 0:
                hits.append(db_i)
    return temp_cache


def report_hits(qy_i, hits):
    """
    Report the number of overlaps b/w query and DATABASE
    """
    print str(qy_i) + "\t" + str(len(hits))


def chrom_check(qy_i, db_i, QUERY, DATABASE, db_cache, hits):
    """
    Check if both files are at the same chromosome.
    If not, we need to fast forward the lagging file
    before we proceed.
    """
    if db_i is None or qy_i.chrom == db_i.chrom:
        return (qy_i, db_i, db_cache, hits)
    # ** The Query ** has switched chroms. We must fast-forward B
    if (qy_i.chrom > db_i.chrom):
        tmp_db_i = db_i
        while (tmp_db_i is not None and (tmp_db_i.chrom < qy_i.chrom)):
            tmp_db_i = get_next(DATABASE)
        return (qy_i, tmp_db_i, [], hits)
    # ** The Datasbase ** has switched chroms. We must fast-forward A 
    # and scan each against the database cache in search 
    # of hits from the previous chrom
    elif (qy_i.chrom < db_i.chrom):
        tmp_qy_i = qy_i
        while (tmp_qy_i is not None and (tmp_qy_i.chrom == qy_i.chrom)):
            db_cache = scan_cache(tmp_qy_i, db_cache, hits)
            report_hits(tmp_qy_i, hits)
            tmp_qy_i = get_next(QUERY)
            hits = []
        return (tmp_qy_i, db_i, [], hits)


def get_next(ivls):
    """
    Get the next interval in a file (ivls).
    Return None if at EOF
    """
    try:
        return ivls.next()
    except StopIteration:
        return None


def sweep(QUERY, DATABASE):
    """
    Sweep through QUERY and DB (interval files) in one pass
    and detect overlaps on the fly.

    In BEDTools parlance, QUERY == A, DB == B
    """
    hits  = []
    db_cache = []

    # grab the first interval from each file
    qy_i = get_next(QUERY)
    db_i = get_next(DATABASE)

    while qy_i is not None:

        # Check if we have changed chromosomes. if so, we need to fast-forward
        # the correct chrom, report remining query overlaps, and update the cache
        (qy_i, db_i, db_cache, hits) = chrom_check(qy_i, db_i, QUERY, DATABASE, db_cache, hits)
        
        # Scan the database's of seen, 
        # yet still active feature for overlaps with the current query
        db_cache = scan_cache(qy_i, db_cache, hits)

        # Keep advancing the database until we'e:
        # 1. reached EOF, 2. Changed chromosomes, or
        # 3. Reached an interval that is AFTER  the query (start > query's end)
        # We add each feature to the cache, and track those that overlap
        while (db_i is not None and qy_i.chrom == db_i.chrom and not after(db_i, qy_i)):
            if (overlaps(qy_i, db_i) > 0):
                hits.append(db_i)
            db_cache.append(db_i)
            db_i = get_next(DATABASE)
        
        # Report the query's overlaps and move on to the next query
        report_hits(qy_i, hits)
        hits = []
        qy_i = get_next(QUERY)


if __name__ == "__main__":
    # what are the BED files
    query_file = sys.argv[1]
    database_file = sys.argv[2]

    # open up the BED files.
    QUERY       = IntervalFile(query_file) # The Query File
    DATABASE    = IntervalFile(database_file) # The Database File
    
    sweep(QUERY, DATABASE)

