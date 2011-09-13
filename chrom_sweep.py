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


def scan_cache(curr_qy, db_cache, hits):
    """
    Scan the cache of "in-play" intervals from the database
    for overlaps.  If curr_qy is "after" (right of)
    a cached intervals, the cached interval can be
    removed from further consideration. We can trust
    this because the files are position sorted.
    """
    if curr_qy is None:
        return db_cache
        
    temp_cache = []
    for curr_db in db_cache:
        if (curr_qy.chrom == curr_db.chrom) and not after(curr_qy, curr_db):
            temp_cache.append(curr_db)
            if overlaps(curr_qy, curr_db) > 0:
                hits.append(curr_db)
    return temp_cache


def report_hits(curr_qy, hits):
    """
    Report the number of overlaps b/w query and database
    """
    print str(curr_qy) + "\t" + str(len(hits))


def chrom_check(curr_qy, curr_db, query, database, db_cache, hits):
    """
    Check if both files are at the same chromosome.
    If not, we need to fast forward the lagging file
    before we proceed.
    """
    if curr_db is None or curr_qy.chrom == curr_db.chrom:
        return (curr_qy, curr_db, db_cache, hits)
    # ** The Query ** has switched chroms. We must fast-forward B
    if (curr_qy.chrom > curr_db.chrom):
        tmp_curr_db = curr_db
        while (tmp_curr_db is not None and (tmp_curr_db.chrom < curr_qy.chrom)):
            tmp_curr_db = get_next(database)
        return (curr_qy, tmp_curr_db, [], hits)
    # ** The Database ** has switched chroms. We must fast-forward A 
    # and scan each against the database cache in search 
    # of hits from the previous chrom
    elif (curr_qy.chrom < curr_db.chrom):
        tmp_curr_qy = curr_qy
        while (tmp_curr_qy is not None and (tmp_curr_qy.chrom == curr_qy.chrom)):
            db_cache = scan_cache(tmp_curr_qy, db_cache, hits)
            report_hits(tmp_curr_qy, hits)
            tmp_curr_qy = get_next(query)
            hits = []
        # catch query up to database
        while (tmp_curr_qy is not None and (tmp_curr_qy.chrom < curr_db.chrom)):
            # hits is empty to reflect the fact that no hits are found in catch-up mode
            report_hits(tmp_curr_qy, hits)
            tmp_curr_qy = get_next(query)
        return (tmp_curr_qy, curr_db, [], hits)


def get_next(ivls):
    """
    Get the next interval in a file (ivls).
    Return None if at EOF
    """
    try:
        return ivls.next()
    except StopIteration:
        return None


def sweep(query, database):
    """
    Sweep through query and DB (interval files) in one pass
    and detect overlaps on the fly.

    In BEDTools parlance, query == A, DB == B
    """
    hits  = []
    db_cache = []

    # grab the first interval from each file
    curr_qy = get_next(query)
    curr_db = get_next(database)

    while curr_qy is not None:

        # Check if we have changed chromosomes. if so, we need to fast-forward
        # the correct chrom, report remining query overlaps, and update the cache
        (curr_qy, curr_db, db_cache, hits) = chrom_check(curr_qy, curr_db, query, database, db_cache, hits)
        
        # Scan the database's of seen, 
        # yet still active feature for overlaps with the current query
        db_cache = scan_cache(curr_qy, db_cache, hits)

        # Keep advancing the database until we'e:
        # 1. reached EOF, 2. Changed chromosomes, or
        # 3. Reached an interval that is AFTER  the query (start > query's end)
        # We add each feature to the cache, and track those that overlap
        while (curr_db is not None and curr_qy.chrom == curr_db.chrom and not after(curr_db, curr_qy)):
            if (overlaps(curr_qy, curr_db) > 0):
                hits.append(curr_db)
            db_cache.append(curr_db)
            curr_db = get_next(database)
        
        # Report the query's overlaps and move on to the next query
        report_hits(curr_qy, hits)
        hits = []
        curr_qy = get_next(query)


if __name__ == "__main__":

    if len(sys.argv) < 3:
        print "Usage:"
        print "chrom_sweep.py [query] [database]"
        sys.exit()

    query_file = sys.argv[1]
    database_file = sys.argv[2]

    # open up the BED files.
    query       = IntervalFile(query_file) # The Query File
    database    = IntervalFile(database_file) # The Database File
    
    sweep(query, database)

