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


def scan_cache(curr_query, db_cache, hits):
    """
    Scan the cache of "in-play" intervals from the DATABASE
    for overlaps.  If curr_query is "after" (right of)
    a cached intervals, the cached interval can be
    removed from further consideration. We can trust
    this because the files are position sorted.
    """
    if curr_query is None:
        return db_cache
        
    temp_cache = []
    for curr_database in db_cache:
        if (curr_query.chrom == curr_database.chrom) and not after(curr_query, curr_database):
            temp_cache.append(curr_database)
            if overlaps(curr_query, curr_database) > 0:
                hits.append(curr_database)
    return temp_cache


def report_hits(curr_query, hits):
    """
    Report the number of overlaps b/w query and DATABASE
    """
    print str(curr_query) + "\t" + str(len(hits))


def chrom_check(curr_query, curr_database, QUERY, DATABASE, db_cache, hits):
    """
    Check if both files are at the same chromosome.
    If not, we need to fast forward the lagging file
    before we proceed.
    """
    if curr_database is None or curr_query.chrom == curr_database.chrom:
        return (curr_query, curr_database, db_cache, hits)
    # ** The Query ** has switched chroms. We must fast-forward B
    if (curr_query.chrom > curr_database.chrom):
        tmp_curr_database = curr_database
        while (tmp_curr_database is not None and (tmp_curr_database.chrom < curr_query.chrom)):
            tmp_curr_database = get_next(DATABASE)
        return (curr_query, tmp_curr_database, [], hits)
    # ** The Datasbase ** has switched chroms. We must fast-forward A 
    # and scan each against the database cache in search 
    # of hits from the previous chrom
    elif (curr_query.chrom < curr_database.chrom):
        tmp_curr_query = curr_query
        while (tmp_curr_query is not None and (tmp_curr_query.chrom == curr_query.chrom)):
            db_cache = scan_cache(tmp_curr_query, db_cache, hits)
            report_hits(tmp_curr_query, hits)
            tmp_curr_query = get_next(QUERY)
            hits = []
        return (tmp_curr_query, curr_database, [], hits)


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
    curr_query = get_next(QUERY)
    curr_database = get_next(DATABASE)

    while curr_query is not None:

        # Check if we have changed chromosomes. if so, we need to fast-forward
        # the correct chrom, report remining query overlaps, and update the cache
        (curr_query, curr_database, db_cache, hits) = chrom_check(curr_query, curr_database, QUERY, DATABASE, db_cache, hits)
        
        # Scan the database's of seen, 
        # yet still active feature for overlaps with the current query
        db_cache = scan_cache(curr_query, db_cache, hits)

        # Keep advancing the database until we'e:
        # 1. reached EOF, 2. Changed chromosomes, or
        # 3. Reached an interval that is AFTER  the query (start > query's end)
        # We add each feature to the cache, and track those that overlap
        while (curr_database is not None and curr_query.chrom == curr_database.chrom and not after(curr_database, curr_query)):
            if (overlaps(curr_query, curr_database) > 0):
                hits.append(curr_database)
            db_cache.append(curr_database)
            curr_database = get_next(DATABASE)
        
        # Report the query's overlaps and move on to the next query
        report_hits(curr_query, hits)
        hits = []
        curr_query = get_next(QUERY)


if __name__ == "__main__":

    if len(sys.argv) < 3:
        print "Usage:"
        print "chrom_sweep.py [QUERY] [DATABASE]"
        sys.exit()

    query_file = sys.argv[1]
    database_file = sys.argv[2]

    # open up the BED files.
    QUERY       = IntervalFile(query_file) # The Query File
    DATABASE    = IntervalFile(database_file) # The Database File
    
    sweep(QUERY, DATABASE)

