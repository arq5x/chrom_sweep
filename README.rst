Motivation
===========
Files of genomic features continue to grow in size.  Traditional methods of detecting overlaps (e.g., R-trees, interval trees, binary search) among genomic intervals require substantial memory for very large datasets.  However, if datasets are sorted by chromosome and start position, one can, in principle, adapt principles from sort-merge and sweep-line algorithms to detect overlaps "on the fly" by marching from "left to right" across a chromosome and catching overlaps as you encounter features.

This is the principle behind chromsweep.  Corner cases arise when handling changes in chromosomes in the two files. This version now handles these cases, save for the caveats laid out below in the "issues" section.

Previous work
=============
Joel Richardson published a similar algorithm called "fjoin".  This approach required that the intervals be sorted by position, irrespctive of chromosome.  Thus, it did not report results in chrom order and required files to be sorted in an atypical way for most datasets.  The following implementation has several advantages:

#. It works with files, streams, and FIFOs
#. It works with BED, GFF, VCF.
#. In principle, it works with BAM;  once this is ported to BEDTools it will.
#. Typical annotation and dataset files are already sorted in exactly the manner required.

Issues
======
#. Requires sorting files with sort -k1,1 -k2,2n (BED files). That is, downloading straight from UCSC is not in an order that works with the chrom_check() function.  Need to investigate this.  Perhaps a misunderstanding of the Python "<" operator.
#. Not sure if it handles cases where one file has entries for a given chrom while the other doesn't.
#. Currently, it only reports the _count_ of overlaps between the Query file (A file in BEDTools parlance) and the Database file (B file).


Future Works
==========
#. N-files.
#. Paralellism.



Example
==========
::

	$  ./chrom_sweep.py 
	Usage:
	chrom_sweep.py [QUERY] [DATABASE]

	$ zcat rmsk.bed.gz | head
	chr1	10000	10468	(CCCTAA)n	1504	+
	chr1	10468	11447	TAR1	3612	-
	chr1	11503	11675	L1MC	437	-
	chr1	11677	11780	MER5B	239	-
	chr1	15264	15355	MIR3	318	-
	chr1	16712	16749	(TGG)n	203	+
	chr1	18906	19048	L2a	239	+
	chr1	19947	20405	L3	652	+
	chr1	20530	20679	Plat_L3	270	+
	chr1	21948	22075	MLT1K	254	+
	
	$ zcat knownGene.bed.gz | head
	chr1	11873	14409
	chr1	11873	14409
	chr1	11873	14409
	chr1	14362	16765
	chr1	14362	19759
	chr1	14362	19759
	chr1	14362	19759
	chr1	14362	24901
	chr1	14362	29370
	chr1	14362	29370
	
	$  ./chrom_sweep.py knownGene.bed.gz rmsk.bed.gz | head
	chr1	11873	14409	0
	chr1	11873	14409	0
	chr1	11873	14409	0
	chr1	14362	16765	2
	chr1	14362	19759	3
	chr1	14362	19759	3
	chr1	14362	19759	3
	chr1	14362	24901	10
	chr1	14362	29370	17
	chr1	14362	29370	17
	
	$ ./chrom_sweep.py knownGene.bed.gz rmsk.bed.gz > chrom_sweep.out

	$ intersectBed -a knownGene.bed.gz -b rmsk.bed.gz -c > bedtools.out
	
	$ diff chrom_sweep.out bedtools.out
	