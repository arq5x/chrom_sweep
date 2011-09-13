Development playground for a sweeping line algorithm that detects overlaps between
position sorted BED/GFF/VCF files.

The concept is much like a sort-merge.  One sweeps a imaginary line from start of the chromosome
to the end and reports overlaps between the query file and the database file.  As such, **essentially zero
memory is used in most cases**.

Example (The first file is the "query" and the second is the "database").  The last column is the number of overlaps.

::

	$ cat idx_test.bed
	chr1	10000	10468	trf	789	+
	chr1	10627	10800	trf	346	+
	chr1	10757	10997	trf	434	+
	chr1	11225	11447	trf	273	+
	chr1	11271	11448	trf	187	+
	chr1	11283	11448	trf	199	+
	chr1	19305	19443	trf	242	+
	chr1	20828	20863	trf	70	+
	chr1	30862	30959	trf	79	+
	chr1	44835	44876	trf	73	+
	chr2	10000	10468	trf	789	+
	chr2	10627	10800	trf	346	+
	chr2	10757	10997	trf	434	+
	chr2	11225	11447	trf	273	+
	chr2	11271	11448	trf	187	+
	chr2	11283	11448	trf	199	+
	chr2	19305	19443	trf	242	+
	chr2	20828	20863	trf	70	+
	chr2	30862	30959	trf	79	+
	chr2	44835	44876	trf	73	+

	$ time ./linesweep.py idx_test.bed idx_test.bed
	chr1	10000	10468	trf	789	+	1
	chr1	10627	10800	trf	346	+	2
	chr1	10757	10997	trf	434	+	2
	chr1	11225	11447	trf	273	+	3
	chr1	11271	11448	trf	187	+	3
	chr1	11283	11448	trf	199	+	3
	chr1	19305	19443	trf	242	+	1
	chr1	20828	20863	trf	70	+	1
	chr1	30862	30959	trf	79	+	1
	chr1	44835	44876	trf	73	+	1
	chr2	10000	10468	trf	789	+	1
	chr2	10627	10800	trf	346	+	2
	chr2	10757	10997	trf	434	+	2
	chr2	11225	11447	trf	273	+	3
	chr2	11271	11448	trf	187	+	3
	chr2	11283	11448	trf	199	+	3
	chr2	19305	19443	trf	242	+	1
	chr2	20828	20863	trf	70	+	1
	chr2	30862	30959	trf	79	+	1
	chr2	44835	44876	trf	73	+	1

	real	0m0.218s
	user	0m0.094s
	sys	0m0.113s
	
	$ time intersectBed -a idx_test.bed -b idx_test.bed -c
	chr1	10000	10468	trf	789	+	1
	chr1	10627	10800	trf	346	+	2
	chr1	10757	10997	trf	434	+	2
	chr1	11225	11447	trf	273	+	3
	chr1	11271	11448	trf	187	+	3
	chr1	11283	11448	trf	199	+	3
	chr1	19305	19443	trf	242	+	1
	chr1	20828	20863	trf	70	+	1
	chr1	30862	30959	trf	79	+	1
	chr1	44835	44876	trf	73	+	1
	chr2	10000	10468	trf	789	+	1
	chr2	10627	10800	trf	346	+	2
	chr2	10757	10997	trf	434	+	2
	chr2	11225	11447	trf	273	+	3
	chr2	11271	11448	trf	187	+	3
	chr2	11283	11448	trf	199	+	3
	chr2	19305	19443	trf	242	+	1
	chr2	20828	20863	trf	70	+	1
	chr2	30862	30959	trf	79	+	1
	chr2	44835	44876	trf	73	+	1

	real	0m0.005s
	user	0m0.001s
	sys	0m0.002s