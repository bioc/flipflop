/**
Common variables are defined here.
See COMMON.CPP for the actural definition.
*/

#ifndef COMMON_H
#define COMMON_H
/**
GLOBAL PARAMETERS
**/

/*
verbosity
*/
extern int n_INST;

extern int VERBOSE;

/*
minimum coverage fraction to calculate minimum coverage cutting threshold
*/
extern double MIN_CVG_FRACTION;

/*
 minimum exon coverage 
*/
extern int MINEXONCVG;


/*
if this is true, only exon boundaries from annotations and gene ranges from annotations are used.
*/
extern bool REFONLY; 

/*
this variable indicates whether we use fixed gene range provided by input
*/
extern bool FIXRANGE;

/*
This variable indicates whether we use provided boundary instead of boundaries from the reads 
*/
extern bool FIXBOUND;


/*
 output the coverage of individual bases into instance file.
*/
extern bool OUTPUT_INDIVIDUAL_COVERAGE; 

/*
the default minimum junction (adjustable by -j option)
*/
extern int DEFAULT_MIN_JUNCTION; 


/*
For the gene range creation:
This is the maximum gap allowed for unconnected reads to be considered as one gene.
If the gap between two reads is higher than this number, it will be separated as two gene range.
*/
extern int MIN_GRANGE_DISTANCE;


/*
Minimum number of reads to be considered in one instance
*/
extern int MIN_GRANGE_READ_CNT;

/*
The maximum number of instance to be output. Default -1 (no limit)
*/
extern int MAX_INSTANCE;

/*
The maximum distance allowed between pairs of paired-end reads.
*/
extern long MAX_PE_DISTANCE;


/*
The maximum number of segments allowed in an instance.
(DEPRECIATED)
*/
extern int MAX_N_SEGS;

/*
Treat paired-end reads as single-end reads.
*/
extern int USE_SINGLE_ONLY;



/*
Parameters used to remove reads that spans too long and spans too many reads.
if a read spans >=MAX_READ_SPAN and covers >=MAX_READ_OVER reads, this read is considered as an error read.
*/
extern int MAX_READ_SPAN;
extern int MAX_READ_OVER;

/*
Parameters used to remove some incorrectly mapped reads.
When reads are mapped to junctions, this parameter checks the length of the first/last segment of this read. If the length is smaller than this parameter, such segment is discarded.
*/
extern int MIN_SEG_FL_OVERLAP;

/* To save the memory usage, replace the read name with some simpler strings */
extern bool REPLACE_READ_NAME;

/* Are the input reads stranded? 1 for positive, -1 for negative, and 0 for unstranded. This parameter is set by -d/--direction option and will be used when writing instances */
extern int STRANDED_RNASEQ;

/* 2015-01-15 ELSA: add option. Boundary should be created based on coverage discrepancies or not */
extern int CVG_CUT;
 
#endif
