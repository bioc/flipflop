/*
*/

#include "common.h"
#include "structdef.h"

/**
GLOBAL PARAMETERS
**/

/*
verbosity
*/
int n_INST=0;

int VERBOSE=0;

/*
minimum coverage fraction to calculate minimum coverage cutting threshold
*/
double MIN_CVG_FRACTION=0.25;

/* 
minimum exon coverage 
*/
int MINEXONCVG=2; 

/*
if this is true, only exon boundaries from annotations and gene ranges from annotations are used.
*/
bool REFONLY=false; 

/*
this variable indicates whether we use fixed gene range provided by input
*/
bool FIXRANGE=false;

/*
This variable indicates whether we use provided boundary instead of boundaries from the reads 
*/
bool FIXBOUND=false;


/*
 output the coverage of individual bases into instance file.
*/
bool OUTPUT_INDIVIDUAL_COVERAGE=false; 


/*
the default minimum junction (adjustable by -j option)
*/
int DEFAULT_MIN_JUNCTION=1; 


/*
For the gene range creation:
This is the maximum gap allowed for unconnected reads to be considered as one gene.
If the gap between two reads is higher than this number, it will be separated as two gene range.
*/
int MIN_GRANGE_DISTANCE=100;

/*
Minimum number of reads to be considered in one instance
*/
int MIN_GRANGE_READ_CNT=4;

/*
The maximum number of instance to be output. Default -1 (no limit)
*/
int MAX_INSTANCE=-1;



/*
The maximum distance allowed between pairs of paired-end reads.
*/
long MAX_PE_DISTANCE=700000;



/*
The maximum number of segments allowed in an instance.
(DEPRECIATED)
*/
int MAX_N_SEGS=500000;


/*
Treat paired-end reads as single-end reads.
*/
int USE_SINGLE_ONLY=0;

/*
Parameters used to remove reads that spans too long and spans too many reads.
if a read spans >=MAX_READ_SPAN and covers >=MAX_READ_OVER reads, this read is considered as an error read.
*/
int MAX_READ_SPAN=50000;
int MAX_READ_OVER=100000;

/*
Parameters used to remove some incorrectly mapped reads.
When reads are mapped to junctions, this parameter checks the length of the first/last segment of this read. If the length is smaller than this parameter, such segment is discarded.
*/
//int MIN_SEG_FL_OVERLAP=10;
int MIN_SEG_FL_OVERLAP=0;



/* To save the memory usage, replace the read name with some simpler strings */
bool REPLACE_READ_NAME=true;

/* Are the input reads stranded? 1 for positive, -1 for negative, and 0 for unstranded. This parameter is set by -d/--direction option and will be used when writing instances */
int STRANDED_RNASEQ=0;
