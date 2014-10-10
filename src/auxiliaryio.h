/**
* Auxiliary file IO
*/
#ifndef AUXILIARYIO_H
#define AUXILIARYIO_H

#include <vector>
#include <string>
#include <map>
#include "common.h"
#include "structdef.h"

using namespace std;

/*
Prepare the aux output file, including parsing the input arguments and open corresponding files
*/
void prepareAuxFile(vector<string> args,
	string prefix,
	bool& outinstance
	);



void closeAuxFile();

/*
ISOINFER OUTPUT files
*/

void write2grange(int rangecounter, range_t& nowrange, string rname,int duprcount);
/*
Write one line in readinfo file.
Used for readSamFile().
*/
void writeoneline2readinfo( string rname, int flag, vector<long>& startpos, vector<long>& endpos,int pos, bool writeendl);
	
int writeoneline2readinfo(Align &al);

int readgenerangefile(string filename, GeneRange& gr);

/**
ISOLASSO aux files
*/

/**
 * Write annotation files (wig: coverage, bedfs: junctions)
 */
void writeannotation(
  string chrname,
  vector<vector<long> >& rpoolstart, vector<vector<long> >&rpoolend,
  int appearpereads,
  int instanceid);
	
	
void writeannotation(
    ReadGroup & rg,
    int instanceid);

void writeCoverage( ReadGroup & rg);
//write to boundfs
//range is inclusive
void writeboundfs(range_t range,int instanceid,string chrname,
		vector<range_t>& allrange );
	

/*
Write to file
*/
int write2rangeandinstance(ReadGroup& rd,range_t& range);

/**
 * Write the constructed instance to output
 * All coordinates should be 1 base.
 * Notice that rpoolstart and rpoolend are inclusive
 */
void write2Instance( range_t& currentrange,
	       	vector<vector<long> >&rpoolstart, vector<vector<long> >&rpoolend,
		int appearpereads,
		int instanceid,
		string chrname,
		int readlen,
		vector<vector<long> >& annostart, vector<vector<long> >& annoend,   //annotations within currentrange
		bool wantcoverage,
		bool refonly
		);

/*read existing annotation file (using .bed format), and incorporate into instance, and cluster them at the same time
  Notice that the start coordinate of bed file is 0-based, and the end coordinate is 1-based.
  so to be compatible with SAM and our program (SAM  uses 1-based, and the range in this program is 1-based, and inclusive [start,end]),
  we need to handle the coordinate carefully
*/
int parseannobed(string infile, 
Annotation& anno
	);

#endif
