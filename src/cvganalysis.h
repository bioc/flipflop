/*
*/

#ifndef CVGANALYSIS_H
#define CVGANALYSIS_H

#include <map>
#include <vector>
#include "structdef.h"

using namespace std;

/** (Replaced by ReadGroup.getConnectedCoverage())
 * Get the coverage from a set of reads, but skipping exons and paired-end read spans are considered covered
 * cvg must be cleaned before calling.
 * Notice: all coordinates in rpoolstart and rpoolend is INCLUSIVE
 */

void getconnectedcoverage(map<long,int>&cvg, vector<vector<long> >&rpoolstart, vector<vector<long> >&rpoolend,int appearpereads);

/**(Replaced by ReadGroup.getCoverage())
 * Get the coverage from a set of reads
 * cvg must be cleaned before calling.
 * Notice: all coordinates in rpoolstart and rpoolend is INCLUSIVE
 */
void getcoverage(map<long,int>&cvg, vector<vector<long> >&rpoolstart, vector<vector<long> >&rpoolend);

/*
Generate breakpoints according to the coverage. 
Break points are generated if the coverage falls below some threshold, 
but if two break points are too close, supress the break point.
*/
int cutcvg(const map<long,int>& cvg, vector<range_t>& cutpoint,
  int threshold, int minadjrange);

/*
Get the statistics of the coverage
May add extra  identical points to cvg
*/
int getCvgStat(map<long,int>&cvg, map<long,int>& allbound,
    map<long, vector<double> >&cvgstat);

/*
 For a given point, get the closest position
 If fail, return -1
*/
long getClosestPoint(long pos,map<long,int>& bound);
  
#endif
