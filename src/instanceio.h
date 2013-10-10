#ifndef INSTANCEIO_H
#define INSTANCEIO_H

#include <map>
#include <vector>
#include "common.h"
#include "structdef.h"

using namespace std;
/**
 * This is used to get all boundaries of given reads
 * minjunccount: minimum junction count
 * In all bound, the boundary will be given. There are 3 types of the boundary: junction boundary (marked by 0) and low coverage boundary (marked by 1), and gene annotation boundary (marked by 2)
 */
void getallbound(map<long,int>& allbound, range_t currentrange,
  vector<vector<long> >& rpoolstart, vector<vector<long> >&rpoolend,
  map<long, vector<double> >& cvgstat, //save the statistics of each segment, including 0:max cvg 1:leftmost cvg 2:rightmost cvg 3: fraction of zero coverage
  bool refonly
);


/**
 * Given a set of boundaries, group reads into types  and orders
 */
void getalltypeorders(map<vector<int>, int> &contenttype, map<vector<int>,int>& typeorder,
		vector<vector<int> >&seriestype, vector<int>& allreadtype,
		map<long,int>& allbound, vector<vector<long> >&rpoolstart, vector<vector<long> >&rpoolend,
		map<long,int>& njuncs // this map stores, for each segment, how many junction reads fall onto that segment
	);



/**
 * Remove single-end or paired-end read which are outside the given range
 */
void removeoutofrangereads(range_t & currentrange,
		vector<vector<long> >&rpoolstart, vector<vector<long> >&rpoolend,
		int appearpereads);



#endif
