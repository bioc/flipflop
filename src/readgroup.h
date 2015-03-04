/*
Define a cluster of reads
*/

#ifndef READGROUP_H
#define READGROUP_H

#include "align.h"
#include "rangeset.h"
#include "commontype.h"
#include <iostream>
#include "samples.h" // 2015-01-15 ELSA


/*
Read group, a set of reads, including the start and end positions
*/

typedef vector<int> type_t;

class ReadGroup{
private:

  const Samples* samples; // 2015-01-15 ELSA

  //vector<Align> allalign;
  /*
  Data: a group of start/end positions, and its directions
  */
  vpos_t start;
  vpos_t end;
  vi direction;

  /* This is the ID of the read. 
     ATTENTION: I do not store read ID here; the cost is too expensive. And the IDs will not split when you call splitByRangeSet() function. Use it at your cost.
  */
  vector<string> readid;

  /*a record of paired-end reads. The same size as start, indicating the index of pairs
  */
  vi pairs;
  /*
  During the search of the pair, store the unpaired position, read names and their index. -1 for non-paired-end reads.
  */
  //map<string,range_t > topair;
  map<long, map<string,long> >topair;
  
  /* This variable forces the paired-end reads to be treated as single-end reads. */
  bool forcesingle;
  
  /* boundary */
  map<long,int> allbound;
  /* Segments */
  vector<range_t> segs;
  /* Segment read counts */
  vector<int> segsCnt;
  /* Segment junction read counts */
  vector<int> segsJCnt;

  /* indicator for valid reads */
  vector<int> validRead;
  /* indicator for valid segments */
  vector<int> validSeg;
  /* indicator for valid types*/
  vector<int> validSGType;

  /*allbound -> segs */
  void bound2Seg();

  /*save the statistics of each segment, including 0:max cvg 1:leftmost cvg 2:rightmost cvg 3: fraction of zero coverage */
  map<long, vector<double> > cvgstat; 

  /* Segmentation reads and types
  */
  vector<type_t> alltype;

  vector<string> rdg_names; // 2015-01-15 ELSA
  map<string, vector<int> > rdg_typecount;

  /* THe number of occurance of each type on reads */
  vector<int> typecount;
  /* The direction of the type. >0/<0/=0 */
  vector<int> typedir;
  /* Recording which type a read belongs to. -1 leads to invalid type */
  vector<int> read2type;
  type_t getType(int n) const; //get the type
  type_t getType(const pos_t& sn,const pos_t& en) const;

  
  /* the set of ranges of the reads in this group */
  RangeSet rg;

  bool isfixrange;
  range_t fixrange;
  
  void setRangeSet(RangeSet& s);
  /* Calculate the range set of this read group */
  void calculateRangeSet();
  
  /* Return the index of reads containing specified segments */
  vector<int> getReadsContainingSegs(int segid);

  static int statReadLen;

  static string statChr;

  /* Orientation of the instance: 1 for +, -1 for -, 0 for unknown */
  int orientation;
  map<long, vector<double> >::const_iterator getSegStats(map<long,int>&cvg, const range_t& seg);
  //void displayAllBound(const std::string& msg) const;
  /*
  Add a paired-end reads. Do not go through the paired-end search process
  */
  int addPair(pos_t &s1, pos_t &e1,int d1, const string & pairname1 , pos_t&s2, pos_t & e2,int d2, const string & pairname2); // 2015-01-15 ELSA
  /* Add one alignment to the current cluster of reads. Do not try to pair it. */
  int addOnly(pos_t &s, pos_t &e,int dir, const string& rgname); // 2015-01-15 ELSA

public:
  ReadGroup():forcesingle(false),isfixrange(false),orientation(0),samples(NULL){} // 2015-01-15 ELSA
  void setSamples(const Samples* samples);

  void setForceSingle(){forcesingle=true;}
  /*
  Access functions
  */
  vpos_t & s(){return start;}
  vpos_t & e(){return end;}
  /* Get the direction of the READS, not the direction of the instance.*/
  vi& getDirection(){return direction;}
  vector<string> & getID(){return readid;}
  
  //vector<Align>& a(){return allalign;}
  
  int getReadLen(int n=0) const;
  
  string getChr() const;
 
  /* Total number of reads */
  int size()const{return (int)start.size();}
  int validSize() const;
  
  /* Total number of paired end reads */
  int peSize()const;
  
  /* clear all the vectors */
  void clear();
  /* clear paired-end search information (topair) */
  void clearPairInfo();
  /*
  Should rewrite this function to 
   change different behaviors of single/pair reads
  If you want to break paired-end reads into single-end reads, set forcesingle=true
  */
  int add(Align & al);

  bool hasPairEnd(){for(int i=0;i<pairs.size();i++) if(pairs[i]!=-1)return true; return false; }

  /* Return the coverage information
     In connected coverage, positions between junctions and pairs of paired-end reads are all counted.
     In point coverage, only the first base is covered (ignoring paired-end type)
     getCvgStatistics() calculates the statistics of each segments according to current segment.
  */
  void getConnectedCoverage(map<long,int>& cvg);
  void getCoverage(map<long,int>& cvg);
  void getCoverage(map<long,int>&cvg, const range_t & crange);
  void getPointCoverage(map<long,int>& cvg);
  void getCvgStatistics();

  /*
  Get the current range
  */
  range_t getRange() const;

  /* In cases of fix range, set up the range of this read group */
  void setRange(range_t s);

  /*
  Split the read group by range sets, as long as the distances between these rages are greater than some threshold
  */
  void splitByRangeSet(vector<ReadGroup> & vr,long minD);
   
  void splitByRangeSet(vector<ReadGroup> & vr,vector<RangeSet> & vs);

  /* Split the read group by different directions. Useful when two opposite genes overlap. */
  void splitByDirection(vector<ReadGroup>& vr,vector<range_t>& vrt );
/**
 * Get all boundaries of given reads
 * minjunccount: minimum junction count
 * In all bound, the boundary will be given. There are 3 types of the boundary: junction boundary (marked by 0) and low coverage boundary (marked by 1), and gene annotation boundary (marked by 2)
 */
  void calculateBound( 
    bool cvgcut, //if new boundary should be generated from coverage cutoff
    int minjunc,
    int minadjrange
  );
  
/*
  Use provided boundary, instead of reads boundary
*/
  void setupBound(vector<range_t>& jbound);
/*
Calculate the types according to the boundary
*/
  void calculateType();
 
/*
  Return segments with at least one read covered
*/
  vector<range_t> getSegs() const;

  void toStream(ostream &out);

  /* If a segment is too weak, mark it as invalid */
  void removeWeakSegs(float minf);
  /*
  Keep only valid segments, types, and reads
  */
  void calculateValidSegs();

  friend ostream& operator<<(ostream& out, ReadGroup & rg);
  
  /* remove reads spanning too many reads 
     If two segments of a read contains >=minspan bases and contains >maxread reads, this read is marked invalid
  */
  void removeTooLongReads(long minspan,int maxread);

  /* The following functions set up the static read length and chromosome name.
     This is useful when some instances have no mapped reads thus are not able to provide chromosome names and read length
  */
  static void setStatReadLen(int n){statReadLen=n;}
  static void setStatChr(string c){statChr=c;}
  /* 
  get the direction information. If all splicing reads are positive/negative, return +/-100.
  Otherwise, return the ratio between positive/negative. The abs value always >1.
  */
  int getDirSum();

  int getDir() const {return orientation;}
  int setDir(int d){int d0=orientation; orientation=d;return d0;}

};




#endif
