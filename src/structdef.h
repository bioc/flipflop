/**
structdef.h defines commonly used structures.
*/
#ifndef STRUCTDEF_H
#define STRUCTDEF_H

#include <vector>
#include <string>
#include <map>
#include "align.h"
#include "readgroup.h"

using namespace std;



/**
Gene range structure
Storing the set of ranges in the chromosome, including 
chromosome name, range start and end
*/
class GeneRange{
protected:
  int index;
  map<string,int> chrexistmap;
  vector<string> chrmap;
  vector<int> chr;
  vector<range_t > range;

  /*For each chromosome, keep a map structure showing the starting position and the first index of a range */
  vector<map<long,int> > rangepos;
  
public:
  /* BASIC FUNCTIONS */
  void clear(){index=0;chr.clear();range.clear();}
  int inc(){index++;return index;}
  bool isEnd(){return index>=range.size();}
  int size(){return range.size();}
  int getIndex()const{return index;}
  void reset(){index=0;}

  // "Reset everything"
  void resetAll( void ) {
     clear();
     chrexistmap.clear();
     chrmap.clear();
     rangepos.clear();
  }

  /*
  Add one record. Return value is the index of the inserted pairs
  ATTENTION: MUST ensure that the insert order is sorted, or apply sort() after finish insertion.
  */
  int push_back(string c, range_t r);

  /*
  Sort the range according to the chromosome and range
  */
  int sort();

  /* Check if the ranges are sorted according to chromosome name and starting position */
  bool check();
  /*
  Sequential function
  */
  /* Access function */
  string getChr(){return chrmap[chr[index]];}
  range_t getRange(){return range[index];}
  /* Move the index to the next chromosome. Return -1 if reaches the end of the queue */
  int moveToNextChr(string nextChr);
  /*
  Retrieve the next range, and save to range_t.
  If reaches the end of the queue, set cr.first=cr.second=-1;
  if reaches the end of the chromosome, set to 0.
  */
  int getNext(string rname,long pos, string& cname, range_t & cr);

  /* Retrieve all ranges which are inside a given range t*/
  int getAllRanges(string rname,range_t t, vector<range_t>& allr);

};
/*
Annotation, a set of known isoforms
*/

class Annotation{
  /*
  New data structure
  The first map stores the <chromosome, range> relationship, and the second stores the <range, isoform> relationship.
  */
  map<string, map<range_t, ReadGroup> > data;
  
  
  /*
  variables tracking the previous clustered records
  */
  string prevs; //previous chromosome
  range_t prevr; //previous range
  ReadGroup prevg; //previous read group
  
  static ReadGroup EmptyReadGroup; // an empty read group

public:
  Annotation(){
    prevr=range_t(-1,-1);
    prevs="";
  }

  void resetAll( void ) {
     data.clear();
     prevr = range_t( -1, -1 );
     prevs = "";
     EmptyReadGroup.clear();
  }

  /*
  Functions
  */
  /*
  Add one record (rs,re) to the pool. 
  Dir is the direction: +1 is positive, and -1 is negative
  This function will automatically cluster the records, so the records must be added in ascending order by chr and 1st element of rs.
  Return 0 for success, -1 for failure (if the records are not added in order).
  */
  int add(string chr, pos_t rs, pos_t re,int dir,string id);

  /*
  Put current clusters into the data. Only explicitly called at the end of the annotation input.
  */
  int cluster();
  /*
  get the number of chromosomes
  */ 
  int getNChrom();
  /*
  get the list of chromosomes
  */
  vector<string> getChrom();
  /*
  get the number of clusters
  */
  int getNCluster();

  /*
  get a vector of clusters, specified by chromosome name
  */
  vector<range_t> getCluster(string chr);

  /*
  Get the corresponding read group from the chromosome name and the range
  */
  ReadGroup getReadGroup(string chr, range_t r);

  /*
  get the number of annotations
  */
  int getNAnnotation();


  /*
  For a given range (current), check if there are overlapping annotations. If no, do nothing; otherwise, extend the range to include annotations.
  Return value:
    0 no range overlaps
    1 updated range
  */
  int checkOverlapRange(range_t& current, string rname);
};

/*
Global structure storing the gene range
*/
extern GeneRange C_GENE_RANGE;

extern Annotation C_RANGE_ANNO;
/*
Global structure storing all the annotation 
*/
extern Annotation C_ANNO;

/*
Global structure storing the junction
*/
extern GeneRange C_JUNCTION;

#endif
