#include <iostream>
#include "instanceio.h"
#include "cvganalysis.h"
#include "common.h"
#include "structdef.h"
#include "rangeset.h"

/**
 * Replaced by corresponding member functions in ReadGroup
 * This is used to get all boundaries of given reads
 * minjunccount: minimum junction count
 * In all bound, the boundary will be given. There are 3 types of the boundary: junction boundary (marked by 0) and low coverage boundary (marked by 1), and gene annotation boundary (marked by 2)
 */
void getallbound(map<long,int>& allbound, range_t currentrange,
    vpos_t& rpoolstart, vpos_t& rpoolend,
    map<long, vector<double> >& cvgstat, //save the statistics of each segment, including 0:max cvg 1:leftmost cvg 2:rightmost cvg 3: fraction of zero coverage
    bool refonly
    ){
  //this variable is used to set if we need to insert boundary for zero coverage points
  bool considerranges=true;
  map<long,int> boundfilter;
  boundfilter[currentrange.first]=DEFAULT_MIN_JUNCTION;
  boundfilter[currentrange.second+1]=DEFAULT_MIN_JUNCTION;
  allbound[currentrange.first]=0;allbound[currentrange.second+1]=0;
  for(int i=0;i<rpoolstart.size();i++){
    if(rpoolstart[i].size()>1){
      for(int j=1;j<rpoolstart[i].size();j++) boundfilter[rpoolstart[i][j]]++;
      for(int j=0;j<rpoolend[i].size()-1;j++) boundfilter[rpoolend[i][j]+1]++;
    }
  }
  //to ensure the number of segments is smaller than maxnsegs:
  //sort the boundary in a descending order, choose a threshold such that 
  //the number of segments does not exceed maxnsegs
  //
  map<long,int>::iterator mitr,mitr2;
  map<long,int>cvg;
  //get the real coverage
  getcoverage(cvg,rpoolstart,rpoolend);

  //if(maxnsegs==-1|| boundfilter.size()<maxnsegs){
    //if maxnsegs is not set, use DEFAULT_MIN_JUNCTION as a filter
  for(mitr=boundfilter.begin();mitr!=boundfilter.end();mitr++){
    if(mitr->second>=DEFAULT_MIN_JUNCTION && !refonly)
      allbound[mitr->first]=0;//set up boundary. ok to reset 2 to 0
  }
    //now, check the coverage (real), insert boundaries if the coverage drops below 0.
    //but do not insert boundaries if the maximum coverage is too low, or if the gap is too short
    if(considerranges==false)return;
    
    vector<range_t> cutpoint;
    int minadjrange=100;//this number should depend on paired-end insert size
    cutcvg(cvg,cutpoint,0,minadjrange);
    
    // now, check adjacent ranges,
    // if they're too close, merge them
    if(cutpoint.size()>0 && !refonly){
      for(int i=0;i<cutpoint.size();i++){
          if(allbound.count(cutpoint[i].first)==0 )
            allbound[cutpoint[i].first]=1; //begin
          if(allbound.count(cutpoint[i].second)==0 )
            allbound[cutpoint[i].second]=1; //begin
      }
    }
    getCvgStat(cvg,allbound,cvgstat);

  //}

}

/**
 * Given a set of boundaries, group reads into types  and orders
 */
void getalltypeorders(map<vector<int>, int> &contenttype, map<vector<int>,int>& typeorder,
    vector<vector<int> >&seriestype, vector<int>& allreadtype,
    map<long,int>& allbound, vector<vector<long> >&rpoolstart, vector<vector<long> >&rpoolend,
    map<long,int>& njuncs // this map stores, for each segment, how many junction reads fall onto that segment
  ){

  //current order, begin with 0
  int currentorder=0;
  njuncs.clear();

  //current type
  vector<int> currenttype;
  map<long,int>::iterator it, it2;

  for(int i=0;i<rpoolstart.size();i++){
    //go over the range
    currenttype.clear();
    for(int j=0;j<rpoolstart[i].size();j++){
      long a=rpoolstart[i][j],b=rpoolend[i][j];
      //in the normal case, we need to check if range boundary is in the bound too.
      //but in some cases, the spliced read may be contained within one region. 
      //in this case, no boundary check is needed.
      if( (j>0 && allbound.count(a)==0) || (j<rpoolstart[i].size()-1 && allbound.count(b+1)==0)){
        //the boundary is not in the list: remove this read
        //or, do nothing
        //currenttype.clear();
        //break;
      }
      it=allbound.begin();
             it2=it;  it2++;
      int exonid=0;
      while(it2!=allbound.end()){
        if(checkoverlap(a,b,it->first,it2->first-1)) 
        {
          it->second++;
          //in cases of spliced read inside one region, only retain one record
          if(currenttype.size()==0 || currenttype.back()!=exonid)
            currenttype.push_back(exonid);
          if(rpoolstart[i].size()>1){
            //this is a junction read
            njuncs[it->first]++;
          }
        }
        it++;
        it2++;
        exonid++;
      }
    }//end j

    //for single-end reads, update the type count
    contenttype[currenttype]++;
    //for paired end reads, we need to assign a type number to it.
    int currenttypeorder=-1;
    if(typeorder.count(currenttype)==0){
      typeorder[currenttype]=currentorder;
      seriestype.push_back(currenttype);
      currentorder++;
    }
    
    //save this type
    allreadtype.push_back(typeorder[currenttype]);

  }//end i


}


/**
 * Remove single-end or paired-end read which are outside the given range
 */
void removeoutofrangereads(range_t & currentrange,
    vector<vector<long> >&rpoolstart, vector<vector<long> >&rpoolend,
    int appearpereads){

  //check the range:
  //if any of the reads fall out of the range, discard it
  vector<int> validread(rpoolstart.size(),1);
  for(int i=0;i<rpoolstart.size();i++){
    if(rpoolstart[i][0]<currentrange.first)validread[i]=0;
    if(rpoolend[i].back()>currentrange.second) validread[i]=0;
  }
  //make backup, delete invalid reads
  vector<vector<long> > bkpoolstart, bkpoolend;
  if(appearpereads==0){
    for(int i=0;i<rpoolstart.size();i++) 
      if (validread[i]==1){
        bkpoolstart.push_back(rpoolstart[i]);
        bkpoolend.push_back(rpoolend[i]);
      }
  }
  else{
    for(int i=0;i<rpoolstart.size();i+=2){
      if(validread[i]==1 && validread[i+1]==1){
        bkpoolstart.push_back(rpoolstart[i]);
        bkpoolstart.push_back(rpoolstart[i+1]);
        bkpoolend.push_back(rpoolend[i]);
        bkpoolend.push_back(rpoolend[i+1]);
      }
    }
  }
  if(bkpoolstart.size()==0){
    //cout<<"Instance "<<instanceid<<" has 0 reads. before has "<<rpoolstart.size()<<" reads.\n";
    //cout<<"range:"<<currentrange.first<<","<<currentrange.second<<endl;
    //cout<<"reads:"<<endl;
    //for(int i=0;i<rpoolstart.size();i++){
    //  for(int j=0;j<rpoolstart[i].size();j++)cout<<rpoolstart[i][j]<<",";
    //  cout<<"...";
    //  for(int j=0;j<rpoolend[i].size();j++)cout<<rpoolend[i][j]<<",";
    //  cout<<endl;
    }

  

  rpoolstart=bkpoolstart;
  rpoolend=bkpoolend;


}


