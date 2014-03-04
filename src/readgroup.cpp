#include <climits>
#include <sstream>
#include "readgroup.h"
#include "rangeset.h"
#include "cvganalysis.h"
#include "common.h"

int ReadGroup::statReadLen=0;

string ReadGroup::statChr="";

void ReadGroup::clear()
{
  //allalign.clear();
  start.clear();
  end.clear();
  direction.clear();
  pairs.clear();
  topair.clear();
  allbound.clear();
  segs.clear();
  segsCnt.clear();
  segsJCnt.clear();
  
  validRead.clear();
  validSeg.clear();
  validSGType.clear();
  cvgstat.clear();
 
  alltype.clear();
  typecount.clear();
  typedir.clear();
  read2type.clear();

  rg.clear();
  clearPairInfo();

  readid.clear();
}

void ReadGroup::clearPairInfo(){
  topair.clear();
}

/*
Add a paired-end reads. Do not go through the paired-end search process
*/
int ReadGroup::addPair(pos_t &s1, pos_t &e1, int d1, pos_t&s2, pos_t & e2, int d2 ){
  addOnly(s1,e1,d1);
  addOnly(s2,e2,d2);
  pairs[pairs.size()-1]=pairs.size()-2;
  pairs[pairs.size()-2]=pairs.size()-1;
  //allalign[allalign.size()-1].qname=allalign[allalign.size()-2].qname;
  return 0;
}

/* Add one alignment to the current cluster of reads. Do not try to pair it. */
int ReadGroup::addOnly(pos_t &s, pos_t &e,int dir){
  start.push_back(s);
  end.push_back(e);
  direction.push_back(dir);
  //allalign.push_back(al);
  validRead.push_back(1);
  pairs.push_back(-1);
  stringstream ss; ss<<(long)start.size();
  //if(REPLACE_READ_NAME){allalign.back().qname=""; ss>>allalign.back().qname; }
  return 0;
}

/*
Add one alignment to the current cluster of reads. If the read is paired-end reads, try to automatically pair them.
*/
int ReadGroup::add(Align & al){
  addOnly(al.s(),al.e(),al.splicedir);
  if(forcesingle || al.isPairedEnd()==false){
    //single-end
  }
  else{
    //paired-end
    if(al.rnext=="=" || al.rnext==al.rname){
      if(al.pos<al.pnext){
        // the next read should appear in future
        //topair[al.qname]=range_t(al.pos,long(start.size()-1));
        topair[al.pnext][al.qname]=long(start.size()-1);
      }
      else{
        if(topair.count(al.pos)>0 ){
          if(topair[al.pos].count(al.qname)>0){
            int previd=topair[al.pos][al.qname];
            pairs[previd]=int(start.size()-1);
            pairs.back()=previd;
            //if(REPLACE_READ_NAME){allalign.back().qname=allalign[previd].qname;}
            topair[al.pos].erase(al.qname);
          }
          //old code
         // range_t nr=topair[al.qname];
         // if( nr.first==al.pnext){
         //   //pairs.push_back(topair[al.qname]);
         //   //pairs.push_back(int(start.size()-1));
         //   int previd=(int)nr.second;
         //   pairs[previd]=int(start.size()-1);
         //   pairs.back()=previd;
         //   topair.erase(al.qname);
         // }
        }
      }
    }
    //remove past topair records
    //map<long,map<string,long> >::iterator mit=topair.begin();
    //while(mit!=topair.end()){
    //  if(mit->first<al.pos) topair.erase(mit++);
    //  else break;
    //}
  }//end if
  return 0;
}


/**
 * Get the coverage from a set of reads, but skipping exons and paired-end read spans are considered covered
 * cvg must be cleaned before calling.
 * Notice: all coordinates in start and end is INCLUSIVE
 */
void ReadGroup::getConnectedCoverage(map<long,int>&cvg){
  cvg.clear();
  map<long,int>::iterator mitr;

  //get the boundaries
  for(int i=0;i<start.size();i++){
    if(pairs[i]==-1){
      cvg[start[i].front()]=0;
      cvg[end[i].back()+1]=0;
    }
    else{
      if(i<pairs[i])
        cvg[start[i].front()]=0;
      else
        cvg[end[i].back()+1]=0;
    }
  }      
  long startp,endp;
  for(int i=0;i<start.size();i++){
    if(pairs[i]==-1){
      startp=start[i].front();
      endp=end[i].back()+1;
    }else{
      if(i<pairs[i]){
        startp=start[i].front();
        endp=end[pairs[i]].back()+1;
      }else{
        startp=endp=-1;
      }
    }
    if(startp!=-1){
      mitr=cvg.find(startp);
      do{
        mitr->second=mitr->second+1;
        mitr++;
      }while(mitr->first!=endp);
    }
  }
  //condense
  int prevv=-1;
  for(mitr=cvg.begin();mitr!=cvg.end();){
    if(prevv!=mitr->second){
      prevv=mitr->second;
      mitr++;
    }
    else{
      cvg.erase(mitr++);
    }
  }
}


/**
 * Get the coverage from a set of reads
 * cvg must be cleaned before calling.
 * Notice: all coordinates in start and end is INCLUSIVE
 */
void ReadGroup::getCoverage(map<long,int>&cvg){
  getCoverage(cvg,range_t(-1,-1));
}

/**
 * Get the coverage from a set of reads for a given range
 * cvg must be cleaned before calling.
 * Notice: all coordinates in start and end is INCLUSIVE
 * The default range is (-1,-1).
 */
void ReadGroup::getCoverage(map<long,int>&cvg, const range_t & crange){
  cvg.clear();
  map<long,int>::iterator mitr;
  long a=crange.first,b=crange.second;
  // if(b<a){
    // cerr<<"Error: incorrect range in getCoverage():["<<a<<","<<b<<"].\n";
  // }
  //get the boundaries
  for(int i=0;i<start.size();i++){
    if(validRead[i]==0 )continue;
    if(b!=-1 && start[i][0]>b )continue;
    if(a!=-1 && end[i].back()<a )continue;
    for(int j=0;j<start[i].size();j++){
      long startp=start[i][j],endp=end[i][j]+1;
      cvg[startp]=0;
      cvg[endp]=0;
    }
  }
  if(a!=-1)cvg[a]=0;
  if(b!=-1)cvg[b]=0;
  //set boundaries
  for(int i=0;i<start.size();i++){
    if(validRead[i]==0)continue;
    if(b!=-1 && start[i][0]>b )continue;
    if(a!=-1 && end[i].back()<a )continue;
    for(int j=0;j<start[i].size();j++){
      long startp=start[i][j],endp=end[i][j]+1;
      mitr=cvg.find(startp);
      while(true){
        mitr->second=mitr->second+1;
        mitr++;
        if(mitr->first==endp)break;
      }
    }
  }
  //condense
  int prevv=-1;
  for(mitr=cvg.begin();mitr!=cvg.end();){
    if(prevv!=mitr->second || mitr->first==a || mitr->first==b){
      prevv=mitr->second;
      mitr++;
    }
    else{
      cvg.erase(mitr++);
    }
  }
  //finally, if the range is valid, keep only values between the range
  if(a!=-1 && b!=-1){
    for(mitr=cvg.begin();mitr!=cvg.end();){
      if(mitr->first>=a && mitr->first<=b){
        mitr++;
      }
      else{
        cvg.erase(mitr++);
      }
    }
  }
}

void ReadGroup::getPointCoverage(map<long,int>& cvg){
  cvg.clear();
  for(int i=0;i<start.size();i++)
    cvg[start[i][0]]++;
}

string ReadGroup::getChr() const{
  return statChr;
}


  /*
  Get the current range
  */
range_t ReadGroup::getRange() const{
  if(isfixrange) return fixrange;
  range_t r(LONG_MAX,-1);
  for(int i=0;i<start.size();i++){
    if(r.first>start[i].front()) r.first=start[i].front();
    if(r.second<end[i].back()) r.second=end[i].back();
  }
  return r;
}

/*
Split the read group by calculating range sets
*/
void ReadGroup::splitByRangeSet(vector<ReadGroup> & vr,long minD){
  vr.clear();
  //iterate the range sets
  vector<RangeSet> vs;
  
  for(int i=0;i<start.size();i++){
    if(validRead[i]==0)continue;
    bool newrange=true;
    //check if overlap with existing range sets
    for(int j=0;j<vs.size();j++){
      if(pairs[i]==-1){
        long cmdist=vs[j].minDist(start[i],end[i]);
        if(cmdist<=minD){
          vs[j].add(start[i],end[i]);
          vr[j].addOnly(start[i],end[i],direction[i]);
          newrange=false;
          break;
        }
      }
      else{
        if(i<pairs[i]){
          if(vs[j].minDist(start[i],end[i])<=minD
            || vs[j].minDist(start[pairs[i]], end[pairs[i]])<=minD
          ){
            vs[j].add(start[i],end[i]);
            vs[j].add(start[pairs[i]],end[pairs[i]]);
            vr[j].addPair(start[i],end[i],direction[i],
              start[pairs[i]],end[pairs[i]],direction[pairs[i]]);
            //vr[j].add(allalign[i]);
            //vr[j].add(allalign[pairs[i]]);
            newrange=false;
            break;
         }
        }
      }
    }//end for

    if(newrange){
      if(pairs[i]==-1 || i< pairs[i]){
        vr.push_back(ReadGroup());
        vr.back().forcesingle=forcesingle;
        vs.push_back(RangeSet());
        if(pairs[i]!=-1){
          vs.back().add(start[i],end[i]);
          vs.back().add(start[pairs[i]],end[pairs[i]]);
          //vr.back().add(allalign[pairs[i]]);
          vr.back().addPair(start[i],end[i],direction[i],
            start[pairs[i]],end[pairs[i]],direction[pairs[i]]);
        }else{
          vs.back().add(start[i],end[i]);
          vr.back().addOnly(start[i],end[i],direction[i]);
        }
       
      }
    }//end if
  }//end for
  if(VERBOSE==1){ 
     // cout<<"Merging "<<vs.size()<<" ranges...\n";
  }   
  //merge rangesets
  vector<int> valid(vs.size(),1);
  bool hasmerged=0;
  while(true){
    hasmerged=0;
    for(int i1=0;i1<vs.size();i1++){
      if(valid[i1]==0)continue;
      for(int i2=i1+1;i2<vs.size();i2++){
        if(valid[i2]==0)continue;
        //check range first
        //if(checkoverlap(vr[i1].getRange(), vr[i2].getRange())==false) continue;
        //check if two range sets overlap
        if(vs[i1].minDist(vs[i2])<=minD){
          //merge
          vs[i1].add(vs[i2]);
          for(int i=0;i<vr[i2].start.size();i++){
            int p=vr[i2].pairs[i];
            if(p==-1)
              vr[i1].addOnly(vr[i2].start[i],vr[i2].end[i],vr[i2].direction[i]);
            else{
              if(i<p){
                vr[i1].addPair(vr[i2].start[i],vr[i2].end[i],vr[i2].direction[i],
                   vr[i2].start[p],vr[i2].end[p],vr[i2].direction[p]);
              }
            }//end if
          }//end for
          valid[i2]=0;
          hasmerged=1;
        }
      }//end i2
    }//end i1
    if(hasmerged==0)break;
  }
  //if the direction is 0, try to split into two subgenes
  int cnv=valid.size();
  for(int i=0;i<cnv;i++){
    if(valid[i]==0)continue;
    vector<ReadGroup> subvr;
    vector<range_t> vrt;
    vr[i].splitByDirection(subvr,vrt);
    if(subvr.size()>0){
      valid[i]=0;
      // cout<<"SUPERRANGESET:"<<vs[i]<<endl;
      for(int j=0;j<subvr.size();j++) {
        vr.push_back(subvr[j]);
        valid.push_back(1);
        //reconstruct the rangeset
        RangeSet rs=vs[i].subRangeSet(vrt[j]);
        // cout<<"SUBRANGESET:"<<rs<<endl;
        vs.push_back(rs);
      }
    }
  }
  // print
  for(int i=0;i<valid.size();i++){
    if(valid[i]!=0){
      if(VERBOSE==1){
        // cout<<"SubInstance "<<i<<", reads:"<<vr[i].size()<<", range:"<<vs[i]<<endl;
      }
      vr[i].setRangeSet(vs[i]);
      vr[i].clearPairInfo();
    }
  }
  vector<ReadGroup>::iterator rgi=vr.begin();
  for(int i=0;i<valid.size();i++){
    if(valid[i]==0) rgi=vr.erase(rgi);
    else rgi++;
  }
  
}

/*
Split the read group by providing range sets
*/
void ReadGroup::splitByRangeSet(vector<ReadGroup> & vr,vector<RangeSet> & vs){
  vr.clear();
  vr.resize(vs.size());
  for(int i=0;i<start.size();i++){
    if(validRead[i]==0){
      continue;
    }
    int naccepted=0;
    //check if overlap with existing range sets
    for(int j=0;j<vs.size();j++){
      if(pairs[i]==-1){
        if(vs[j].overlapLen(start[i],end[i])==getReadLen(i)){
          vr[j].addOnly(start[i],end[i],direction[i]);
          naccepted++;
          //break; //If the range is non-overlapping, uncomment this break
        }else{
        }
      }
      else{
        int p=pairs[i];
        if(i<p){
          if(vs[j].overlapLen(start[i],end[i])==getReadLen(i)
            || vs[j].overlapLen(start[p], end[p])==getReadLen(p)
          ){
            vr[j].addPair(start[i],end[i],direction[i],
              start[p],end[p],direction[p]);
            naccepted++;
            //break;   //If the range is non-overlapping, uncomment this break
         }else{
         }
        }
      }
    }
    if(naccepted==0){
      if(start[i].size()>1){
        // cout<<"Junction read skipped: "; for(int j=0;j<start[i].size();j++) cout<<"["<<start[i][j]<<","<<end[i][j]<<"] "; cout<<endl;
      }
    }
  }//end for
  // print
  for(int i=0;i<vr.size();i++){
      vr[i].setRangeSet(vs[i]);
      vr[i].setRange(vs[i].toSingleRange());
      // cout<<"SubInstance "<<i<<", reads:"<<vr[i].size()<<", range:"<<vs[i]<<endl;
  }
}

/* 
Split the read group into two sub groups by checking the direction of the junctions.
For a read group of two overlapping opposite genes, this function will split them
Return value:
  vr: the splitted read group. If not splitted, this is empty
  vrt: the corresponding range of the two read groups. the size is equal to vr.
*/
void ReadGroup::splitByDirection(vector<ReadGroup>& vr,vector<range_t>& vrt){
  //check the number of positive/negative junctions
  //we assume that the gene can be separated into 2 groups of pos/neg junctions
  int npos=0;int nneg=0;
  range_t posranges(-1,-1),negranges(-1,-1);
  int lastdir=0;
  int nbadjunc=0;
  //here, we try to cluster the positive and negative junction reads
  for(int i=0;i<start.size();i++){
    if(validRead[i]==0 || start[i].size()<=1 || direction[i]==0)continue;
    if(lastdir!=0 && direction[i]!=lastdir){
      //for those junctions that are opposite to the last junction, count them as 'bad'
      nbadjunc++;
      if(nbadjunc>6)return;
      if(posranges.first!=-1 && negranges.first!=-1) continue;
    }
    if(direction[i]>0){
      npos++; 
      if(posranges.first==-1 || start[i][0]<posranges.first) posranges.first=start[i][0];
      if(posranges.second==-1 || end[i].back()>posranges.second) posranges.second=end[i].back();
    }else{
      nneg++;
      if(negranges.first==-1 || start[i][0]<negranges.first) negranges.first=start[i][0];
      if(negranges.second==-1 || end[i].back()>negranges.second) negranges.second=end[i].back();
    }
    lastdir=direction[i];
  }
  if(npos<2 || nneg<2)return; 
  // cout<<"SPLITBYDIRECTION: npos="<<npos<<",nneg="<<nneg<<", posrange=["<<posranges.first<<","<<posranges.second<<"], negrange=["<<negranges.first<<","<<negranges.second<<"].\n";
  //check if both ranges are valid and don't overlap?
  if(posranges.first!=-1 && negranges.first!=-1){
    if(posranges.first<=negranges.second && negranges.first <=posranges.second) {
       //they do overlap; further actions required
       return;
    }
    //extend ranges?
   //analyze the coverage information
    range_t mergarea(posranges.second,negranges.first);
    if(posranges.second>negranges.first) mergarea=range_t(negranges.second,posranges.first);
    map<long,int> cvg;
    // cout<<"Checking range ["<<mergarea.first<<","<<mergarea.second<<"]...\n";
    getCoverage(cvg,mergarea);
    //If there are 0 ranges
    long zerol=-1,zeror=-1;
    bool reach0=false;
    for(map<long,int>::iterator cvgit=cvg.begin();cvgit!=cvg.end();cvgit++) {
      if(cvgit->second==0 && zerol==-1) zerol=cvgit->first;
      if(cvgit->second==0) reach0=true;
      if( cvgit->second>0 && reach0) {zeror=cvgit->first; reach0=false;}
    }
    // modify the range
    if(zerol!=-1 && zeror!=-1){
    }else{
      zerol=mergarea.first/2+mergarea.second/2;
      zeror=zerol;
    }

    if(posranges.first>zeror && zerol>negranges.second  ){posranges.first=zeror; negranges.second=zerol;}
    if(negranges.first>zeror && zerol>posranges.second  ){negranges.first=zeror; posranges.second=zerol;}

    //split by range set
    // cout<<"Split range by direction: ["<<posranges.first<<","<<posranges.second<<"] and ["<<negranges.first<<","<<negranges.second<<"]\n";
    vector<RangeSet>  vs(2);
    vs[0].add(posranges); vs[1].add(negranges);
    vrt.push_back(posranges);
    vrt.push_back(negranges);
    splitByRangeSet(vr,vs);
  }
}

/* Split the read group into different subgroups, using provided ranges. 
Notice that the ranges are stored as readgroup.
void ReadGroup::splitByRangeSet(vector<ReadGroup> & vr,ReadGroup& prange,double rangeonly){
  //setting up different range sets
  vector<RangeSet> annors(prange.size());
  vpos_t & ans=prange.s();
  vpos_t & ane=prange.e();
  for(int i=0;i<prange.size();i++){
    annors[i].add(ans[i],ane[i]);
  }
  //split by current range set
  splitByRangeSet(vr,annors);
  //If reads are further required to stay within the segment regions of the ranges?
  if(!rangeonly){
  }
}
*/


void ReadGroup::setRange(range_t  s){
  fixrange=s;
  isfixrange=true;
}
void ReadGroup::setRangeSet(RangeSet & s){
  rg=s;
}
  
void ReadGroup::calculateRangeSet(){
  rg.clear();
  for(int i=0;i<start.size();i++){
    rg.add(start[i],end[i]);
  }
}

void ReadGroup::setupBound(vector<range_t>& jbound){
  range_t currentrange=getRange();
  allbound[currentrange.first]=1;
  allbound[currentrange.second+1]=1;
  for(int i=0;i<jbound.size();i++){
    allbound[jbound[i].first]=1;
    allbound[jbound[i].second]=1;
  }
}
 
/**
 * Get all boundaries of given reads
 * minjunc: minimum junction count
 * minadjrange: for coverage search, the minimum distance between ranges
 * In all bound, the boundary will be given. There are 3 types of the boundary: junction boundary (marked by 0) and low coverage boundary (marked by 1), and gene annotation boundary (marked by 2)
 */
void ReadGroup::calculateBound( 
    bool cvgcut, //if new boundary should be generated from coverage cutoff
    int minjunc,
    int minadjrange
    ){
  //this variable is used to set if we need to insert boundary for zero coverage points
  allbound.clear();
  vpos_t & rpoolstart=start;
  vpos_t & rpoolend=end;
  map<long,int> boundfilter;
  range_t currentrange=getRange();
  boundfilter[currentrange.first]=minjunc;
  boundfilter[currentrange.second+1]=minjunc;
  for(int i=0;i<rpoolstart.size();i++){
    if(rpoolstart[i].size()>1){
      //Ranges are [start,end)
      for(int j=1;j<rpoolstart[i].size();j++) boundfilter[rpoolstart[i][j]]++;
      for(int j=0;j<rpoolend[i].size()-1;j++) boundfilter[rpoolend[i][j]+1]++;
    }
  }
  //
  map<long,int>::iterator mitr,mitr2;
  for(mitr=boundfilter.begin();mitr!=boundfilter.end();mitr++){
    if(mitr->second>=minjunc)
      allbound[mitr->first]=0;//set up boundary. ok to reset 2 to 0
  }
    //now, check the coverage (real), insert boundaries if the coverage drops below 0.
    //but do not insert boundaries if the maximum coverage is too low, or if the gap is too short
  if(cvgcut==false)return;
  map<long,int>cvg;
  //get the real coverage
  getCoverage(cvg);
 
  vector<range_t> cutpoint;
  cutcvg(cvg,cutpoint,1,minadjrange);
    
  // now, check adjacent ranges,
  // if they're too close, merge them
  //int mincutseglen=statReadLen;
  int mincutseglen=2; //change here compared to isolasso. Replace 0 by 2 for possible sub-exons created by coverage cut-off. 
  if(cutpoint.size()>0){
    for(int i=0;i<cutpoint.size();i++){
        if(allbound.count(cutpoint[i].first)==0 
          && getClosestPoint(cutpoint[i].first, allbound)> mincutseglen )
          allbound[cutpoint[i].first]=1; //begin
        if(allbound.count(cutpoint[i].second+1)==0
         && getClosestPoint(cutpoint[i].second+1, allbound)> mincutseglen )
          allbound[cutpoint[i].second+1]=1; //begin
    }
  }
  getCvgStat(cvg,allbound,cvgstat);
}
/*
Convert boundary to segments.
Please keep the segments continuous. Otherwise, there will be problems.
*/
void ReadGroup::bound2Seg(){
  map<long,int>::iterator mit,mit2;
  mit=allbound.begin(); mit2=mit; mit2++;
  while(mit2!=allbound.end()){
    segs.push_back(range_t(mit->first,mit2->first-1));
    mit++;
    mit2++;
  }
  validSeg=vector<int>(segs.size(),1);
}

type_t ReadGroup::getType(int n)const{
  const pos_t &sn=start[n];
  const pos_t &en=end[n];
  return getType(sn,en);
}
type_t ReadGroup::getType(const pos_t& sn, const pos_t& en)const{
  type_t ret;
  long totallen=0;
  for(int i=0;i<sn.size();i++){
    for(int j=0;j<segs.size();j++){
      long d=getoverlaplength(sn[i],en[i],segs[j].first,segs[j].second);
      totallen+=d;
      if(d>0){
        if( d<MIN_SEG_FL_OVERLAP){
         if(i==0 &&   en[i]>segs[j].second 
             && segs[j].second>sn[i] && sn[i]>segs[j].first ) 
           continue;  
         if( i==sn.size()-1  && sn[i]<segs[j].first 
            && segs[j].first<en[i] && en[i]<segs[j].second ) continue;
        }
        ret.push_back(j);
      }
    }
  }
  int readlen=getReadLen(0);//getreadlen(n)
  if(readlen!=totallen){
    //if it's not completely overlap, this is an invalid type
    //cerr<<"Warning: read "<<n<<" are not completely overlap.\n";
    //cerr<<"Read: size "<<sn.size()<<":"<<endl;
    //for(int i=0;i<sn.size();i++) cerr<<"["<<sn[i]<<","<<en[i]<<"],"; cerr<<endl;
    //cerr<<"Segs: size "<<segs.size()<<":"<<endl;
    //for(int i=0;i<segs.size();i++) cerr<<"["<<segs[i].first<<","<<segs[i].second<<"],"; cerr<<endl;
  }
  
  return ret;
}

/*
Calculate the types according to the boundary
*/
void ReadGroup::calculateType(){
  bound2Seg();
  //print segs
  //cout<<"SEGS:"<<endl; for(int i=0;i<segs.size();i++) cout<<"["<<segs[i].first<<","<<segs[i].second<<"] "; cout<<endl;
  alltype.clear();typecount.clear();read2type.clear();
  map<type_t,int> knowntype;
  //also count the number of reads for each segments
  segsCnt=vector<int>(segs.size(),0);
  segsJCnt=segsCnt;
  for(int i=0;i<start.size();i++){
    type_t cty=getType(i);
    if(cty.size()>1){
      //print type
      //cout<<"JUNCTION TYPE:"; for(int j=0;j<cty.size();j++) cout<<cty[j]<<" ";cout<<endl;
    }
    //count segments
    for(int j=0;j<cty.size();j++){
        segsCnt[cty[j]]++;
        if(start[i].size()>1) segsJCnt[cty[j]]++;
    }
    //count types
    int tn=-1;
    if(knowntype.count(cty)>0){
    }else{
      int c=(int)knowntype.size();
      knowntype[cty]=c;
      alltype.push_back(cty);
      typecount.push_back(0);
      typedir.push_back(0);
    }
    tn=knowntype[cty];
    typecount[tn]++;
    typedir[tn]+=direction[i];
    read2type.push_back(tn);
  }//end for
  //set up valid SGType
  validSGType=vector<int>(alltype.size(),1);
  for(int i=0;i<alltype.size();i++){
    if(alltype[i].size()==0)validSGType[i]=0;
  }
}

int ReadGroup::validSize() const{
  int n=0;
  for(int i=0;i<validRead.size();i++)
    if(validRead[i]==1) n++;
  return n;
}


/* Remove segments that are too weak 
   But don't remove those that are included in the junction
*/
void ReadGroup::removeWeakSegs(float minf){
  vi isinjunc(segs.size(),0);
  //check if segments are included in the junctions?
  for(int i=0;i<alltype.size();i++){
    if(alltype[i].size()>1){
      if(alltype[i][0]!=alltype[i][1]-1)
        isinjunc[alltype[i][0]]++;
      if(alltype[i].back()!=alltype[i][alltype[i].size()-2]+1)
        isinjunc[alltype[i].back()]++;
      if(alltype[i].size()>2){
        for(int j=1;j<alltype[i].size()-1;j++){
          int k=alltype[i][j];
          isinjunc[k]++;
        }
      }
    }
  }
  for(int j=0;j<segs.size();j++){
    float maxcvg=0;
    if(segsCnt[j]==0){
      // if(VERBOSE==1){ 
        // cout<<"Invalid seg "<<j<<" because of 0 read count.\n";
      // }
      validSeg[j]=0;
    }
    //if this is not covered by junctions (only available for given junctions
    //get the coverage of flanking exons
    if(j>0 && cvgstat.count(segs[j-1].first)>0)
      maxcvg=cvgstat[segs[j-1].first][4];
    if(j<segs.size()-1 && cvgstat.count(segs[j+1].first)>0){
      float m2=cvgstat[segs[j+1].first][4];
      if(m2>maxcvg) maxcvg=m2;
    }
    if(cvgstat.count(segs[j].first)>0 ){
      if(cvgstat[segs[j].first][4]<=maxcvg*minf){
        if(isinjunc[j]<1){ //add another condition that this segment is not used by any junctions
          // if(VERBOSE==1){
	    // cout<<"Invalid seg "<<j<<" because of too small coverages ("<<cvgstat[segs[j].first][4]<<"/"<<maxcvg<<").\n";
	  // }
	  validSeg[j]=0;
        }
      }
    }
  }

}


/*
Determine if  segments are valid or not.
*/
void ReadGroup::calculateValidSegs(){
  bool haschange=false;
  do{
    //Recalculate SGType 
    for(int i=0;i<alltype.size();i++){
      //if(alltype[i].size()==0)validSGType[i]=0;
      for(int j=0;j<alltype[i].size();j++){
        if(validSeg[alltype[i][j]]==0){
          validSGType[i]=0;
          break;
        }
      }//end for j
    }//end for i
    //recalculate valid reads based on SGType
    for(int j=0;j<read2type.size();j++){
      if(validSGType[read2type[j]]==0)
        validRead[j]=0;
    }
    for(int i=0;i<segsCnt.size();i++) segsCnt[i]=0;
    //recalculate segsCnt
    for(int i=0;i<start.size();i++){
      if(validRead[i]==0)continue;
      int ttype=read2type[i];
      for(int j=0;j<alltype[ttype].size();j++)
        segsCnt[alltype[ttype][j]]++;
    }
    haschange=false;
    //invalidate segs with 0 counts
    for(int j=0;j<segs.size();j++){
      if(validSeg[j]==0)continue;
      if(segsCnt[j]==0){
        validSeg[j]=0;
        haschange=true;
      }
    }
  }while(haschange==true);
  getCvgStatistics();
}

void ReadGroup::getCvgStatistics(){
  //get the coverage stat
  map<long,int>cvg;
  getCoverage(cvg);
  cvgstat.clear();
  getCvgStat(cvg,allbound,cvgstat);
}

/* Return the index of reads containing specified segments */
vector<int> ReadGroup::getReadsContainingSegs(int segid){
  map<int,int> typeshavingreads;
  for(int i=0;i<alltype.size();i++){
    for(int j=0;j<alltype[i].size();j++)
      if(alltype[i][j]==segid)
        typeshavingreads[i]=1;
  }
  vector<int> ret;
  for(int i=0;i<read2type.size();i++)
    if(typeshavingreads.count(read2type[i])>0)
      ret.push_back(i);
  return ret;
}
 
/*
  Return segments with at least one read covered.
  Call calculateValidSegs() before calling this.
*/
vector<range_t> ReadGroup::getSegs()const{
  vector<range_t> r;
  for(int i=0;i<segs.size();i++)
    if(validSeg[i]>0) r.push_back(segs[i]);
  return r;
}


int ReadGroup::getReadLen(int n) const{
  int l=0; if(n>=start.size())return l;
  for(int i=0;i<start[n].size();i++)
    l+=end[n][i]-start[n][i]+1;
  return l;

}

void ReadGroup::toStream(ostream &out)const {
  range_t currentrange=getRange();
  char cdir='.';
  if(getDir()==1) cdir='+';
  if(getDir()==-1)cdir='-';
  out<<"Boundary\t"<<getChr()<<"\t"<<currentrange.first<<"\t"<<currentrange.second<<"\t"<<cdir<<endl;
  if(size()>0)
    out<<"ReadLen\t"<<getReadLen()<<endl;
  else
    out<<"ReadLen\t"<<statReadLen<<endl;
  
  //segments and its statistics
  int nsg=0;
  for(int i=0;i<segs.size();i++) if (validSeg[i]==1)nsg++;
  out<<"Segs\t"<<nsg<<endl;
  for(int i=0;i<segs.size();i++){
    if(validSeg[i]==0)continue;
    out<<segs[i].first<<"\t"<<segs[i].second<<"\t"<<(segs[i].second-segs[i].first+1)<<"\t"<<segsCnt[i];
    //stat
    //vector<double>& stat=cvgstat[vs[i].first];
    //for(int j=0;j<stat.size();j++) out<<"\t"<<stat[j];
    //out<<cvgstat[0][0];
    map<long, vector<double> >::const_iterator st=cvgstat.find(segs[i].first);
    for(int j=0;j<st->second.size();j++)out<<"\t"<<st->second.at(j);
    out<<endl;
  }
  //References
  ReadGroup rg=C_ANNO.getReadGroup(getChr(),currentrange);
  out<<"Refs\t"<<rg.size()<<endl;
  //output references separately
  //due to compatibility, don't output direction now
  for(int i=0;i<rg.size();i++){
    type_t crt=getType(rg.start[i],rg.end[i]);
    // if(VERBOSE==1){
      // cout<<"CRT: "; for(int j=0;j<crt.size();j++) cout<<crt[j]<<" ";cout<<endl;
      // cout<<"Coordinates:"; 
      // for(int j=0;j<rg.start[i].size();j++) cout<<rg.start[i][j]<<"-"<<rg.end[i][j]<<","; cout<<endl;
    // }
    //output crt
    type_t crtbin(nsg);
    for(int j=0;j<crt.size();j++) crtbin[crt[j]]=1;
    for(int j=0;j<crtbin.size();j++){
      out<<crtbin[j];
      if(j!=crtbin.size()-1)out<<" ";
      else out<<"\t";
    }
    //direction and ID
    if(rg.direction[i]==1)
      out<<"+";
    else
      out<<"-";
    out<<"\t";
    if(i<rg.readid.size())out<<rg.readid[i];
    else out<<getChr()<<"_"<<currentrange.first<<"_"<<currentrange.second<<"_"<<i;
    out<<endl;
  }

  //Reads
  out<<"Reads\t"<<validSize()<<endl;
  //SG Types
  //set up a map for valid segs
  map<int,int> validmap;
  int nvalid=0;
  for(int i=0;i<validSeg.size();i++){
    if(validSeg[i]!=0){
      validmap[i]=nvalid;
      nvalid++;
    }else
      validmap[i]=-1;
  }
  //then, write to stream
  //supress output containing invalid segments
  map<int,int> validSGmap;
  int nvalsg=0;
  for(int i=0;i<alltype.size();i++) {
    if(validSGType[i]!=0){
      validSGmap[i]=nvalsg;
      nvalsg++;
    }
    else
      validSGmap[i]=-1;
  }
  out<<"SGTypes\t"<<nvalsg<<endl;
  for(int i=0;i<alltype.size();i++){
    if(validSGType[i]==0){
      continue;
    }
    
    vector<int> exonbin(segs.size(),0);
    for(int j=0;j<alltype[i].size();j++) {
      exonbin[alltype[i][j]]=1;
    }
    //output
    for(int j=0;j<exonbin.size();j++)
      if(validSeg[j]==1) 
        out<<exonbin[j]<<" ";
    int cdir=0;
    if(typedir[i]>0)cdir=1;
    if(typedir[i]<0)cdir=-1;
    out<<typecount[i]<<"\t"<<cdir<<endl;
  }
  //PETypes
  map<range_t,pos_t>allpetypes;
  int nPE=0;
  for(int i=0;i<pairs.size();i++){
    if(pairs[i]!=-1 && i<pairs[i]){
      long fod=read2type[i],sod=read2type[pairs[i]];
      if(validSGType[fod]==0 || validSGType[sod]==0)
        continue;
      nPE++;
      long dist=start[pairs[i]].front()-end[i].back();
      allpetypes[range_t(fod,sod)].push_back(dist);
    }
  }
  out<<"PETypes\t"<<nPE<<"\t"<<allpetypes.size()<<endl;
  for(map<range_t,pos_t>::iterator apitr=allpetypes.begin();
            apitr!=allpetypes.end();
            apitr++){
      //In instance file, index+1 for all SGTypes
      //also, the index of SGTypes are chaged due to some invalid SGTypes
      //cluster the PETypes here to avoid over 30k numbers in a line
      map<long,int> pedistmap;
      for(int i=0;i<(apitr->second).size();i++)
        pedistmap[(apitr->second)[i]]++;
      
      out<<validSGmap[(apitr->first).first]+1<<"\t"
        <<validSGmap[(apitr->first).second]+1<<"\t"<<pedistmap.size()<<endl;
      for(map<long,int>::iterator pedit=pedistmap.begin();pedit!=pedistmap.end();pedit++) out<<pedit->first<<":"<<pedit->second<<" ";
      //for(int i=0;i<(apitr->second).size();i++)
      //  out<<(apitr->second)[i]<<" ";
      out<<endl;
  }
  //Coverage
  //for each type, count the number of reads for each position
  // ELSA modif 28fev14: do not write Coverage in the instance file, comment below:
//   vector<map<long,int> >typecvg(alltype.size());
//   for(int i=0;i<start.size();i++){
//     if(validRead[i]==0)continue;
//     typecvg[read2type[i]][start[i][0]]++;
//   }
// 
//   out<<"Coverage\t"<<nvalsg<<"\t"<<validSize()<<endl;
//   int counter=0;
//   for(int i=0;i<alltype.size();i++){
//     if(validSGType[i]==0)continue;
//     //count the number of positions having a specific number of reads
//     map<int,int> posstat;
//     map<long,int>& ttc=typecvg[i];
//     for(map<long,int>::iterator ittc=ttc.begin();ittc!=ttc.end();ittc++)
//       posstat[ittc->second]++;
//     out<<counter++<<"\t"<<posstat.size()<<endl;
//     for(map<int,int>::iterator ipss=posstat.begin();ipss!=posstat.end();ipss++)
//       out<<ipss->first<<","<<ipss->second<<"\t";
//     out<<endl;
//   }
  
}


ostream& operator<<(ostream& out,const ReadGroup & rg){
  rg.toStream(out);
  return out;
}

int ReadGroup::peSize()const{
  int n=0;
  for(int i=0;i<pairs.size();i++){
    if(pairs[i]!=-1) n++;
  }
  return n/2;
}

/*
  Called by removeTooLongReads().
  Check the number of reads between [start,end].
  If the number is greater than threshold, return true;
  else return false.
*/
bool isExceedThreshold(map<long,int> &cvg, long x, long y, int threshold){
  map<long,int>::iterator tk=cvg.lower_bound(x);
  int cnt=0;
  while(tk!=cvg.end() && tk->first<y){
    cnt+=tk->second;
    if(cnt>threshold) return true;
    tk++;
  }
  return false;
}
 
/* remove reads spanning too many reads 
   If two segments of a read contains >=minspan bases and contains >maxread reads, this read is marked invalid
*/
void ReadGroup::removeTooLongReads(long minspan,int maxread){
  map<long,int> cvg;
  //get the point coverage
  getPointCoverage(cvg);
  for(int i=0;i<start.size();i++){
    if(validRead[i]==0)continue;
    for(int j=1;j<start[i].size();j++){
      long span=start[i][j]-end[i][j-1];
      if(span<minspan)continue;
      if(isExceedThreshold(cvg,end[i][j-1],start[i][j],maxread)==true){
        validRead[i]=0;
        break;
      }
    }//end for j
    if(pairs[i]!=-1 && i<pairs[i]){
      long span=start[pairs[i]][0]-end[i].back();
      if(span<minspan)continue;
      if(isExceedThreshold(cvg,end[i].back(),start[pairs[i]][0],maxread)==true){
        validRead[i]=0;validRead[pairs[i]]=0;
      }
    }
  }//end for i
  for(int i=0;i<pairs.size();i++)
    if( validRead[i]==1 && pairs[i]!=-1 && validRead[pairs[i]]==0)
      validRead[i]=0;
}

int ReadGroup::getDirSum(){
  int r=0;
  int npos=0,nneg=0;
  map<range_t,range_t> splicingreads; // key: splicing positions; value: (nreads, positive-negative reads)
  for(int i=0;i<start.size();i++){
    if(validRead[i]==0 || start[i].size()==1 || direction[i]==0)continue;
    for(int j=1;j<start[i].size();j++){
      long k0=end[i][j-1];
      long k1=start[i][j];
      range_t rtpos(k0,k1);
      splicingreads[rtpos].first++;
      splicingreads[rtpos].second+=direction[i];
    }
  }
  int dir=0;
  //cout<<"DIR data:\n";
  int ntotal=0;
  for(map<range_t,range_t>::iterator srit=splicingreads.begin();
        srit!=splicingreads.end();srit++){
    //cout<<(srit->first).first<<"-"<<(srit->first).second<<":"<<(srit->second).first<<"-"<<(srit->second).second<<endl;
    int cr=(srit->second).second;
    dir+=cr;
    ntotal+=(srit->second).first;
  }
  if(ntotal==0)return 0;
  if(dir>0)dir=1;
  if(dir<0)dir=-1;
  //the number of bad reads and the number of bad junctions
  int nbadreads=0;
  int nbadjunctions=0;
  int nuniquejunction=0;//the bad junction that do not overlap with other junctions
  for(map<range_t,range_t>::iterator srit=splicingreads.begin();
        srit!=splicingreads.end();srit++){
    int currentdir=(srit->second).second;
    if(currentdir>0)currentdir=1;
    if(currentdir<0)currentdir=-1;
    if( dir*currentdir<0){
      if((srit->second).first>1){
        nbadreads+=(srit->second).first;
        nbadjunctions++;
        //check if this junction overlaps with other junctions?
        bool isunique=true;
        range_t tr1=srit->first;
        for(map<range_t,range_t>::iterator krit=splicingreads.begin();
          krit!=splicingreads.end();krit++){
           range_t tr2=krit->first;
           if(tr2==tr1)continue;
           if(tr1.first<=tr2.second && tr2.first<=tr1.second){//they overlap
             isunique=false;
             break;
           }
        }
        if(isunique) nuniquejunction++;
      }
    }
  }
  if(nuniquejunction>0){
    //cerr<<"nunique junction>0.\n";
    return 0;
  }
  return dir;
}

