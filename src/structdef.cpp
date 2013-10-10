#include <sstream>
#include <iostream>
#include <set>
#include <map>
#include "common.h"
#include "structdef.h"


/*
This is the structure to store gene range  and annotation
*/

/* the clustering range of gene range. Since gene ranges allow overlapping, these ranges are clustered into a larger range and stored in C_GENE_RANGE */
GeneRange C_GENE_RANGE;
/* provided gene ranges. each range may contain multiple segments */
Annotation C_RANGE_ANNO;

/* Reference annotation */
Annotation C_ANNO;

/* Provided junction boundaries */
GeneRange C_JUNCTION;

ReadGroup Annotation::EmptyReadGroup;

/*
Modify function
*/
int GeneRange::push_back(string c, range_t r){
  if(chrexistmap.count(c)==0){
    chrexistmap[c]=(int)chrmap.size();
    chrmap.push_back(c);
    rangepos.push_back(map<long,int>());
  }  
  int chrid=chrexistmap[c];
  chr.push_back(chrid);
  range.push_back(r);
  if(rangepos[chrid].count(r.first)==0){
    rangepos[chrid][r.first]=(int)chr.size()-1;
  }
  return chr.size()-1;
}

/* Sort and remove duplicates*/
int GeneRange::sort(){
  map< pair<string, long>, range_t > alls;
  for(int i=0;i<chr.size();i++)
    alls[pair<string,long>(chrmap[chr[i]],range[i].first)]=range[i];//(pair<string,range_t>(chrmap[chr[i]],range[i]));

  map<pair<string,long>,range_t >::iterator it;
  int ind=0;
  for(int i=0;i<rangepos.size();i++)
    rangepos[i].clear();
  chr.resize(alls.size());
  range.resize(alls.size());
  for(it=alls.begin();it!=alls.end();it++,ind++){
    chr[ind]=chrexistmap[((*it).first).first];
    range_t crange=(*it).second;
    range[ind]=crange;
    int chrid=chr[ind];
    if(rangepos[chrid].count(crange.first)==0){
      rangepos[chrid][crange.first]=ind;
    }
  }
  return 0;
}

/* Check if the ranges are sorted according to chromosome name and starting position */
bool GeneRange::check(){
  if(size()==0)return true;
  reset();
  string prevc=getChr();
  range_t prevr=getRange();
  while(true){
    inc();
    if(getIndex()>=size())break;
    string cc=getChr();
    range_t cr=getRange();
    if(cc<prevc || (cc==prevc && cr.first<prevr.first)){
      // cerr<<"Error: gene_range.check() failed. Prev "<<prevc<<":"<<prevr.first<<"-"<<prevr.second<<", Current "<<cc<<":"<<cr.first<<"-"<<cr.second<<endl;
      return false;
    }
    if(prevc=="" || cc==""){
      // cerr<<"Error: empty chromosome name.\n";
      return false;
    }
    prevc=cc;
    prevr=cr;
  }
  return true;
}

/*
Move the index to the next chromosome. Return -1 if reaches the end of the queue
*/
int GeneRange::moveToNextChr(string nextChr){
        if(chrexistmap.count(nextChr)==0)return -1;
        int nextchrid=chrexistmap[nextChr];
  while( !isEnd() && chr[index]!=nextchrid){
    //move to the next available range
    inc();
  }
  if(isEnd())return -1;
  else return index;
}

/*
Retrieve the next range, and save to range_t.
Return value:

  If normal case (rname=current chromosome name):
    set cr, and return 1.

  If reaches the end of the queue:
    set cr.first=cr.second=-1 and return -1;

  If reaches the end of the chromosome:
    do not increase, set range to 0, return 0.

  If rname implies a new chromosome name that is beyond the current chromosome name, 
    or the current range end position is smaller than the current position (pos):
      Set the cname/cr as the new  chromosome name/range, and return 2.
*/
int GeneRange::getNext(string rname,long pos, string& cname,range_t & cr){
  //reads are in advance; increase the gene range until it reaches the end
  if(isEnd() ){
    cname="";
    cr.first=cr.second=-1;
    return -1;
  }
  if(getChr()<rname ){
    cname=getChr();cr=getRange();
    inc();
    return 2;
  }
  if(getChr()==rname){
    cname=rname;
    //cr=getRange();
    cr=getRange();
    inc();
    if(cr.second<pos) return 2;
    else return 1;
  }
  else{
    //getChr()>rname: set to 0
    cname=getChr();
    cr.first=cr.second=0;
    return 0;
  }
}

/* Retrieve all ranges which are inside a given range t*/
int GeneRange::getAllRanges(string rname,range_t t, vector<range_t> &allr){
  if(chrexistmap.count(rname)==0)return -1;
  int mapcode=chrexistmap[rname];
  map<long,int> & crp=rangepos[mapcode];
  
  map<long,int>::iterator tk=crp.lower_bound(t.first);
  
  if(tk==crp.end())return 0;
  int ri1=tk->second;
  int ri2=ri1;
  while(tk!=crp.end() && tk->first<t.second){
    ri2=tk->second;
    tk++;
  }
  //cout<<"range.size()="<<range.size()<<",ri1-ri2:"<<ri1<<","<<ri2<<endl;
  for(int i=ri1;i<=ri2;i++){
    if(range[i].first>=t.first && range[i].second<t.second)
      allr.push_back(range[i]);
  }
  return 0;
}






 /*
 (OUTDATED)
  Add one record (rs,re) to the pool. 
  This function will automatically cluster the records, so the records must be added in ascending order by chr and 1st element of rs.
  Return 0 for success, -1 for failure (if the records are not added in order).
  int Annotation::add(string chr, pos_t rs, pos_t re){
     //If it's necessary to put the current group into record
     if(chr!=prevs || rs.front()>prevr.second){
        //write a new record
        map<range_t,pos_t > & crange=cluster[prevs];
        if(prevr.first!=-1){
          crange[prevr]=rangegroup;
        }
        rangegroup.clear();
        prevr.first=rs.front();
        prevr.second=re.back()-1;
    }
    else{
      if(rs.front()<prevr.first){
        return -1;
      }
    }
    if(chr!=""){
      ReadGroup & rg=allan[chr];
      rg.s().push_back(rs);
      rg.e().push_back(re);
      rangegroup.push_back(rg.s().size()-1);//save current 
      if(re.back()-1>prevr.second)prevr.second=re.back()-1;
    }
    prevs=chr;
    return 0;
  
  }
  */

/*
Put current clusters into the data. Only explicitly called at the end of the annotation input.
*/
int Annotation::cluster(){
  if(prevg.s().size()>0){ // jump the first one
    data[prevs][prevr]=prevg;
  }else{
  }
  //reset ReadGroup
  prevg.clear();
  prevr.first=prevr.second=-1;
  return 0;
}

/*
Add one record (rs,re) to the pool. 
This function will automatically cluster the records, so the records must be added in ascending order by chr and 1st element of rs.
Return 0 for success, -1 for failure (if the records are not added in order).
To finish the input, use cluster().
*/
int Annotation::add(string chr, pos_t rs, pos_t re,int dir,string id){
   //If it's necessary to put the current group into record
   if( chr>prevs || rs.front()>prevr.second){
      //write a new record
      cluster();
      //reset range
      prevr.first=rs.front();
      prevr.second=re.back()-1;
  }
  else{
    if(rs.front()<prevr.first || chr<prevs){
      return -1;
    }
  }
  //update current pool
  prevg.s().push_back(rs);
  prevg.e().push_back(re);
  prevg.getDirection().push_back(dir);
  prevg.getID().push_back(id);
  //prevg.add(rs,re,dir);
  if(re.back()-1>prevr.second)prevr.second=re.back()-1;
  prevs=chr;
  return 0;

}

/*
get the number of chromosomes
*/ 
int Annotation::getNChrom(){ return data.size(); }

vector<string> Annotation::getChrom(){
  vector<string> v;
  for(map<string, map<range_t, ReadGroup> >::iterator it=data.begin();
    it!=data.end();it++)
    v.push_back(it->first);
  return v;
}

  
/*
get the number of clusters
*/
int Annotation::getNCluster(){ 
  int n=0; 
  for(map<string, map<range_t, ReadGroup> >::iterator it=data.begin();it!=data.end();it++)
    n+=(it->second).size();
  return n;
}
/*
get a vector of clusters, specified by chromosome name
*/
vector<range_t> Annotation::getCluster(string chr){
  vector<range_t> v;
  if(data.count(chr)==0) return v;
  map<range_t,ReadGroup>& cr=data[chr];
  for(map<range_t,ReadGroup>::iterator it=cr.begin();
    it!=cr.end();it++)
    v.push_back(it->first);
  return v;
} 

/*
Get the corresponding read group from the chromosome name and the range
*/
ReadGroup Annotation::getReadGroup(string chr, range_t r){
  if(data.count(chr)==0) return EmptyReadGroup;
  if(data[chr].count(r)==0) return EmptyReadGroup;
  return data[chr][r];
}



/*
get the number of annotations
*/
int Annotation::getNAnnotation(){
  int n=0; 
  for(map<string, map<range_t, ReadGroup> >::iterator it=data.begin();it!=data.end();it++){
    for(map<range_t,ReadGroup>::iterator it2=(it->second).begin();
        it2!=(it->second).end();it2++)
    n+=(it2->second).s().size();
  }
  return n;
}

/*
For a given range (current), check if there are overlapping annotations. If no, do nothing; otherwise, extend the range to include annotations.
Return value:
  0 no range overlaps
  1 updated range
*/
int Annotation::checkOverlapRange(range_t& current, string rname){
  if(data.count(rname)==0) return -1;
  map<range_t, ReadGroup> & crg=data[rname];
  for(map<range_t,ReadGroup>::iterator it=crg.begin();
    it!=crg.end();it++){
    range_t  crt=it->first;
    if (crt.first>current.second) return 0;
    if (crt.first>=current.first 
      && crt.first<=current.second){
      //overlapping
      if(crt.second>current.second)
        current.second=crt.second;
      return 1;
    }
  }
  return 0;
}

