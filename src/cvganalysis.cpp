#include <iostream>
#include "cvganalysis.h"

inline long abs2(long a){
	return (a>0?a:-a);
}



/**
 * Get the coverage from a set of reads, but skipping exons and paired-end read spans are considered covered
 * cvg must be cleaned before calling.
 * Notice: all coordinates in rpoolstart and rpoolend is INCLUSIVE
 */
void getconnectedcoverage(map<long,int>&cvg, vector<vector<long> >&rpoolstart, vector<vector<long> >&rpoolend,int appearpereads){
  map<long,int>::iterator mitr;

  //get the boundaries
  if(appearpereads==0)
    for(int i=0;i<rpoolstart.size();i++){
      cvg[rpoolstart[i].front()]=0;
      cvg[rpoolend[i].back()+1]=0;
    }
  else
    for(int i=0;i<rpoolstart.size();i+=2){
      cvg[rpoolstart[i].front()]=0;
      cvg[rpoolend[i+1].back()+1]=0;
    }

  //set boundaries
  long startp,endp;
  if(appearpereads==0)
    for(int i=0;i<rpoolstart.size();i++){
      startp=rpoolstart[i].front();
      endp=rpoolend[i].back()+1;
      mitr=cvg.find(startp);
      while(true){
        mitr->second=mitr->second+1;
        mitr++;
        if(mitr->first==endp)break;
      }
    }
  else
    for(int i=0;i<rpoolstart.size();i+=2){
      startp=rpoolstart[i].front();
      endp=rpoolend[i+1].back()+1;
      mitr=cvg.find(startp);
      while(true){
        mitr->second=mitr->second+1;
        mitr++;
        if(mitr->first==endp)break;
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
 * Notice: all coordinates in rpoolstart and rpoolend is INCLUSIVE
 */
void getcoverage(map<long,int>&cvg, vector<vector<long> >&rpoolstart, vector<vector<long> >&rpoolend){
  map<long,int>::iterator mitr;

  //get the boundaries
  for(int i=0;i<rpoolstart.size();i++){
    for(int j=0;j<rpoolstart[i].size();j++){
      long startp=rpoolstart[i][j],endp=rpoolend[i][j]+1;
      cvg[startp]=0;cvg[endp]=0;
    }
  }
  //set boundaries
  for(int i=0;i<rpoolstart.size();i++){
    for(int j=0;j<rpoolstart[i].size();j++){
      long startp=rpoolstart[i][j],endp=rpoolend[i][j]+1;
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
    if(prevv!=mitr->second){
      prevv=mitr->second;
      mitr++;
    }
    else{
      cvg.erase(mitr++);
    }
  }
}

/*
Generate breakpoints according to the coverage. 
Break points are generated if the coverage falls below some threshold, 
but if two break points are too close, supress the break point.
*/
int cutcvg(const map<long,int>& cvg, vector<range_t>& cutpoint,
  int threshold, int minadjrange){

   map<long,int>::const_iterator mitr;
    vector<range_t >cvgranges;
    range_t rentrange;
    bool comeintononzerorange=false;
    for(mitr=cvg.begin();mitr!=cvg.end();mitr++){
      if(comeintononzerorange==false){ 
        //ranges with low coverage; ignore until a higher range is reached.
        if(mitr->second>=threshold){
          //encounter a non-zero one
          rentrange.first=mitr->first;
          rentrange.second=mitr->first+1;
          cvgranges.push_back(rentrange);
          comeintononzerorange=true;
        }
        else{
          //do nothing, continue
        }
      }
      else{
        int nind=cvgranges.size()-1;//point to the last position
        //non-zero range
        if(mitr->second<threshold){
          //drop to zero range
          comeintononzerorange=false;
        }
        else{
          //if(mitr->first>currentrange.second){//exceed the current range, stop searching
          //  cvgranges[nind].second=currentrange.second+1;
          //  break;
          // }
         // cvgranges[nind].second=mitr->first;
          cvgranges[nind].second=mitr->first;
        }
      }
    }//end mitr for

    // now, check adjacent ranges,
    // if they're too close, merge them
    if(cvgranges.size()>0){
      bool prevjump=false;
      for(int i=0;i<cvgranges.size();i++){
        if(cutpoint.size()>0 && cvgranges[i].first-cutpoint.back().second<=minadjrange)
          cutpoint.back().second=cvgranges[i].second;
        else
          cutpoint.push_back(cvgranges[i]);
      }
    }
    return 0;

}


/*
Get the statistics of the coverage
May add extra  identical points to cvg
*/
int getCvgStat(map<long,int>&cvg, map<long,int>& allbound,
    map<long, vector<double> >&cvgstat){
    map<long,int>::iterator mitr,mitr2;

    //then, get some statistics of the ranges
    for(mitr=allbound.begin();mitr!=allbound.end();mitr++){
      //if the coverage data does not include the boundary point, insert the boundary point to coverage, 
      //and set it to the nearest value
      if(cvg.count(mitr->first)==0){
        pair<map<long,int>::iterator, bool> inspot=cvg.insert(pair<long,int>(mitr->first,0));

        mitr2=inspot.first;
        if(mitr2->first!=(cvg.begin())->first){
          map<long,int>::iterator mitr3=mitr2;
          mitr3--;
          mitr2->second=mitr3->second;
        }
      }
    }
    //iterate these ranges, get statstics
    mitr=allbound.begin();mitr2=mitr;mitr2++;

    while(mitr2!=allbound.end()){
      //identify corresponding positions in cvg
      map<long,int>::iterator cvgtr1=cvg.find(mitr->first);
      map<long,int>::iterator cvgtr2=cvg.find(mitr2->first);
      // if(cvgtr1==cvg.end() || cvgtr2==cvg.end()){
        // cerr<<"Error: cannot find values in cvg map.\n";
      // }
      double startcvg=cvgtr1->second;
      double zerofrac=0;
      map<long,int>::iterator iit=cvgtr1,iit2=cvgtr1;
      iit2++;
      double maxcvg=-1;
      double endcvg=0;
      double meancvg=0;
      while(iit!=cvgtr2){
        if(iit->second>maxcvg)maxcvg=iit->second;
        if(iit->second==0)zerofrac+=(iit2->first-iit->first);
        meancvg+=(iit2->first-iit->first)*iit->second;
        endcvg=iit->second;
        iit++;iit2++;
        if(iit2==cvg.end())break;
      }
      //if(iit->second>maxcvg)maxcvg=iit->second;

      vector<double> val(5,0);
      val[0]=maxcvg;val[1]=startcvg;val[2]=endcvg;val[3]=zerofrac/(cvgtr2->first-cvgtr1->first);
      val[4]=meancvg/(cvgtr2->first-cvgtr1->first);
      cvgstat[mitr->first]=val;

      //increase pointer
      mitr++;mitr2++;
    }


  return 0;
}

/*
 For a given point, get the closest position
 If fail, return -1
*/
long getClosestPoint(long pos,map<long,int>& bound){
  if(bound.size()==0)return -1;
  map<long,int>::iterator mit=bound.lower_bound(pos);
  if(mit!=bound.end()){
    long cm=abs2(mit->first-pos);
    while(mit!=bound.begin()){
      mit--;
      if(abs2(mit->first-pos)<cm)cm=abs2(mit->first-pos);
    }
    return cm;
  }else{
    return abs2(pos-bound.rbegin()->first);
  }
  
}



