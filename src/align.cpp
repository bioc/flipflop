#include <iostream>
#include <sstream>
#include <cstdio>

#include "align.h"


ostream& operator<<(ostream& out, Align& a){
  range_t r=a.getRange();
  out<<a.qname<<"\t"<<a.rname<<":"<<r.first<<"-"<<r.second<<"\t";
  for(int i=0;i<a.start.size();i++) out<<a.start[i]<<"-"<<a.end[i]<<",";
  
  return out;
}


int Align::getReadLen()const{
  int n=0;
  for(int i=0;i<start.size();i++)
    n+=end[i]-start[i]+1;
  return n;
}

 /*
  parse fields from a line of string
  Return 0 if success, -1 if any error occurs.
  */
int Align::parse(string oneline){
    stringstream ss(oneline);
    ss>>qname>>flag>>rname>>pos>>mapq>>cigar>>rnext>>pnext>>plen;
    start.clear();
    end.clear();
    parsecigar();
    string twostr;
    ss>>twostr;
    ss>>twostr;
    int pos=ss.tellg();
    string rest=oneline.substr(pos+1);
    //parse XS:A:
    if(start.size()==1)splicedir=0;
    else{
      int xspos=rest.find("XS:A:");
      if(xspos!=string::npos){
        char x='.';
        string tx=rest.substr(xspos);
        sscanf(tx.c_str(),"XS:A:%c",&x);
        if(x=='+')splicedir=1;
        if(x=='-')splicedir=-1;
      }
    }
    //parse NM:i:
    int nmpos=rest.find("NM:i:");
    if(nmpos!=string::npos){
      string tx=rest.substr(nmpos);
      sscanf(tx.c_str(),"NM:i:%d",&nmismatch);
    }
    //parse NM:i:
    int nhpos=rest.find("NH:i:");
    if(nhpos!=string::npos){
      string tx=rest.substr(nhpos);
      sscanf(tx.c_str(),"NH:i:%d",&nhits);
    }
   
    if (ss.fail()) return -1;
    return 0;
}

/*
Return the current range of the reads
*/
range_t Align::getRange(bool forcesingle){
  if(start.size()==0)
    return range_t(0,0);
  if(isPairedEnd() && !forcesingle){
    if(pos<pnext)
      return range_t(start.front(),start.front()+plen-1);
    else
      return range_t(pnext, end.back());
  }
  else
    return range_t(start.front(),end.back());
}


/**
 * Parese CIGAR strings
 * notice that the range should be interpreted as [startpos, endpos], not [startpos,endpos).
 * Handle I and D Cigar characters
 * Handle S and X and = Cigar characters (soft-clipping, match, mismatch) # ELSA TO DO - test more - 
 */
void Align::parsecigar( ){
  pos_t & startpos=start; 
  pos_t &endpos=end;
  stringstream ss2(cigar);
  long tpos=pos;
  
  bool islastn=false;
  while(!ss2.eof()){
    int nof;
    char t;
    ss2>>nof;
    if(ss2.eof())break;
    ss2>>t;
    if(t=='M'){
      if(islastn==false){
        startpos.push_back(tpos);
        endpos.push_back(tpos+nof-1);
      }
      else{
        endpos.back()+=nof;
      }
      islastn=true; // the last segment is 'M'
      tpos=tpos+nof;
    }
    else if(t=='N'){
      tpos=tpos+nof;
      islastn=false; //the last segment is 'N'
    }
    else if(t=='I'){
    }
    else if(t=='D'){
      if(endpos.size()==0){
        // cerr<<"Error: cannot process CIGAR string "<<cigar<<endl;
      }
      else{
        endpos.back()+=nof;
        tpos=tpos+nof;
      }
    }
    else if(t=='S' || t=='X' || t=='='){
       // cerr<<"Encounter the string "<<t<<" in CIGAR "<<cigar<<endl;
       startpos.push_back(tpos);
       endpos.push_back(tpos+nof-1);
       tpos=tpos+nof;
       // cerr<<"String S at the beginning:"<<tpos<<"et:"<<islastn<<endl;
    }
    else{
      // cerr<<"Error: not supported CIGAR character "<<t<<" of CIGAR "<<cigar<<endl;
    }

  }//end of ss2 loop

}

