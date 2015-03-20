#include <fstream>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include "common.h"
#include "auxiliaryio.h"
#include "cvganalysis.h"
#include "instanceio.h"
#include "structdef.h"
#include "bedio.h"
#include "FileSplitter.h"

/**
FILE IOS
**/

ofstream wigfs;  //wigfs, the connected  coverage
ofstream realwigfs;//realwigfs, the real coverage
ofstream bedfs; //bedfs, all junctions
ofstream juncfs;//juncfs, summary of junctions
ofstream boundfs;//boundfs, the boundary of an instance

ofstream iofs; //out instance file

ofstream ofs; //read_info file
ofstream bofs; //bound file string
ofstream gofs; //gene range file

FileSplitter *fileSplitter;


/* Instance ID */
//int n_INST=0;

/*
Prepare the aux output file, including parsing the input arguments and open corresponding files
*/
void prepareAuxFile(vector<string> args,
  string prefix,
  bool& outinstance
  ){

  
  string wigfile="/dev/null";
  string bedfile="/dev/null";
  string realwigfile="/dev/null";
  string juncfile="/dev/null";
  string boundfile="/dev/null";


  //The prefix of the generated files
  for(int i=0;i<args.size();i++){
    if(i<args.size()-1 && (string(args[i])=="-o" || string(args[i])=="--prefix") ){
      prefix=string(args[i+1]);
    }
  }
  if(prefix=="-") prefix="isolasso";
  string outReadInfoFile="/dev/null"; 
  string outBoundFile="/dev/null"; 
  string outGeneRangeFile="/dev/null"; //output files
  string outInstanceFile=prefix+".instance"; //output instance file (optional)
  
  for(int i=0;i<args.size();i++){
    if(string(args[i])=="-a" ||string(args[i])=="--annotation" ){
      wigfile=prefix+".wig";
      bedfile=prefix+".bed";
      realwigfile=prefix+".real.wig";
      juncfile=prefix+".junc.bed";
      boundfile=prefix+".bound.bed";
      // cout<<"Wig file:"<<wigfile<<endl;
      // cout<<"Bed file:"<<bedfile<<endl;
      // cout<<"Real Wig file:"<<realwigfile<<endl;
      // cout<<"Junction file:"<<juncfile<<endl;
      // cout<<"Boundary file:"<<boundfile<<endl;

      wigfs.open(wigfile.c_str());
      if(!wigfs.is_open()){
	// cerr<<"Error opening wig file "<<wigfile<<endl;
	//return -1;
      }
      //header line of annotations
      wigfs<<"track type=bedGraph name=Read_Cvg description="<<prefix<<"\n";
      realwigfs.open(realwigfile.c_str());
      if(!realwigfs.is_open()){
	// cerr<<"Error opening wig file "<<realwigfile<<endl;
	//return -1;
      }
      //header line
      realwigfs<<"track type=bedGraph name=Read_Coverage description="<<prefix<<"\n";
      bedfs.open(bedfile.c_str());
      if(!bedfs.is_open()){
	// cerr<<"Error opening bed file "<<bedfile<<endl;
	//return -1;
      }
      //header line
      bedfs<<"track name=Junctions description="<<prefix<<"\n";


      juncfs.open(juncfile.c_str());
      if(!juncfs.is_open()){
	// cerr<<"Error opening bed file "<<juncfile<<endl;
	//return -1;
      }
      juncfs<<"track name=Junction_Summary description="<<prefix<<"\n";
  
      boundfs.open(boundfile.c_str());
      if(!boundfs.is_open()){
	// cerr<<"Error opening bed file "<<boundfile<<endl;
	//return -1;
      }
      boundfs<<"track name=Instance_Boundary description="<<prefix<<"\n";
    }

    if(string(args[i])=="-n" || string(args[i])=="--isoinfer"){
      outReadInfoFile=prefix+".readinfo"; 
      outBoundFile=prefix+".bound"; 
      outGeneRangeFile=prefix+".generange"; //output files

      // ISOINFER related files
      ofs.open(outReadInfoFile.c_str());
      if(!ofs.is_open()){
	// cerr<<"Error opening output file "<<outReadInfoFile<<endl;
	//return -1;
      }
      bofs.open(outBoundFile.c_str());
      if(!bofs.is_open()){
	// cerr<<"Error opening boundary output file "<<outBoundFile<<endl;
	//return -1;
      } 
      gofs.open(outGeneRangeFile.c_str());
      if(!gofs.is_open()){
	// cerr<<"Error opening gene range file "<<outGeneRangeFile<<endl;
	//return -1;
      }
    }
    if(string(args[i])=="-i"){//this option removed
      outInstanceFile=prefix+".instance";
      // cout<<"Instance file: "<<outInstanceFile<<endl;
    }
  }

  //output instance
  if(outInstanceFile!=""){
    outinstance=true;
  }
  if(outinstance==true){
    iofs.open(outInstanceFile.c_str());
    if(!iofs.is_open()){
      // cerr<<"Error opening instance file "<<outInstanceFile<<endl;
      outinstance=false;
    }

    fileSplitter=new FileSplitter(outInstanceFile, iofs); // 2015-02-25
    fileSplitter->startWritingHeader();

    //writing args to out instance file
    iofs<<"@CMD:";
    for (int i=0;i<args.size();i++)
      iofs<<args[i]<<" ";
    iofs<<endl;
    iofs<<"@Fields in Segs section:\n";
    iofs<<"@0: Segment start, 1: Segment end, 2: Segment length, 3:Reads falling onto this segment, 4: max coverage, 5: coverage on left, 6: coverage on right, 7: percentage of bases with 0 coverage, 8: mean coverage\n";

    fileSplitter->endWritingHeader(); // 2015-02-25

  }
}

void closeAuxFile(){

  iofs.close();

  gofs.close();

  bofs.close();

  ofs.close();

  wigfs.close();
  realwigfs.close();
  bedfs.close();
  juncfs.close();
  boundfs.close();


}


/*
ISOINFER OUTPUT files
*/

void write2grange(int rangecounter, range_t& nowrange, string rname,int duprcount){
      //write to generange file
      gofs<<rangecounter++<<"\t"<<rname<<"\t";
      gofs<<"+"<<"\t"; // don't know direction, use + instead
      gofs<<nowrange.first-1<<"\t"<<nowrange.second;//1st is 0-base, 2nd is 1-base
      gofs<<"\t"<<duprcount;
      gofs<<endl;
}


/**
Write one line information to readinfo 
*/
void write2os(ostream& ofs,
  string rname, int flag, 
  vector<long>& startpos, 
  vector<long>& endpos,
  bool writeendl=true){
    ofs<<rname<<"\t";
    if(flag & 0x0010)//forward
      ofs<<"+\t";
    else
      ofs<<"-\t";
    for(int i=0;i<startpos.size();i++){
      ofs<<startpos[i]-1; //.SAM is 1 based, but .read_info start coordinates are 0 based.
      if(i!=startpos.size()-1)
        ofs<<",";
      else
        ofs<<"\t";
    }
    for(int i=0;i<endpos.size();i++){
      ofs<<endpos[i];
      if(i!=endpos.size()-1)
        ofs<<",";
      else
        ofs<<"\t";
    }
    if(writeendl)
      ofs<<endl;

}

/*
write bound information to ostream boundfs
Notice that sam file is 1-based, but bound file should be 0-based.
*/
void writeBound(ostream& ofs, string rname, vector<long>&startpos, vector<long>& endpos){
    for(int i=1;i<startpos.size();i++){
      //Should be 0-based, so minus 1
      ofs<<rname<<"\t"<<"+"<<"\t"<<startpos[i]-1<<"\t"<<0<<endl;
    }
    for(int i=0;i<endpos.size()-1;i++){
      //Here, 1-based exon endpos is exactly 0-based intron start.
      ofs<<rname<<"\t"<<"+"<<"\t"<<endpos[i]<<"\t"<<1<<endl;
    }

}


/* To be replaced by a new version
Write one line in readinfo file.
Used for readSamFile().
*/

int writeoneline2readinfo(Align &al){
  write2os(ofs,al.rname,al.flag,al.s(),al.e(),true);
  writeBound(bofs,al.rname,al.s(),al.e());
  return 0;
}

/**
ISOLASSO aux files
*/


void writeCoverage( ReadGroup & rg){
    //connected coverage
  map<long,int>cvg;
  map<long,int>::iterator mitr,mitr2;
  string chrname=rg.getChr();
  if(rg.size()==0)return;
  //getcoverage(cvg,rpoolstart,rpoolend);
  if(wigfs.is_open()){
    rg.getConnectedCoverage(cvg);
    //write to wigfs
    mitr=cvg.begin();mitr2=mitr;mitr2++;
    while(mitr2!=cvg.end()){
      wigfs<<chrname<<"\t"<<mitr->first-1<<"\t"<<mitr2->first-1<<"\t"<<mitr->second<<endl; //the wig file should be 1-base, but it seems that in UCSC genome browser, the wig file is 0-base?
      mitr++;
      mitr2++;
    }
  }
  //real coverage
  cvg.clear();
  if(realwigfs.is_open()){
    rg.getCoverage(cvg);
    //write to wigfs
    mitr=cvg.begin();mitr2=mitr;mitr2++;
    while(mitr2!=cvg.end()){
      realwigfs<<chrname<<"\t"<<mitr->first-1<<"\t"<<mitr2->first-1<<"\t"<<mitr->second<<endl;
      mitr++;
      mitr2++;
    }
  }
}

void writeannotation(
    ReadGroup & rg,
    int instanceid){
  
  string chrname=rg.getChr();
  //write to file
  int nj=0;
  vpos_t &rpoolstart=rg.s();
  vpos_t &rpoolend=rg.e();
  vi& dir=rg.getDirection();

  bedrec rec;
  rec.chr=chrname;
  rec.score=1;
  rec.dir='+';
  rec.color="255,0,0";

  //junctions
  if(bedfs.is_open()){
    for(int i=0;i<rpoolstart.size();i++){
      if(rpoolstart[i].size()>1){//junction reads
        rec.segstart=rpoolstart[i];
        rec.segend=rpoolend[i];
        rec.start=rpoolstart[i][0];
        rec.end=rpoolend[i].back();
        rec.nsegs=rec.segstart.size();
        stringstream ss; ss<<(nj++);string snj; ss>>snj;
        stringstream ss2; ss2<<instanceid; string sid; ss2>>sid;
        rec.name=string("Inst")+sid+string("_Junc")+snj;;
        if(dir[i]>=0)rec.dir='+'; else rec.dir='-';
        writeBed(bedfs,rec);
        
      }
    }//end for
  }//end if

  map<range_t,int>sumjunc;
  map<range_t,int> dirjunc;
  map<range_t,int>::iterator jmit;
  if(juncfs.is_open()){
    for(int i=0;i<rpoolstart.size();i++){
      if(rpoolstart[i].size()>1){//junction reads
        //get stat
        for(int j=0;j<rpoolstart[i].size()-1;j++){
          long starti=rpoolstart[i][j+1];
          long endi=rpoolend[i][j];
          int idir=dir[i];
          sumjunc[range_t(endi,starti)]++;
          dirjunc[range_t(endi,starti)]+=(idir>=0?1:-1);
        }
      }
    }
    //write 
    int width=20;
    nj=0;
    for(jmit=sumjunc.begin();jmit!=sumjunc.end();jmit++){
        
      long b=(jmit->first).first; //end of 1st exon, 1-base inclusive
      long c=(jmit->first).second;//beginning of 2nd exon, 1-base inclusive
      long a=b-width+1;// beginning of 1st exon, 1-base inclusive
      long d=c+width-1;// end of 2nd exon, 1-base inclusive
      int djdr=dirjunc[jmit->first]; 
     
      juncfs<<chrname<<"\t"<<a-1<<"\t"<<d<<"\t";//chrom, chromstart, chromend
      juncfs<<"Inst"<<instanceid<<"_J"<<nj++<<"_C"<<jmit->second<<"\t";//4th, name
      juncfs<<jmit->second<<"\t";//5th, score
      if(djdr>0)
        juncfs<<"+\t";//6th,orientation
      else
        juncfs<<"-\t";//6th,orientation
      juncfs<<a-1<<"\t"<<d<<"\t";//7,8
      if(djdr>0)
        juncfs<<"255,0,255\t";//9th, color
      else
        juncfs<<"0,255,255\t";//9th, color
      juncfs<<2<<"\t";//10th, blockcounts
      //11th, sizes
      juncfs<<width<<","<<width<<"\t";
      //12th, starts
      juncfs<<0<<","<<c-a;
      juncfs<<endl;

    }
  }
}


//write to boundfs
//range is inclusive
void writeboundfs(range_t range,int instanceid,string chrname,
    //map<long,int>&allbound, vector<int>& validsegs
    vector<range_t>& allrange,
    int score, int dir
  ){
  bool plotexons=true;
  char bd='.';
  if(dir>0)bd='+';
  if(dir<0)bd='-';
  if(boundfs.is_open()){
    //map<long,int>::iterator mitr=allbound.begin();
    //map<long,int>::reverse_iterator rmitr=allbound.rbegin();
    long startr=allrange.front().first;
    long endr=allrange.back().second;


    int nbks=allrange.size();
    //for(int i=0;i<validsegs.size();i++)if(validsegs[i]==1)nbks++;
    //map<long,int>::iterator mitr2;
    //mitr=allbound.begin();mitr2=mitr;mitr2++;
    //int nt=0;
    vector<long> startcod,lencod;
    for(int i=0;i<allrange.size();i++){
      startcod.push_back(allrange[i].first);
      lencod.push_back(allrange[i].second-allrange[i].first+1);
    }
    //while(mitr2!=allbound.end()){
    //  if(validsegs[nt]==1){
    //    startcod.push_back(mitr->first);
    //    lencod.push_back(mitr2->first-mitr->first);
    //  }
    //  mitr++;
    //  mitr2++;
    //  nt++;
    //}

    boundfs<<chrname<<"\t"<<startr-1<<"\t"<<endr<<"\t";//chrom, chromstart, chromend
    boundfs<<"Inst"<<instanceid<<"_range"<<"\t";//4th, name
    boundfs<<score<<"\t";//5th, score
    boundfs<<bd<<"\t";//6th,orientation
    boundfs<<startr-1<<"\t"<<endr<<"\t";//7,8
    boundfs<<"0,255,0\t";//9th, color

    boundfs<<nbks<<"\t";//10th, blockcounts
    //11th, sizes
      //boundfs<<(range.second-range.first)<<"\t";
    for(int i=0;i<lencod.size();i++){
      boundfs<<lencod[i];
      if(i<lencod.size()-1)boundfs<<",";
      else boundfs<<"\t";
    }
    //12th, starts
      //boundfs<<0;
    for(int i=0;i<startcod.size();i++){
      boundfs<<startcod[i]-startcod[0];
      if(i<startcod.size()-1)boundfs<<",";
      else boundfs<<"\t";
    }
    boundfs<<endl;
    
    if(plotexons==true && startcod.size()>1){
      //plot exons boundary
      for(int n=0;n<startcod.size();n++){
        //only plot adjacent segs
        bool needplot=false;
        if(n+1<startcod.size() && startcod[n]+lencod[n]==startcod[n+1])needplot=true;
        if(n>0 && startcod[n-1]+lencod[n-1]==startcod[n])needplot=true;
        if(needplot==false)continue;
        boundfs<<chrname<<"\t"<<startcod[n]-1<<"\t"<<startcod[n]-1+lencod[n]<<"\t";//chrom, chromstart, chromend
        boundfs<<"Inst"<<instanceid<<"_segs_"<<n<<"\t";//4th, name
        boundfs<<1<<"\t";//5th, score
        boundfs<<bd<<"\t";//6th,orientation
        boundfs<<startcod[n]-1<<"\t"<<startcod[n]-1+lencod[n]<<"\t";//7,8
        if(n%2==0)
          boundfs<<"255,255,0\t";//9th, color
        else
          boundfs<<"255,0,255\t";//9th, color

        boundfs<<1<<"\t";//10th, blockcounts
        //11th, sizes
        boundfs<<lencod[n]<<"\t";
        //12th, starts
        boundfs<<0;
        boundfs<<endl;
      }
    }
  }
}

/**
 Writing a cluster of reads to instance
 REFONLY is not yet implemented.
*/
int write2rangeandinstance(ReadGroup& rd,range_t& range){
  if(rd.size()<MIN_GRANGE_READ_CNT )
    return -1;
  
  vector<ReadGroup> subg;
  // if(VERBOSE==1){
  //   cout<<"Coming in, Reads: "<<rd.size()<<", range: ["<<range.first<<","<<range.second<<"], ReadGroup range: ["<<rd.getRange().first<<","<<rd.getRange().second<<"]"<<endl;
  // }
  writeCoverage(rd);
  //split
  if(FIXRANGE){
    ReadGroup crn=C_RANGE_ANNO.getReadGroup(rd.getChr(),range);
    if(crn.size()==0){
      int r=rd.getDirSum();
      subg.push_back(rd);//do not try to split the large readgroup into smaller ones
      subg.back().setDir(r);
    }
    else{
      //setting up different range sets
      vector<RangeSet> annors(crn.size());
      vpos_t & ans=crn.s();
      vpos_t & ane=crn.e();
      for(int i=0;i<crn.size();i++){
        annors[i].add(ans[i],ane[i]);
      }
      //split by current range set
      rd.splitByRangeSet(subg,annors);
      //use the annotation's direction
      for(int i=0;i<crn.size();i++){
        int curdir=crn.getDirection()[i];
        //cout<<"CURRENTDIR: "<<curdir<<endl;
        subg[i].setDir(curdir);
      }
    }
  }
  else{
    //remove too long reads; this parameter should be large for human/mouse organism
    rd.removeTooLongReads(MAX_READ_SPAN,MAX_READ_OVER);
    rd.splitByRangeSet(subg,MIN_GRANGE_DISTANCE);

    //determine the direction
    for(int i=0;i<subg.size();i++){
      int r;
      if(STRANDED_RNASEQ!=0){
        r=STRANDED_RNASEQ;
      }else{
        r=subg[i].getDirSum();
      }
      subg[i].setDir(r);
    }
  }

  //enumerate all boundary
  for(int i=0;i<subg.size();i++){
    writeannotation(subg[i],n_INST);
    n_INST++;
    // cout<<"Writing Instance "<<n_INST<<",range "<<subg[i].getChr()<<":"<<subg[i].getRange().first<<"-"<<subg[i].getRange().second
    //  <<" Reads: "<<subg[i].size()<<",paired:"<<subg[i].peSize()<<endl;
    if(FIXBOUND){
      //with boundary fixed
      vector<range_t> jrange;
      C_JUNCTION.getAllRanges(subg[i].getChr(),subg[i].getRange(),jrange);
      subg[i].setupBound(jrange);
      subg[i].calculateType();
      subg[i].getCvgStatistics();
    }else{
      if(subg[i].size()<MIN_GRANGE_READ_CNT)continue;
      subg[i].calculateBound(true,DEFAULT_MIN_JUNCTION,MIN_GRANGE_DISTANCE);
      subg[i].calculateType();
      subg[i].removeWeakSegs(MIN_CVG_FRACTION);
      subg[i].calculateValidSegs();
      if(subg[i].validSize()<MIN_GRANGE_READ_CNT)continue;
    }
    vector<range_t> segrange=subg[i].getSegs();
    writeboundfs(subg[i].getRange(),n_INST,subg[i].getChr(),segrange,subg[i].size(),subg[i].getDir());
    fileSplitter->startWritingInstance(); // 2015-01-25
    iofs<<"Instance\t"<<n_INST<<endl;
    subg[i].toStream(iofs);
    fileSplitter->endWritingInstance(/* info ? */); // 2015-01-25
  }
  
  return 0;
};



/**
 * Write the constructed instance to output
 * All coordinates should be 1 base.
 * Notice that rpoolstart and rpoolend are inclusive
 * Replaced by Instance<< function.
 */
void write2Instance(range_t& currentrange,
           vector<vector<long> >&rpoolstart, vector<vector<long> >&rpoolend,
    int appearpereads,
    int instanceid,
    string chrname,
    int readlen,
    vector<vector<long> >& annostart, vector<vector<long> >& annoend,   //annotations within currentrange
    bool wantcoverage,
    bool refonly
    ){
  ofstream &out=iofs;
  //cout<<"write2Inst, rpoolsize: "<<rpoolstart.size()<<" and "<<rpoolend.size()<<endl;
  //go over the reads, find out all junctions, and separate the genome into several different segments
  map<long,int> allbound;
  //the type of every boundary; 0 for junction, 1 for low coverage
  map<long,int> boundtype;
  //use the length of the 1st read as the length of the read
  //int readlen=0;
  //for(int i=0;i<rpoolstart[0].size();i++){
  //  readlen+=rpoolend[0][i]-rpoolstart[0][i]+1;
  //}

  removeoutofrangereads(currentrange,rpoolstart,rpoolend,appearpereads);
  if(rpoolstart.size()==0){
    // cout<<"Instance "<<instanceid<<"("<<chrname<<":"<<currentrange.first<<"-"<<currentrange.second<<") is empty. Ignore it.\n";
    return;
  }

  map<long,int>::iterator it, it2;

  //this variable is used to collect statistics of different segments. The following values in vector<double> are: 0, max cvg, 1, leftmost cvg, 2, rightmost cvg, 3, fraction of 0 coverage
  map<long,vector<double> >cvgstat;
  //incorporate boundaries from annotation
  //notice that start and end boundaries are also added here
  //use this at caution!
  //if(false)
  for(int i=0;i<annostart.size();i++){
    if(annostart[i].size()>0){
      for(int j=0;j<annostart[i].size();j++) {
        int val=annostart[i][j];
        if(allbound.count(val)==0)
          allbound[val]=2;//mark as special junctions
      }
      for(int j=0;j<annostart[i].size();j++) {
        int val=annoend[i][j]+1;
        if(allbound.count(val)==0)
          allbound[val]=2;//mark as special junctions from the annotation
      }
    }
  }
  //get the boundary
  getallbound(allbound,currentrange,rpoolstart,rpoolend,cvgstat,refonly);
  boundtype=allbound;
  //reset all values of allbound to 0
  for(it=allbound.begin();it!=allbound.end();it++)
    it->second=0;

  int nExons=allbound.size()-1;

  //consider ref annotations as special reads, and get their type
  map<vector<int>,int> reftype, reforder;
  vector<vector<int> > refseriestype;
  vector<int> allreftype;
  //get all types of the annotations
  map<long,int>refbound=allbound;
  //njuncs is used to store the number of junction reads falling onto one segment
  //emptyjuncs is a dummy variable
  map<long,int> njuncs,emptyjuncs;
  getalltypeorders(reftype,reforder,refseriestype,allreftype,refbound,annostart,annoend,emptyjuncs);

  //after that, count the number of reads falling into each segments, and calculate the fingerprint for each read
  map<vector<int>,int> contenttype;
  //and, assign a type number for it. We will write these types in the order of this number
  map<vector<int>,int> typeorder;
  //save the type according to its order, from 0 to n
  vector<vector<int> > seriestype;
  //save the type of all the reads
  vector<int>  allreadtype;

  //get all types and their orders
  getalltypeorders(contenttype,typeorder,seriestype,allreadtype,allbound,rpoolstart,rpoolend,njuncs);

  //get seg lengths
  vector<int> seglengths;
  it=allbound.begin();it2=it;it2++;
  while(it2!=allbound.end()){
    seglengths.push_back(it2->first-it->first);
    it++;it2++;
  }
  

  //check segs with non-zero reads
  vector<int> validsegs(nExons,0);
  int nvalidsegs=0;
  int validid=0;
  
  it=allbound.begin();it2=it;it2++;
  //also, iterate refbound too, retain segments which appear in ref isoforms
  map<long,int>::iterator rit=refbound.begin();

  while(it2!=allbound.end()){
    if(it->second!=0 || rit->second!=0){
      if(seglengths[validid]<readlen 
        &&( boundtype[it->first]==1 || boundtype[it2->first]==1 )
        && rit->second==0 // a new condition that this segment does not appear in ref isoforms
        && njuncs[it->first]==0 //a new condition that this segment contains no junction reads
      ){ 
          //encounter a short range, check if the range is due to the coverage cut

      }
      else{
        validsegs[validid]=1;
        nvalidsegs++;
      }
    }
    validid++;
    it++;
    it2++;
    rit++;
  }
  if(nvalidsegs==0){
    // cout<<"0 segs, return..."<<endl;
    return;//don't write instance with 0 segs
  }


  //annotation: write boundary of instances
  //writeboundfs(currentrange,instanceid,chrname,allbound,validsegs);



  //finally, write this data to output
  out<<"Instance\t"<<instanceid<<endl;
  out<<"Boundary\t"<<chrname<<"\t"<<currentrange.first<<"\t"<<currentrange.second<<endl;
  out<<"ReadLen\t"<<readlen<<endl;


  out<<"Segs\t"<<nvalidsegs<<endl;
  // cout<<"Segs\t"<<nvalidsegs<<"\tBoundary\t"<<chrname<<":"<<currentrange.first<<"-"<<currentrange.second<<endl;
  //start and end of segments
  it=allbound.begin();it2=it;it2++;
  validid=0;
  while(it2!=allbound.end()){
    if(validsegs[validid]==1){
      out<<it->first<<"\t"<<it2->first-1<<"\t"<<(it2->first-it->first)<<"\t"<<it->second;//range and length
      vector<double>& currentstat=cvgstat[it->first];
      if(currentstat.size()==0){
        //cerr<<"Error: not enough fields in cvg stat.\n";
      }
      else{
        //statistics, including: max cvg, leftmost cvg, rightmost cvg, zero fraction cvg, and average cvg
        for(int i=0;i<currentstat.size();i++)
          out<<"\t"<<currentstat[i];
        //out<<"\t"<<currentstat[0]<<"\t"<<currentstat[1]<<"\t"<<currentstat[2]<<"\t"<<currentstat[3];
      }
      out<<endl;//[it->first, it2->first-1]
    }
    it++;
    it2++;
    validid++;
  }
  //ref 
  out<<"Refs\t"<<refseriestype.size()<<endl;
  for(int i=0;i<refseriestype.size();i++){
    vector<int> exonbin(nExons,0);
    for(int j=0;j<refseriestype[i].size();j++) exonbin[refseriestype[i][j]]=1;
    for(int j=0;j<nExons;j++){
      if(validsegs[j]==1)
        out<<exonbin[j]<<" ";
      // else{
        // if(exonbin[j]==1){
          // cerr<<"Error: incorrect valid segs...\n";
        // }
      // }
    }
    out<<endl;
  }

  //reads
  out<<"Reads\t";
  if(appearpereads==1)
    out<<rpoolstart.size()/2<<endl;
  else if(appearpereads==0)
    out<<rpoolstart.size()<<endl;
  // else
    // cerr<<"Error: appearpereads is "<<appearpereads<<endl;
  
  //SG types
  //we need to re-arrange SGType orders, if some SGType contains invalid segments
  //reshuforder is a map storing <pre-order,post-order>.
  map<int,int>reshuforder;
  int reshufid=0;
  vector<int> validsgtypes(seriestype.size(),0);
  
  for(int i=0;i<seriestype.size();i++){
    bool hasinvalidsegs=false;
    for(int j=0;j<seriestype[i].size();j++){
      if(validsegs[seriestype[i][j]]==0){
        //ignore this type and jump
        hasinvalidsegs=true;
        break;
      }

    }
    if(hasinvalidsegs){
      validsgtypes[i]=0;
    }
    else{
      validsgtypes[i]=1;
      reshuforder[i]=reshufid;
      reshufid++;
    }

  }
  out<<"SGTypes\t"<<reshufid<<endl;//adjusted # of SGTypes
  for(int i=0;i<seriestype.size();i++){
    if (validsgtypes[i]==0)continue;
    vector<int> exonbin(nExons,0);
    for(int j=0;j<seriestype[i].size();j++){
             exonbin[seriestype[i][j]]=1; 
    }
    for(int j=0;j<nExons;j++){
      if(validsegs[j]==1)
               out<<exonbin[j]<<" ";
    }
    out<<contenttype[seriestype[i]]<<endl; //# of reads falling onto this type
  }
  //PE types
  if(appearpereads==0){
    out<<"PETypes\t0"<<endl;
  }
  else{
    //group the pe reads by orders
    map< pair<int,int>, vector<long> > allpetypes;
    pair<int,int> thispair=make_pair<int,int>(0,0);
    for(int i=0;i<rpoolstart.size();i+=2){
      //ATTENTION: change the order from 0 to 1 here.
      //out<<typeorder[allreadtype[i]]+1<<"\t"<<typeorder[allreadtype[i+1]]+1<<"\t"<<(rpoolstart[i+1][0]-rpoolend[i][rpoolend[i].size()-1])<<endl;
      //allreadtype stores the type order
      int fod=allreadtype[i],sod=allreadtype[i+1];
      if(validsgtypes[fod]==0 || validsgtypes[sod]==0)continue;
      // if(reshuforder.count(fod)==0 || reshuforder.count(sod)==0){
        // cerr<<"Error: incorrect reshuffle count id!!!\n";
      // }
      thispair.first=reshuforder[fod]+1;
      thispair.second=reshuforder[sod]+1;
      allpetypes[thispair].push_back((rpoolstart[i+1].front()-rpoolend[i].back()));
    }
    out<<"PETypes\t"<<rpoolstart.size()/2<<"\t"<<allpetypes.size()<<endl;
    map<pair<int,int>, vector<long> >::iterator apitr=allpetypes.begin();
    while(apitr!=allpetypes.end()){
      out<<(apitr->first).first<<"\t"<<(apitr->first).second<<"\t"<<(apitr->second).size()<<endl;
      for(int i=0;i<(apitr->second).size();i++){
        out<<(apitr->second)[i]<<" ";
      }
      out<<endl;
      apitr++;
    }
  }
  if(wantcoverage){
    //the individual base coverage should count # of bases having # of reads.
    //the information is stored in 2*x lines, where x is the number of SG Types.
    //For every two lines, the first line is the number of records of the second line
    //the second line stores pairs (m,n), where m is the # of reads falling into one base, 
    //and n is the number of bases having m reads on it.
    //NOTICE: for paired-end reads, use the 1st nt of the 2nd read as start position 
    out<<"Coverage\t"<<reshufid<<"\t"<<rpoolstart.size()<<endl;
    vector<vector<long> > allstartrange(seriestype.size());
    int step=1;
    int pst=0;
    if(appearpereads==1){step=2;pst=1;}
    for(int i=pst;i<rpoolstart.size();i+=step){
      allstartrange[allreadtype[i]].push_back(rpoolstart[i][0]); //storing starting position
    }
    //count occurences
    map<int,int> tcounts;
    for(int i=0;i<allstartrange.size();i++){
      if (validsgtypes[i]==0)continue;
      tcounts.clear();
      long prevp=-1;
      int ctr=1;
      for(int j=0;j<allstartrange[i].size();j++){
        if(allstartrange[i][j]!=prevp){
          if(prevp!=-1)tcounts[ctr]++;
          ctr=1;
          if(allstartrange[i][j]<prevp){
            // cerr<<"Error: in output individual coverage, the read must be sorted."<<endl;
            // cerr<<"Start range: "<<j<<" :"<<allstartrange[i][j]<<"-"<<prevp<<endl;
            // for(int k=0;k<rpoolstart.size();k++){
              // cerr<<rpoolstart[k][0]<<" ";
            // }
            // cerr<<endl;
            // exit(-1);
          }
        }
        else
          ctr++;
        prevp=allstartrange[i][j];
      }
      if(prevp!=-1)tcounts[ctr]++;
      out<<reshuforder[i]<<"\t"<<tcounts.size()<<endl;
      for(map<int,int>::iterator tit=tcounts.begin();tit!=tcounts.end();tit++){
        out<<tit->first<<","<<tit->second<<"\t";
      }
      out<<endl;
    }

  }//end OUTPUT_INDIVIDUAL_COVERAGE
  
  

}

/* Read gene range file.
Only the first 4 fields for each line is used.
Other fields are ignored. */
int readgenerangefile(string filename, GeneRange & gr){
  gr.clear();
  ifstream ifs(filename.c_str());
  if(!ifs.is_open()){
    // cerr<<"Error opening gene range file "<<filename<<"."<<endl;
    return -1;
  }
  //read gene range file
  string oneline;
  int nread=0;
  while(true){
    getline(ifs,oneline);
    if(oneline[0]=='#')continue;//ignore comment lines
    if(ifs.eof())break;
    stringstream ss(oneline);
    int id;
    string chr;
    char ori;
    long starti,endi;
    //ss>>id>>chr>>ori>>starti>>endi;
    ss>>chr>>starti>>endi>>ori;
    gr.push_back(chr,make_range_t(starti,endi));
    nread++;
  }
  // cout<<"Reading "<<gr.size()<<" records.\n";
  ifs.close();
  return 0;
}

/*read existing annotation file (using .bed format), and incorporate into instance, and cluster them at the same time
  Notice that the start coordinate of bed file is 0-based, and the end coordinate is 1-based.
  so to be compatible with SAM and our program (SAM  uses 1-based, and the range in this program is 1-based, and inclusive [start,end]),
  we need to handle the coordinate carefully
*/
int parseannobed(string infile, Annotation& anno ){
  ifstream annobedifs;
  annobedifs.open(infile.c_str());
  if(!annobedifs.is_open()){
    // cerr<<"Error opening annotation bed file "<<infile<<endl;
    return -1;
  }

  string oneline;
  int linecount=0;
  while(true){
    getline(annobedifs,oneline);
    if(annobedifs.eof())break;
    linecount++;
    //parse bed file
    vector<string>fields;
    int current=-1,next=-1;
    while((next=oneline.find('\t',current+1))!=-1){
      fields.push_back(oneline.substr(current+1,next-current-1));
      current=next;
    }
    fields.push_back(oneline.substr(current+1));
    //it should be 12 fields
    if(fields.size()!=12){
      // cerr<<"Error: not enough field (12) at line "<<linecount<<", size="<<fields.size()<<endl;
      continue;
    }
    long startp,endp;
    int nsegments;
    stringstream ss1(fields[1]), ss2(fields[2]), ss3(fields[9]),ss4(fields[10]),ss5(fields[11]);
    ss1>>startp;ss2>>endp;ss3>>nsegments;
    //convert 0-base to 1-base;
    startp+=1;
    endp+=1;//notice that here endp is exclusive in 1-base
    //parse the 11th and 12th field
    pos_t startpos(nsegments),endpos(nsegments);
    char tmpc;
    for(int i=0;i<nsegments;i++){
      long len,relstart;
      ss4>>len;ss5>>relstart;
      ss4>>tmpc;ss5>>tmpc;
      startpos[i]=startp+relstart;
      endpos[i]=startp+relstart+len-1;//inclusive
    }
    if(endpos.back()!=endp-1){
      // cerr<<"Error: incorrect field 10 and field 11 at line "<<linecount<<endl;
      continue;
    }
    //direction
    int dir=0;
    if(fields[5]=="+") dir=1;
    if(fields[5]=="-") dir=-1;

    //save to record
    if(anno.add(fields[0],startpos,endpos,dir,fields[3])==-1){
      // cerr<<"Error: the BED file is not sorted at line "<<linecount<<endl;
      return -1;
    }
  }//end while(true)

  annobedifs.close();
  //add one dummy record to make annotation cluster the remaining records
  anno.cluster();

  //some statistics
  // cout<<"Annotation chromosomes: ";
  vector<string> alln=anno.getChrom();
  // for(int i=0;i<alln.size();i++) cout<<alln[i]<<" "; cout<<endl;


  // cout<<"Clusters:\n";
  // cout<<"Total clusters: "<<anno.getNCluster()<<", total ref in clusters: "<<anno.getNAnnotation()<<endl;

  return 0;
}



