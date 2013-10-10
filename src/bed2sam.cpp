/**
 * This program is used to convert .bed file (e.g., flux output) into .sam file (e.g., tophat output).
 * ATTENTION: the name of the read must follow some patterns. It may be flux read conventions, or, like 'uc007aev.1_e_11652_chr1_3638394'
 */

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <vector>
#include <set>
#include <ctime>
#include <cstdlib>

using namespace std;


bool issingle=true;

/*
 * Sample rate
 */
double samplerate=1.0;

bool needsort=false;

long nrecord=0;

int fixreadlen=-1;

bool nofluxname=false;

//record iso orientation
map<string,char>isoorientation;

//one bed record
struct bedrec{
	string name; //read name
	
	long start; //mapping start position
	long end; //mapping end position
	bool isjunction; //If this is a junction read
	string cigar;//CIGAR string
	char ori; //oritentation
	int readlen; //length of the read
	long pairstart;//for pair-end reads only, the pair-end starting position
	char mateori; //for pair-end reads only, strand of the pair
	bool isfirst; //for pair-end reads only, indicate if this read is the first read of the pair or not
};

//Parse 10, 11 and 12th output of .bed file, and use it to calculate CIGAR string.
string getcigar(int nblock,string blocksizes,string blockstarts, int& readlen){
	stringstream bz(blocksizes);
	stringstream bst(blockstarts);
	char tmp; //used to separate ','
	stringstream cigar;
	int prevstart=-1;
	int prevsize=-1;
	readlen=0;
	for(int i=0;i<nblock;i++){
		int onesize,onestart;
		bz>>onesize;
		bst>>onestart;
		if(i!=nblock-1){
			bz>>tmp;//parse ','
			bst>>tmp;
		}
		if(i>=1){
			cigar<<(onestart-prevstart-prevsize)<<"N";
		}
		cigar<<onesize<<"M";
		prevstart=onestart;
		prevsize=onesize;
		readlen+=onesize;

	}
	return cigar.str();
}

/**
 * Get the flag field for SAM output
 */
int getflag(bedrec& onerec){
	int flag=0;

	if(issingle==false) {
		flag |= 0x0001; //is paired
		flag |= 0x0002; //is mapped in a proper pair
		if(onerec.mateori=='-')
			flag |=0x0020; //strand of the mate
		if(onerec.isfirst)
			flag |=0x0040; //1st read
		else
			flag |=0x0080; //2nd read
	}
	if(onerec.ori=='-'){
		flag |= 0x0010;	//strand
	}

	return flag;
	
}

/**
	For flux name convention, get the transcript id from read name. transcript id lies between the 2nd and the 3rd :.
*/
string getisonamefromname(string inname){
	int colonstart=-1;
	int colonend=-1;
	int coloncount=0;
	if(nofluxname){
		int bs=inname.find('_');
		return inname.substr(0,bs);
	}
	for(int i=0;i<inname.size();i++){
		if(inname[i]==':'){
			coloncount++;
			if(coloncount==2)colonstart=i+1;
			if(coloncount==3){colonend=i;break;}
		}
	}
	return inname.substr(colonstart,colonend-colonstart);
}


void writeonerecord2sam(ofstream &out, bedrec& rec, string chrname){
		static int errcount=0;
		if(rec.start<=0)return;
		out<<rec.name<<"\t"; //1st ,QNAME
		int flag=getflag(rec);//2nd, flag
		out<<flag<<"\t";
		out<<chrname<<"\t";//3rd, chromosome name
		out<<((rec).start)<<"\t"; //4th, start position, 1 based
		int score=255;
		out<<score<<"\t";	//5th, mapping quality
		out<<(rec).cigar<<"\t"; //6th, CIGAR string
		if(issingle){
			out<<"*\t0\t"; //for single-end , 7th and 8th should be star and 0.
		}
		else{
			out<<"=\t"<<(rec).pairstart<<"\t"; //for pair-end, 7th is '=' and 8th is the pair-end mapping position.
		}
		out<<0<<"\t"; //9th, ISIZE, 0 right now
		out<<"*"<<"\t"; //10th, SEQ, unavailable right now
		out<<"*"<<"\t"; //11th, QUAL, unavailable right now
		out<<"NM:i:0\t";
		if((rec).isjunction){
			char xs;
			//THIS IS A BUG VERSION:
			//According to cufflinks, this field must exist for junction reads. This indicates which strand the RNA that produced this read came from.
			//If this is the first read (0x0040) and forward (0x0020, mate reverse), or this is the second read (0x0080) and reverse (0x0010).
			//if( ( (flag & 0x0040 ) & (flag & 0x0020) ) || ( (flag & 0x0080) & (flag & 0x0010)))
			//	xs='+';
			//else
			//	xs='-';
			//out<<"XS:A:"<<xs<<"\t"<<"NS:i:0";
			string isoname=getisonamefromname(rec.name);
			if(isoorientation.count(isoname)==0){
				// if(errcount<10)cerr<<"Warning: no records "<<isoname<<" found in isoorientation.txt.\n";
				errcount++;
			}
			else{
				out<<"XS:A:"<<isoorientation[isoname]<<"\t"<<"NS:i:0";
			}

		}
		out<<"\n";

		nrecord++;
		// if(nrecord%1000000==1)cout<<"Writing "<<nrecord<<" records...\n";


}

/**
 * Write results to .sam file
 */
void write2samfile(ofstream & out, multimap<long, bedrec>& sortmap, string chrname){
	if(needsort==false)return;
	// cout<<"Writing records of chromosome "<<chrname<<", records: "<<sortmap.size()<<"...\n";
	multimap<long,bedrec>::iterator it;
	for(it=sortmap.begin();it!=sortmap.end();it++){
		if(fixreadlen==-1 || fixreadlen==(it->second).readlen)
			writeonerecord2sam(out,it->second,chrname);
	}
}
/**
 * Parse one line, store the result to onerec ,and return the chromosome name
 */
string parseoneline(string oneline, bedrec& onerec){
		stringstream ss(oneline);
		string chrname;
		ss>>chrname; //1st, chromosome name
		long lstart,lend;
		ss>>lstart;
		ss>>lend; //2nd, 3rd, mapping start and end
		string readname;
		ss>>readname; //4th, read id
		int tmp;
		ss>>tmp; //5th, not used by flux
		char ori;
		ss>>ori; //6th, orientation, should be + or -
		string tmps;
		ss>>tmps; //jump 7th, 8th, and 9th field
		ss>>tmps; //jump 7th, 8th, and 9th field
		ss>>tmps; //jump 7th, 8th, and 9th field
		int nblock;
		ss>>nblock; //10th, the number of blocks one read maps to
		string blocksizes;
		ss>>blocksizes; //11th, block sizes
		string blockstarts;
		ss>>blockstarts; //12th, block starts

		//write to bedrec
		onerec.name=readname;
		onerec.start=lstart+1; //ATTENTION: here it changes from 0-based to 1-based
		onerec.end=lend;
		onerec.ori=ori;
		onerec.cigar=getcigar(nblock,blocksizes,blockstarts,onerec.readlen);
		if(nblock>1)
			onerec.isjunction=true;
		else
			onerec.isjunction=false;

		return chrname;
}

/*
 * Read input out stream for chromosome desired chromosome, until it encounter a new chromosome name.
 * If the current chromosome name is desiredchr, it will place all the record into the sortmap.
 * If a new chromosome is detected, it will stop and go back to the place before such chromosome is detected.
 */

int readinput(ifstream & ifs,ofstream& ofs,
		 multimap<long, bedrec>& sortmap,
	       	string desiredchr, string& newchr){
	streampos pos;
	//Read .bed files
	bedrec onerec, onerec2;
	onerec.pairstart=0;
	onerec2.pairstart=0;
	string currentchrname="";
	while(true){
		string oneline,oneline2;
		pos=ifs.tellg();//remember the position before
		getline(ifs,oneline);
		if(ifs.eof()){
			newchr="";
			break;
		}
		if(issingle){
			//for single-end
			//parse .bed inputs
			string tchr;
			tchr=parseoneline(oneline,onerec);
			if(tchr=="polyA")continue;
			currentchrname=tchr;
			//detect chromosome change, write to file
			if(currentchrname!=desiredchr){
				//cout<<"Switching to chromosome "<<currentchrname<<"...\n";
				newchr=currentchrname; //remember next possible chromosome
				ifs.seekg(pos); //retreat
				break;
			}
			else{
				//insert to map
				if(rand()*1.0/RAND_MAX<=samplerate){
					if(needsort)
						sortmap.insert(pair<long, bedrec>(onerec.start,onerec));
					else{
						if(fixreadlen==-1 || onerec.readlen==fixreadlen)
							writeonerecord2sam(ofs,onerec,currentchrname);
					}
				}
			}
		}
		else{
			//for pair-end read, read two line at the same time.
			getline(ifs,oneline2);
			//parse
			string tchr;
			tchr=parseoneline(oneline,onerec);
			
			string anotherchr;
			anotherchr=parseoneline(oneline2,onerec2);
			if(tchr=="polyA" || anotherchr=="polyA")continue;
			currentchrname=tchr;
			if(anotherchr!=currentchrname){
				// cerr<<"Error: pair-end reads have different chromosome name for "<<currentchrname<<" and "<<anotherchr<<endl;
			}

			//exchange start position
			onerec.pairstart=onerec2.start;
			onerec2.pairstart=onerec.start;

			onerec.isfirst=true;
			onerec2.isfirst=false;

			onerec.mateori=onerec2.ori;
			onerec2.mateori=onerec.ori;

			//detect chromosome change, write to file
			if(currentchrname!=desiredchr){
				//cout<<"Switching to chromosome "<<currentchrname<<"...\n";
				newchr=currentchrname;
				ifs.seekg(pos);
				break;
			}
			else{
				//insert to map
				if(rand()*1.0/RAND_MAX<=samplerate){
					if(needsort){
						sortmap.insert(pair<long, bedrec>(onerec.start,onerec));
						sortmap.insert(pair<long, bedrec>(onerec2.start,onerec2));
					}
					else{
						if(fixreadlen==-1 || (onerec.readlen==fixreadlen && onerec2.readlen==fixreadlen)){
							writeonerecord2sam(ofs,onerec,currentchrname);
							writeonerecord2sam(ofs,onerec2,currentchrname);
						}
					}
				}
			}


		}//end if(issingle)
	}
	return 0;

}


void printhelp(){
		// cout<<"This program is used to convert .bed files (e.g., flux output) into .sam file (e.g., tophat output).\n";
		// cout<<"Usage: bed2sam {OPTIONS} [.bed file1] [.bed file2] ...\n";
		// cout<<"Options:\n";
		// cout<<"-s	The read is single-ended.\n";
		// cout<<"-p	The read is pair-ended.\n";
		// cout<<"-o [.sam file name] Output file name. Default is bed2samout.sam at current directory\n";
		// cout<<"-r rate	Sample rate. Randomly select part of the reads to output. The range of rate should be [0,1].\n";
		// cout<<"-n	Sort all records before writing to sam file. (Not recommended; use Unix's sort command instead).\n";
		// cout<<"-t [isoorentation file] Specify the location of isoorientation.txt.\n";
		// cout<<"-l [fixed read length] Only output records with fixed read length\n";
		//cout<<"-k	The read name is not FLUX-like, but like 'uc007aev.1_e_11652_chr1_3638394'.\n";
}

int main(int argc, char* argv[]){
	if(argc<2){
		printhelp();
		// cerr<<"Error: no .bed file found.\n";
		return -1;
	}
	srand(time(NULL));
	string ofsstr="bed2samout.sam";
	vector<ifstream*> allinputs;
	string isoorientationfile="isoorientation.txt";
	for(int i=1;i<argc;i++){
		if(string(argv[i])=="-s"){
			issingle=true;
			// cout<<"Using single end reads...\n";
		}
		else if(string(argv[i])=="-p"){
			issingle=false;
			// cout<<"Using pair end reads...\n";
		}
		else if(string(argv[i])=="-n"){
			needsort=true;
			// cout<<"Warning: you turn the sort option on. This is not recommended. Use Unix's sort command instead.\n";
		}
		else if(string(argv[i])=="-o"){
			if(i+1<argc){
				ofsstr=argv[i+1];
				i+=1;
			}
			else{
				printhelp();
				return -1;
			}
		}
		else if(string(argv[i])=="-r"){
			stringstream ss;
			if(i+1<argc){
				ss<<(string(argv[i+1]));
				double nrate=-1;
				ss>>nrate;
				if(0<=nrate && nrate<=1)
					samplerate=nrate;
				else{
					// cerr<<"Error: sample rate should be [0,1].\n";
				}
				i+=1;
			}
			else{
				printhelp();
				return -1;
			}


		}
		else if(string(argv[i])=="-t"){
			isoorientationfile=string(argv[i+1]);
			// cout<<"isoorientation.txt path: "<<isoorientationfile<<endl;
			i+=1;
		}
		else if(string(argv[i])=="-l"){
			stringstream ss;
			if(i+1<argc){
				ss<<string(argv[i+1]);
				ss>>fixreadlen;
				// cout<<"Fix read length: "<<fixreadlen<<endl;
			}
			else{
				printhelp();
				return -1;
			}
			i+=1;
		}
		else if(string(argv[i])=="-k"){
			nofluxname=true;
		}
		else{
			//input sam files
			ifstream* ifsp=new ifstream(argv[i]);
			if(!ifsp->is_open()){
				// cerr<<"Error opening file "<<argv[i]<<"...."<<endl;
			}
			else{
				// cout<<"Using file "<<argv[i]<<"...\n";
				allinputs.push_back(ifsp);
			}

		}
	}

	//read iso orientation file
	ifstream oriifs(isoorientationfile.c_str());
	if(!oriifs.is_open()){
		// cerr<<"Error opening isoorientation.txt file "<<isoorientationfile<<endl;
	}
	else{
		string oneline;
		while(true){
			getline(oriifs,oneline);
			if(oriifs.eof())break;
			stringstream ss(oneline);
			string isoname;
			char ori;
			ss>>isoname>>ori;
			isoorientation[isoname]=ori;
			
		}
		// cout<<"Total isoorientation records: "<<isoorientation.size()<<endl;
	
	}
	
	
	ofstream ofs(ofsstr.c_str());
	if(!ofs.is_open()){
		// cerr<<"Error opening output file "<<ofsstr<<"...."<<endl;
		return -1;
	}
	// cout<<"Output file: "<<ofsstr<<endl;
	// cout<<"Sample rate: "<<samplerate<<endl;

	//writing head information for sam
	ofs<<"@HD	VN:1.0	SO:sorted\n";
	ofs<<"@PG	bed2sam	VN:1.0.13	CL:";
	for(int i=0;i<argc;i++)
		ofs<<argv[i]<<" ";
	ofs<<endl;

	multimap<long, bedrec> sortmap;
	
	string currentchrname="";
	string futurechrname;
	set<string> possiblechrnames;
	
	while(true){
		possiblechrnames.clear();
		for(int i=0;i<allinputs.size();i++){
			readinput( *(allinputs[i]), ofs, sortmap,  currentchrname,futurechrname);
			// cout<<"Next :"<<futurechrname<<endl;
			if(futurechrname!="")
				possiblechrnames.insert(futurechrname);
		}
		if(currentchrname!="")
			write2samfile(ofs,sortmap,currentchrname);
		if(possiblechrnames.size()==0){
			break; //all files finish
		}
		else{
			currentchrname=*(possiblechrnames.begin());
			//cout<<"Current chr name: "<<currentchrname<<endl;
			sortmap.clear();
		}

	}

	ofs.close();
	for(int i=0;i<allinputs.size();i++){
		(allinputs[i])->close();
		delete allinputs[i];
	}
	return 0;
}
