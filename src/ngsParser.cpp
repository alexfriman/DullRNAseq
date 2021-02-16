#include <fstream>
#include <cstring>
#include <vector> 
//#include <math.h> 
#include "ngsParser.h"
using namespace std;
void readNSG(char* fname, string* portions, int clusterSize){
	unsigned long int readnum=0;
	int phase=0;
	fstream f;
	printf("Using FASTQ: %s \n",fname);
	f.open(fname, fstream::in);
	string tmpLine="";
	string line;
	while(getline(f, line)) {
		if (phase==0) {
			if (line.at(0)!='@'){
				printf("PROBLEM line: %s\n",line.c_str());
				printf("Something is wrong with ngsParser phase\n");
				abort();
			}
			//saving data
			portions[readnum%clusterSize]+=tmpLine;
			tmpLine="";
			readnum++;
		}
		tmpLine+="\n"+line;
		phase=(phase+1)%4;
    }
	f.close();
}

/*float OldIlluminaErrorCodes(char symbol){
	int Q=int(symbol)-64;
	float coef=pow(10,-0.1*Q);
	float p=coef/(1+coef);
	return p;
}*/

vector<string> parseNGS(string rawreads, float acceptedError, int rank){
	vector<string> reads;
	int phase=0;
	string buffer;
	string seqbuffer;
	string scorebuffer;
	string goodSeqbuffer="";
	size_t Nposition;
	size_t NSkip;
	size_t skip=0;
	size_t nextline=rawreads.find("\n",0);
	while  (nextline!=string::npos){
		if (phase==1){ 
			buffer=rawreads.substr(skip,nextline-skip);
			seqbuffer=buffer;
			NSkip=0;
			Nposition=seqbuffer.find("N");
			
			while (Nposition!=string::npos){
				//there is N in the sequence
				if (Nposition-NSkip>MIN_SEQ) reads.push_back(seqbuffer.substr(NSkip, Nposition));
				NSkip=Nposition+1;
				Nposition=seqbuffer.find("N",NSkip);
			}
			if (seqbuffer.size()-NSkip>MIN_SEQ) reads.push_back(seqbuffer.substr(NSkip));
		}
		skip=nextline+1;
		nextline=rawreads.find("\n",skip);
		phase=(phase+1)%4;
	}
	printf("Rank %d processed fastQ of %lu symbols, found %lu reads\n",rank,rawreads.size(),reads.size());
	return reads;
}