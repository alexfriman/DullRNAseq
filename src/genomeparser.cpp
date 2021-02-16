#include <fstream>
#include <cstring>
#include <vector> 
#include "genomeparser.h"
using namespace std;



string readGenome(char* fname){
	fstream f;
	printf("Reading fastA: %s\n",fname);
	f.open(fname, fstream::in);
	string wholedata;
	string line;
	while(getline(f, line)) {
		if (line.at(0)=='>'){
			line="^"+line+"\n"; //Marking beginning of the genes with ^> marker
		}
        wholedata+=line;
    }
	f.close();
	return wholedata;
}

vector<gene> parseGenome(char* genome, int rank){
	//parse into lines
	vector<gene> allthegenes;
	char *token;
	string tmpString;
	size_t found;
	const char* s = "^>";
	token = strtok(genome, s);
	while( token != NULL ) {
	  tmpString=token; //assignment constructor
	  found=tmpString.find("\n");
	  if (found!=string::npos){
	  allthegenes.push_back({tmpString.substr(0,found),tmpString.substr(found)});
	  }else{
		tmpString+=strtok(NULL, s);
		found=tmpString.find("\n");
		allthegenes.push_back({tmpString.substr(0,found),tmpString.substr(found)});
	 }
	  
      token = strtok(NULL, s);
	}
   printf("Rank %d found %lu genes in the fastA\n",rank,allthegenes.size());
   return allthegenes;
}