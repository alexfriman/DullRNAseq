#include <vector> 
#include "aligner.h"
//#include "genomeparser.h"
using namespace std;
void hitCounter(vector<gene> Genome, vector<string> Reads, unsigned long int* counts, int rank){
	size_t hitPosition;
	for (long unsigned int gn=0; gn<Genome.size(); gn++){
		for (long unsigned int sn=0; sn<Reads.size(); sn++){
			if (Genome[gn].seq.size()>=Reads[sn].size()){
				hitPosition=Genome[gn].seq.find(Reads[sn]);
				if (hitPosition!=string::npos){ 
					counts[gn]++;
				}
			}
			
		}
		if (gn%500==0) printf("Alignment: Rank %d processed %lu/%lu genes\n", rank, gn, Genome.size());
	}
	
}