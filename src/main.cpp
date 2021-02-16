#include <stdio.h>
#include <string>
#include "mpi.h"
//#include "genomeparser.h"
#include "ngsParser.h"
#include "aligner.h"
#define GENOME_SIZE_TAG 0
#define GENOME_TAG 1
#define READS_SIZE_TAG 2
#define READS_TAG 3
#define READS_NUMBER_TAG 4
#define COUNTS_TAG 5
using namespace std;
int main(int argc, char *argv[])
{
	int rank;
	int size;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	if (argc<3){
		if (rank==0) printf("Usage:\nmpirun -np /threads_number/ DullRNAseq genome.fasta reads.fastq [reads2.fastq ...]\n");
		MPI_Finalize();
		return 0;
	}
	basic_string<char>::size_type actualsize;
	char *rawgenome;
	string myReads;
	string rawGenomestring;
	if (rank==0){
		rawGenomestring=readGenome(argv[1]);
		actualsize=rawGenomestring.size()+1;
	}
	MPI_Bcast(&actualsize,1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
	rawgenome=new char[actualsize];
	if (rank==0) strcpy(rawgenome, rawGenomestring.c_str());
	rawGenomestring.erase();
	MPI_Bcast(rawgenome, actualsize, MPI_CHAR, 0, MPI_COMM_WORLD);		
	printf("Rank %d, size of fastA %lu\n",rank, actualsize);
	vector<gene> allGenes=parseGenome( rawgenome, rank);
	delete rawgenome;
	string::size_type Readsbuffersize;
	if (rank==0){
			string *allReads = new string[size];
			for (int i=0; i<size;i++) allReads[i]="";
			for (int an=2; an<argc; an++) readNSG(argv[an], allReads, size);
			for (int i=0; i<size;i++) allReads[i]=allReads[i].substr(1);//removing the first \n
			myReads=allReads[0];
			//allReads[0].erase();
			for (int i=1; i<size;i++){ 
				//sending length of the reads	
				Readsbuffersize=allReads[i].size();
				MPI_Send(&Readsbuffersize,1, MPI_UNSIGNED_LONG,  i,READS_SIZE_TAG, MPI_COMM_WORLD);
			}
			for (int i=1; i<size;i++){ 
				Readsbuffersize=allReads[i].size();
				char* sendReadsBuffer = new char[Readsbuffersize+1];
				strcpy(sendReadsBuffer, allReads[i].c_str());
				MPI_Send(sendReadsBuffer,Readsbuffersize+1, MPI_CHAR,  i,READS_TAG, MPI_COMM_WORLD);
				printf("Rank 0 has sent %lu byte of fastQ to rank %d\n",Readsbuffersize, i);
				allReads[i].erase();
				delete sendReadsBuffer;
			}
	}else{
		//recieving reads
		MPI_Recv(&Readsbuffersize,1, MPI_UNSIGNED_LONG, 0, READS_SIZE_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		printf("Rank %d size of fastQ %lu\n", rank,Readsbuffersize);
		char* myReadsBuffer= new char[Readsbuffersize+1];
		MPI_Recv(myReadsBuffer,Readsbuffersize+1, MPI_CHAR, 0, READS_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		myReads=myReadsBuffer;
		delete myReadsBuffer;
	}
	vector<string> vectorReads=parseNGS(myReads, 0, rank);
	unsigned long int* myCounts = new unsigned long int[allGenes.size()];
	memset(myCounts, 0, allGenes.size()*sizeof(*myCounts));
	hitCounter(allGenes, vectorReads, myCounts, rank);
	//sending back reads number and counts
	if (rank==0){
		basic_string<char>::size_type totalReads=vectorReads.size();
		basic_string<char>::size_type recieverReads;
		unsigned long int* reciverCounts = new unsigned long int[allGenes.size()];
		for (int i=1; i<size;i++){ 
			MPI_Recv(&recieverReads,1, MPI_UNSIGNED_LONG, i, READS_NUMBER_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			totalReads+=recieverReads;
			MPI_Recv(reciverCounts,allGenes.size(), MPI_UNSIGNED_LONG, i, COUNTS_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			for (long unsigned int gn=0; gn<allGenes.size(); gn++) myCounts[gn]+=reciverCounts[gn];
		}
		printf("Total reads: %lu\n",totalReads);
		printf("Gene name : reads : lenght : reads per kilobase per milloin of TOTAL reads\n");
		float normalizedReads;
		for (long unsigned int gn=0; gn<allGenes.size(); gn++){ 
			normalizedReads=float(myCounts[gn])*1000000000/(allGenes[gn].seq.size()*totalReads);
			printf("%s\t %lu \t %lu \t %f\n",allGenes[gn].name.c_str(),myCounts[gn],allGenes[gn].seq.size(), normalizedReads);
		}
	}else{
		basic_string<char>::size_type myReadsNumb=vectorReads.size();
		MPI_Send(&myReadsNumb,1, MPI_UNSIGNED_LONG,  0,READS_NUMBER_TAG, MPI_COMM_WORLD);
		MPI_Send(myCounts,allGenes.size(), MPI_UNSIGNED_LONG,  0,COUNTS_TAG, MPI_COMM_WORLD);
	}
	MPI_Finalize();
	return 0;
}
