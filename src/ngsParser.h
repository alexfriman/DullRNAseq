#include <vector> 
void readNSG(char* fname, std::string* portions, int clusterSize);
#define MIN_SEQ 10
std::vector<std::string> parseNGS(std::string rawreads, float acceptedError, int rank);