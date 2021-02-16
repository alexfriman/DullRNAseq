#include <vector> 
#include <string>
using namespace std;
std::string readGenome(char* fname);
struct gene {
    std::string name;
    std::string seq;
};
std::vector<gene> parseGenome(char* genome, int rank);
