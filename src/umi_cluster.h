#include <Rcpp.h>
#include <string>
#include <algorithm>
#include <vector>
#include <map>
using namespace Rcpp;
using namespace std;

int selectSum(std::map<std::string,int> count,std::vector<std::string> id);
DataFrame shareNeighbor(std::vector<std::string> index, List neighbor,NumericVector count);
