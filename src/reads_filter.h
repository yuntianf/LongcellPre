#include <Rcpp.h>
#include <vector>
#include <string>
#include <algorithm>

using namespace Rcpp;
using namespace std;

int iso2_dis(std::string a,std::string b,
                    const std::string split,
                    const std::string sep);

DataFrame isos_dis(const std::vector<std::string> isoforms,const int thresh = 10,
                     const std::string split = "|",const std::string sep = ",");

NumericVector size_filter_cpp(NumericVector size,double ratio = 0.1);

static inline void parse_sites(const char* p, std::vector<double>& v);

Rcpp::NumericVector isos_len_cpp(Rcpp::CharacterVector isos);

