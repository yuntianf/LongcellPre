#include <Rcpp.h>
#include <string>
#include <algorithm>
#include <vector>
#include <unordered_map>
using namespace Rcpp;
using namespace std;

string getSubString(const string& strValue,const string& startChar,const string& endChar);
bool isin(const string element,const StringVector vec);

map<std::string, int> isoform_count(vector<string> isoform,const string sep);
vector<string> splice_site_cpp(const string& input, const string& delimiters);
unordered_map<string, int> splice_site_count_cpp(vector<string> isoform,
                                                   const string split,const string sep);


List splice_site_table_cpp(vector<string> isoform,
                           const string split = "|",const string sep = ",",
                           const int splice_site_thresh = 10);

LogicalMatrix matrix_xor(IntegerMatrix mat);
