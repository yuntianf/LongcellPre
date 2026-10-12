#include <Rcpp.h>
#include <vector>
#include <string>
#include <utility>
#include <set>
#include <algorithm>
using namespace Rcpp;
using namespace std;

pair<int,int> exonstr2bin(string exon,const string& delimiters);
vector<pair<int,int> >isostr2bins(string isoform,const string& delimiters);
int bin_sum(vector<pair<int,int> > bins);
string paste(vector<string> s,string sep = "|");
vector<string> exon_status(vector<string> exons,string split = "|");

string bin2exonid(const string bin_str,int status,
                  vector<int> start, vector<int> end,
                  vector<string> exon_id,
                  int mid_bias = 0,int end_bias = 10,int end_overlap = 10,
                  string nonsense_label = "N",string split = "|",string sep = ",");
map<string, string> bins2exonids(vector<string> bins,vector<int> status,
                                 vector<int> start, vector<int> end,
                                 vector<string> exon_id,
                                 int mid_bias = 0,int end_bias = 10,
                                 int end_overlap = 10,
                                 string nonsense_label = "N",
                                 string split = "|",string sep = ",");
List isos2exonids_index(vector<string> isoform,
                        vector<int> start, vector<int> end,
                        vector<string> exon_id,
                        int mid_bias = 0,int end_bias = 10,
                        int end_overlap =10,
                        string nonsense_label = "N",
                        string split = "|",string sep = ",");
