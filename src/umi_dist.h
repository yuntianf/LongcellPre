#include <Rcpp.h>
#include<string>
#include <algorithm>
#include<vector>
#include<set>
#include<math.h>
#include<cmath>
#include <stack>

using namespace std;
using namespace Rcpp;

std::vector<std::string> flatten(std::vector<std::vector<std::string> > const &vec);
std::vector<std::string> str_split(std::string s, std::string split);
std::vector<int> isoform2sites(std::string iso,
                                      const std::string split = "|",const std::string sep = ",");

int iso_len(std::vector<int> sites);
int bin2_intersect(int a_start,int a_end,
                          int b_start,int b_end);
std::vector<int> sites_chop(std::vector<int> sites, int start,int end);
int iso2_mid_dist(std::string a,std::string b,
                    const std::string split = "|",
                    const std::string sep = ",");
NumericMatrix iso_mid_dist(std::vector<std::string> iso_set,int thresh = 80,
                       std::string split = "|",std::string sep = ",");

NumericVector iso2_mid_diff(std::string a,std::string b,
                            const int end_bias,
                  const std::string split,
                  const std::string sep);
NumericMatrix isoset_mid_diff(std::vector<std::string> iso_set1,
                              std::vector<std::string> iso_set2,
                              const int thresh = 3,
                              const double overlap_thresh = 0.5,
                              const int end_bias = 200,
                              std::string split = "|",
                              std::string sep = ",");

int needle(std::string A,std::string B,
                  int match_score = 1,
                  int mismatch_score = -1,
                  int gap_score = -1);

int minEditDis(std::string seq1, std::string seq2, int k = 10);

std::vector<std::vector<double> > umi_graph_table(std::vector<std::string> umi,
                                               std::vector<std::string> isoform,
                                               std::vector<int> count,
                                               int sim_thresh = 5,int iso_thresh = 80,
                                               std::string split = "|",std::string sep = ",");
