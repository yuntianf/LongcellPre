#include "seqio.h"
#include <Rcpp.h>
#include <iostream>
#include <string>
#include <regex>
#include <algorithm>
#include <vector>

using namespace Rcpp;
using namespace std;
using namespace klibpp;

int baseCount(std::string seq, char base);
std::string reverseComplement(const std::string& sequence);


bool polyADetect(const std::string& seq, int bin = 20, int count = 15, char base = 'A');
size_t polyARm(std::string seq, const int polyA_len = 10);

/*
std::vector<std::string> extractSoftclip(std::string seq,StringVector mark,NumericVector count,
                                         const std::string strand,int toolkit,
                                         const int search_len = 55,const int polyA_bin = 20,
                                         const int polyA_count = 15,const int polyA_len = 10);
 */

DataFrame extractTagFastq(const char* fastq_path,const char* out_path,
                                         const std::string adapter,const int toolkit,
                                         const int window, const int step,
                                         const int left_flank, const int right_flank, const bool drop,
                                         const int polyA_bin, const int polyA_base_count,const int polyA_len);

std::vector<std::string> extractTag(KSeq record, const std::string adapter,const int toolkit,
                                      const int window, const int step,
                                      const int left_flank, const int right_flank, const bool drop,
                                      const int polyA_bin, const int polyA_base_count,const int polyA_len);

int strSlideSearch(std::string seq,const std::string adapter,
                   const int window = 12, const int step = 3,const bool first = true);

std::vector<std::string> strSubset(std::string str,const int window, const int step);

std::string replicate(std::string mode,int times);
