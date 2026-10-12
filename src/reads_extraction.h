#include <Rcpp.h>
#include <iostream>
#include <string>
#include <regex>
#include <algorithm>

using namespace Rcpp;
using namespace std;


bool within(int posA,int posB,const int bin = 2);
int seqEnd(int start,StringVector mark,NumericVector count);

List cigarProcess(std::string cigar);
bool cigarcheck(StringVector mark,NumericVector count);

int isoformEnd(NumericMatrix isoform, std::string strand);
std::string isoform2string(NumericMatrix isoform,std::string sep = "|");
NumericMatrix extractIsoform(int read_start,
                             StringVector mark, NumericVector count,
                             int refer_start, int refer_end,
                             int flank = 200);
NumericMatrix chopOutBound(NumericMatrix isoform,NumericMatrix annotation,string strand);
NumericMatrix splicesiteCorrect(NumericMatrix isoform,NumericMatrix annotation,
                                string strand,int bin = 2);

DataFrame extractReads(std::vector<std::string> seq,std::vector<std::string> cigar,
                       NumericVector pos,
                       NumericMatrix annotation,const std::string strand,
                       int toolkit,
                       int end_flank = 200,
                       int splice_site_bin = 2);
