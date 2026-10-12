#pragma once
#include <vector>
#include <string>
#include <set>
#include <algorithm>
#include<math.h>
#include<cmath>
#include <stdlib.h>
#include <utility>
#include <iostream>
#include <fstream>
#include <Rcpp.h>
using namespace std;
using namespace Rcpp;

double mean(std::vector<double> num);
double var(std::vector<double> num);

double update_sigma(std::vector<double> num, double sigma_start);
double update_prob(std::vector<double> m, double n);

std::vector<int> kmer_include(std::string seq, set<std::string> dic);

std::vector<std::vector<int> > barcodes_cos_vec(std::vector<std::string> barcodes, set<std::string> dic);
double cos_sim(std::vector<int> a, std::vector<int> b);
bool cos_sim_comp(std::pair <int, double> cos1, std::pair <int, double> cos2);

std::vector<std::string> barcode_cand_cos(std::string seq, std::vector<std::string> barcodes, set<std::string> dic, std::vector<std::vector<int> > index,
    int top = 8, double thresh = 0.25);
std::vector<int> pos_filter(std::vector<double> start, std::vector<int> edit);
DataFrame barcodeMatch(std::vector<std::string> seq, std::vector<std::string> barcodes,
                       double mu=20, double sigma=10, double sigma_start=10, int k=8, int batch=100,
                       int top=8, double cos_thresh=0.25, double alpha=0.05, int edit_thresh=5,
                       const int UMI_len = 10, const int flank = 1);
std::vector<std::vector<std::string> > NeighborExtract(std::vector<std::string> reads,std::vector<int> id,
                                                       std::vector<int> pos,
                                                       const int UMI_len, const int flank,const int bar_len);
