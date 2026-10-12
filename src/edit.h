#pragma once
#include<string>
#include <algorithm>
#include<vector>
#include<set>
#include<math.h>
#include<cmath>
#include<fstream>
#include <stdlib.h>
#include <utility>
using namespace std;


std::vector<std::string> reverse(std::vector<std::string> s);
set<std::string> kmer(std::vector<std::string> s, int k,int step);
std::vector<std::string> kmer(std::string s, int k, int step);
int editDist(std::string word1, std::string word2);
std::pair<int,int> traceback(std::vector<std::vector<int>> dp);
std::pair<int, int> minEditDist(std::string seq, std::string barcode);
