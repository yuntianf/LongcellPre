#include "umi_dist.h"
#include "reads_filter.h"



int iso2_dis(std::string a,std::string b,
             const std::string split,
             const std::string sep){
  if(a == b){
    return(0);
  }
  std::vector<int> sites_a = isoform2sites(a, split,sep);
  std::vector<int> sites_b = isoform2sites(b, split,sep);

  int a_size = sites_a.size();
  int b_size = sites_b.size();

  int i =0,j = 0, intersect = 0;
  while(i < a_size/2){
    if(sites_a[2*i] > sites_b[2*j+1]){
      if(j == b_size/2-1){
        break;
      }
      j++;
    }
    else if(sites_a[2*i+1] < sites_b[2*j]){
      i++;
    }
    else{
      intersect += bin2_intersect(sites_a[2*i],sites_a[2*i+1],
                                  sites_b[2*j],sites_b[2*j+1]);
      i++;
    }
  }

  int diff = iso_len(sites_a) + iso_len(sites_b) - 2*intersect;
  return(diff);
}
// [[Rcpp::export]]
DataFrame isos_dis(const std::vector<std::string> isoforms,const int thresh,
             const std::string split,const std::string sep){
  int n = isoforms.size();

  std::vector<int> node1;
  std::vector<int> node2;
  std::vector<int> dis;
  for(int i =0;i < n;i++){
    for(int j = i;j < n;j++){
      int temp = iso2_dis(isoforms[i],isoforms[j],split,sep);
      if(temp <= thresh){
        node1.push_back(i+1);
        node2.push_back(j+1);
        dis.push_back(temp);
      }
    }
  }

  DataFrame out = DataFrame::create(Named("node1") = node1 , Named("node2") = node2, Named("dis") = dis);
  return(out);
}

// [[Rcpp::export]]
NumericVector size_filter_cpp(NumericVector size,double ratio){
  int n = size.size();
  double left = 0, right = sum(size);
  int id = 0;
  double diff = 100;
  for(int i = 0;i < n;i++){
    left += size[i];
    right -= size[i];
    double thresh = left-right*ratio/(1-ratio);
    //Rcout << left<< "-" << right*ratio/(1-ratio) << endl;
    if(diff > abs(thresh)){
      diff = abs(thresh);
      id = i;
    }
    if(thresh >= 0){
      break;
    }
  }
  NumericVector weight(n,0.0);
  weight[size > size[id]] = 1;
  weight[size == size[id]] = (sum(size <= size[id])-double(id)-1.0)/sum(size == size[id]);
  return(weight);
}

// Parse a single isoform string into numeric vector of sites
static inline void parse_sites(const char* p, std::vector<double>& v) {
  v.clear();
  double cur = 0.0;
  bool in_num = false;
  bool seen_dot = false;
  double frac = 0.1;

  for (; *p; ++p) {
    char c = *p;
    if (c == 32) continue; // skip spaces

    if ((c >= 48 && c <= 57)) { // digit
      if (!in_num) { in_num = true; cur = 0.0; seen_dot = false; frac = 0.1; }
      if (!seen_dot) {
        cur = cur * 10.0 + (c - 48);
      } else {
        cur += (c - 48) * frac;
        frac *= 0.1;
      }
    } else if (c == 46) { // dot
      if (!in_num) { in_num = true; cur = 0.0; }
      if (!seen_dot) { seen_dot = true; frac = 0.1; }
    } else {
      if (in_num) { v.push_back(cur); in_num = false; }
      // any non-numeric char (e.g., sep or split) is a delimiter
    }
  }
  if (in_num) v.push_back(cur);
}

// [[Rcpp::export]]
Rcpp::NumericVector isos_len_cpp(Rcpp::CharacterVector isos) {
  const int n = isos.size();
  Rcpp::NumericVector out(n);
  std::vector<double> buf; buf.reserve(16);

  for (int i = 0; i < n; ++i) {
    if (isos[i] == NA_STRING) { out[i] = NA_REAL; continue; }
    const char* p = Rf_translateCharUTF8(isos[i]);
    parse_sites(p, buf);
    const int m = (int)buf.size();
    double s = 0.0;
    for (int j = 0; j + 1 < m; j += 2) s += (buf[j + 1] - buf[j]); // by pairs
    out[i] = s + m / 2.0;  // + n/2 like your iso_len()
  }
  return out;
}
