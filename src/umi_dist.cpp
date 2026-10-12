#include "umi_dist.h"
#include "edit.h"

#include <unordered_map>
#include <array>
#include <cstdint>
#include <limits>

// ##### isoform cluster #####
std::vector<std::string> flatten(std::vector<std::vector<std::string> > const &vec){
  std::vector<std::string> flattened;
  for (auto const &v: vec) {
    flattened.insert(flattened.end(), v.begin(), v.end());
  }
  return flattened;
}


std::vector<std::string> str_split(std::string s, std::string split){
  std::vector<std::string> out;

  size_t pos = 0;
  std::string token;
  while ((pos = s.find(split)) != std::string::npos) {
    token = s.substr(0, pos);
    out.push_back(token);
    s.erase(0, pos + split.length());
  }
  if(s.size() > 0){
    out.push_back(s);
  }

  return(out);
}

// [[Rcpp::export]]
std::vector<int> isoform2sites(std::string iso,
                                      const std::string split,
                                      const std::string sep){
  std::vector<std::string> exons = str_split(iso,split);

  std::vector<std::vector<std::string> > bins;
  int len = exons.size();
  for(int i = 0;i < len;i++){
    bins.push_back(str_split(exons[i],sep));
  }
  std::vector<std::string> sites_chr = flatten(bins);

  std::vector<int> sites_int;
  transform(sites_chr.begin(), sites_chr.end(), back_inserter(sites_int),
            [](const std::string& str) { return static_cast<int>(std::stod(str)); });

  return(sites_int);
}


int iso_len(std::vector<int> sites){
  int sites_size = sites.size();

  if(sites_size % 2 > 0){
    cout << "the length of sites should be an even" << endl;
    stop("Isoform coordinates must contain complete start/end pairs.");
  }
  int len = 0;
  for(int i = 0;i < sites_size/2;i++){
    int sites_s = sites[i*2];
    int sites_e = sites[i*2 + 1];

    len= len+sites_e-sites_s+1;
  }
  return(len);
}


int bin2_intersect(int a_start,int a_end,
                          int b_start,int b_end){
  int len = 0;
  if(a_start <= b_end && a_end >= b_start){
    len = min({a_end,b_end})-max({a_start,b_start}) + 1;
  }
  return(len);
}

// [[Rcpp::export]]
std::vector<int> sites_chop(std::vector<int> sites, int start,int end){
  int sites_size = sites.size();

  if(sites_size % 2 > 0){
    cout << "the length of sites should be an even" << endl;
    stop("Isoform coordinates must contain complete start/end pairs.");
  }

  std::vector<int> chop;
  for(int i = 0;i < sites_size/2;i++){
    int sites_s = sites[i*2];
    int sites_e = sites[i*2 + 1];
    if(sites_e <= start || sites_s >= end){
      continue;
    }
    else{
      int temp = max({sites_s,start});
      chop.push_back(max({sites_s,start}));
      chop.push_back(min({sites_e,end}));
    }
  }
  return(chop);
}

// [[Rcpp::export]]
int iso2_mid_dist(std::string a,std::string b,
                    const std::string split,
                    const std::string sep){
  if(a == b){
    return(0);
  }
  std::vector<int> sites_a = isoform2sites(a, split,sep);
  std::vector<int> sites_b = isoform2sites(b, split,sep);

  int a_size = sites_a.size();
  int b_size = sites_b.size();

  int start = max({sites_a[0],sites_b[0]});
  int end = min({sites_a[a_size - 1],sites_b[b_size - 1]});

  if(start >= end){
    return(-1);
  }

  std::vector<int> chop_a = sites_chop(sites_a, start, end);
  std::vector<int> chop_b = sites_chop(sites_b, start, end);

  int chop_a_size = chop_a.size();
  int chop_b_size = chop_b.size();

  int intersect = 0;
  int next_start = 0;
  for(int i =0;i < chop_a_size/2;i++){
    int j = next_start;
    bool flag = false;

    while(j < chop_b_size/2){
      int a_start = chop_a[2*i];
      int a_end = chop_a[2*i + 1];
      int b_start = chop_b[2*j];
      int b_end = chop_b[2*j + 1];

      int inter = bin2_intersect(a_start,a_end,b_start,b_end);

      if(inter > 0){
        if(!flag){
          next_start = j;
          flag = true;
        }
        intersect += inter;
      }
      else{
        if(flag){
          break;
        }
      }
      j++;
    }
  }

  int dis = iso_len(chop_a) + iso_len(chop_b) - 2*intersect;
  return(dis);
}

NumericMatrix iso_mid_dist(std::vector<std::string> iso_set,
                           int thresh,
                           std::string split,std::string sep){
  int iso_size = iso_set.size();

  NumericMatrix iso_dis ((iso_size+1)*iso_size/2,2);
  int id = 0;

  for(int i = 0;i < iso_size;i++){
    for(int j = i;j < iso_size;j++){
      if(i == j){
        iso_dis(id,0) = i;
        iso_dis(id,1) = j;
        id = id + 1;
      }
      else{
        int dis = iso2_mid_dist(iso_set[i],iso_set[j],
                           split,sep);
        //cout << i << "-" << j<<":" << dis << endl;
        if(dis >=0 && dis <= thresh){
          iso_dis(id,0) = i;
          iso_dis(id,1) = j;
          id = id+ 1;
        }
      }
    }
  }
  return(iso_dis(Range(0,id-1),_));
}

// [[Rcpp::export]]
NumericVector iso2_mid_diff(std::string a,std::string b,
                            const int end_bias,
                 const std::string split,
                 const std::string sep){
  NumericVector out(2);
  if(a == b){
    out = {0,1};
    return(out);
  }
  std::vector<int> sites_a = isoform2sites(a, split,sep);
  std::vector<int> sites_b = isoform2sites(b, split,sep);

  int a_size = sites_a.size();
  int start = sites_a[0];
  int end = sites_a[a_size - 1];

  std::vector<int> chop_b = sites_chop(sites_b, start, end);

  int b_size = chop_b.size();

  if(b_size == 0){
    out = {-1,-1};
    return(out);
  }

  if(start <= chop_b[0] && start >= chop_b[0]-end_bias &&
    chop_b[0] <= sites_a[1]){
    start = chop_b[0];
  }
  if(end >= chop_b[b_size - 1] && end <= chop_b[b_size - 1]+end_bias &&
     chop_b[b_size - 1] >= sites_a[a_size - 2]){
    end = chop_b[b_size - 1];
  }

  if(start >= end){
    out = {-1,-1};
    return(out);
  }


  //Rcout << start << " " << end << endl;

  std::vector<int> chop_a = sites_chop(sites_a, start, end);

  int chop_a_size = chop_a.size();
  int chop_b_size = chop_b.size();

  /*
  for (int i = 0; i < chop_a_size; ++i) {
    Rcpp::Rcout << chop_a[i] << " ";
  }
  Rcout<< endl;
  for (int i = 0; i < chop_b_size; ++i) {
    Rcpp::Rcout << chop_b[i] << " ";
  }
  Rcout<< endl;
   */

  int intersect = 0;
  int next_start = 0;
  for(int i =0;i < chop_a_size/2;i++){
    int j = next_start;
    bool flag = false;

    while(j < chop_b_size/2){
      int a_start = chop_a[2*i];
      int a_end = chop_a[2*i + 1];
      int b_start = chop_b[2*j];
      int b_end = chop_b[2*j + 1];

      int inter = bin2_intersect(a_start,a_end,b_start,b_end);

      if(inter > 0){
        if(!flag){
          next_start = j;
          flag = true;
        }
        intersect += inter;
      }
      else{
        if(flag){
          break;
        }
      }
      j++;
    }
  }

  // cout << iso_len(chop_a) << " " << iso_len(chop_b) << " " << intersect << endl;

  double dis = iso_len(chop_a) + iso_len(chop_b) - 2*intersect;
  double ratio = double(intersect)/double(iso_len(sites_b));

  out = {dis,ratio};
  return(out);
}

// [[Rcpp::export]]
NumericMatrix isoset_mid_diff(std::vector<std::string> iso_set1,
                              std::vector<std::string> iso_set2,
                              const int thresh,
                              const double overlap_thresh,
                              const int end_bias,
                              std::string split,
                              std::string sep){
  int iso1_size = iso_set1.size();
  int iso2_size = iso_set2.size();

  std::vector<std::array<double,4> > accepted;
  for(int i = 0; i < iso1_size; ++i){
    if((i & 127) == 0) checkUserInterrupt();
    for(int j = 0; j < iso2_size; ++j){
      NumericVector value = iso2_mid_diff(iso_set1[i],iso_set2[j],end_bias,split,sep);
      const double dis = value[0], ratio = value[1];
      if(dis >= 0 && dis <= thresh && ratio >= overlap_thresh)
        accepted.push_back({double(i),double(j),dis,ratio});
    }
  }
  if(accepted.empty()) {
    NumericMatrix out(1,4);
    out(0,0) = -1;
    return out;
  }
  NumericMatrix out(accepted.size(),4);
  for(size_t i = 0; i < accepted.size(); ++i)
    for(size_t j = 0; j < 4; ++j) out(i,j) = accepted[i][j];
  return out;
}


// ##### UMI distance #####
// ### needleman ###
int needle(std::string A,std::string B,
           int match_score, int mismatch_score, int gap_score){
  const size_t n = A.size(), m = B.size();
  std::vector<int> previous(m + 1), current(m + 1);
  for (size_t j = 0; j <= m; ++j) previous[j] = static_cast<int>(j) * gap_score;
  for (size_t i = 1; i <= n; ++i) {
    current[0] = static_cast<int>(i) * gap_score;
    for (size_t j = 1; j <= m; ++j) {
      const int score = A[i-1] == B[j-1] ? match_score : mismatch_score;
      current[j] = std::max({previous[j-1] + score, previous[j] + gap_score,
                             current[j-1] + gap_score});
    }
    previous.swap(current);
  }
  return previous[m];
}

// ### edit distance ###
int minEditDis(std::string seq1, std::string seq2, int k) {
  std::vector<std::string> seq1_kmer = kmer(seq1, k,1);
  std::vector<std::string> seq2_kmer = kmer(seq2, k,1);

  int kmer_size = seq1_kmer.size();

  int edit = seq1.size();
  for(int i = 0;i < kmer_size;i++){
    for(int j = 0;j < kmer_size;j++){
      int dis = editDist(seq1_kmer[i],seq2_kmer[j]);
      if(edit > dis){
        edit = dis;
      }
      if(dis == 0){
        break;
      }
    }
  }
  return(edit);
}


// Scores are cached only as an optimization. Collisions cause recomputation,
// never a changed score or a dropped edge. Storage is bounded (~1 MiB), local
// to this graph, and released on return. No read/UMI-index limit is imposed.
class UmiScoreCache {
  struct Entry { uint64_t key; int score; };
  std::vector<Entry> slots;
public:
  UmiScoreCache() : slots(65536, {std::numeric_limits<uint64_t>::max(), 0}) {}
  int score(uint32_t a, uint32_t b, const std::vector<std::string>& unique) {
    if (a == b) return unique[a].size();
    if (a > b) std::swap(a, b);
    const uint64_t key = (static_cast<uint64_t>(a) << 32) | b;
    const size_t slot = ((key ^ (key >> 33)) * UINT64_C(0xff51afd7ed558ccd)) & (slots.size()-1);
    Entry& entry = slots[slot];
    if (entry.key != key) { entry.key = key; entry.score = needle(unique[a], unique[b]); }
    return entry.score;
  }
};

// [[Rcpp::export]]
std::vector<std::vector<double> > umi_graph_table(std::vector<std::string> umi,
                              std::vector<std::string> isoform,
                              std::vector<int> count,
                              int sim_thresh,int iso_thresh,
                              std::string split,std::string sep){
  const size_t n = umi.size();
  if (isoform.size() != n || count.size() != n)
    stop("The sizes of isoforms, UMIs and counts must match.");
  std::vector<std::vector<double> > edges;
  if (!n) return edges;
  std::map<std::string, std::vector<size_t> > by_isoform;
  std::unordered_map<std::string, uint32_t> ids;
  std::vector<std::string> unique;
  std::vector<uint32_t> correspondence;
  std::vector<std::array<unsigned int, 256> > composition;
  std::array<bool, 256> alphabet_used = {};
  for (size_t i = 0; i < n; ++i) {
    if (umi[i].empty() || count[i] < 0 || count[i] == NA_INTEGER)
      stop("UMIs must be non-empty and counts must be non-negative integers.");
    by_isoform[isoform[i]].push_back(i);
    auto found = ids.find(umi[i]);
    if (found == ids.end()) {
      const uint32_t id = unique.size();
      ids.emplace(umi[i], id);
      unique.push_back(umi[i]);
      correspondence.push_back(id);
      std::array<unsigned int, 256> frequencies = {};
      for (unsigned char base : umi[i]) { ++frequencies[base]; alphabet_used[base] = true; }
      composition.push_back(frequencies);
    } else correspondence.push_back(found->second);
  }
  UmiScoreCache cache;
  std::vector<size_t> alphabet;
  for (size_t base = 0; base < 256; ++base) if (alphabet_used[base]) alphabet.push_back(base);
  Function getOption = Environment::base_env()["getOption"];
  const double max_edges = as<double>(getOption("longcell.max_graph_edges", 5000000.0));
  if (!R_finite(max_edges) || max_edges < 1)
    stop("longcell.max_graph_edges must be a finite positive number.");
  uint64_t comparisons = 0;
  // map iteration retains the original lexicographic isoform-pair order.
  // Stream compatible pairs instead of allocating a square isoform matrix.
  for (auto first = by_isoform.begin(); first != by_isoform.end(); ++first) {
    checkUserInterrupt();
    for (auto second = first; second != by_isoform.end(); ++second) {
      if (first != second) {
        const int distance = iso2_mid_dist(first->first, second->first, split, sep);
        if (distance < 0 || distance > iso_thresh) continue;
      }
      for (size_t a : first->second) for (size_t b : second->second) {
        if ((++comparisons & 16383) == 0) checkUserInterrupt();
        const uint32_t ia = correspondence[a], ib = correspondence[b];
        if (ia != ib) {
          unsigned int common = 0;
          for (size_t base : alphabet)
            common += std::min(composition[ia][base], composition[ib][base]);
          // Admissible upper bound for match=1, mismatch=gap=-1.
          // Reject only pairs whose exact score cannot pass the threshold.
          if (2 * static_cast<int>(common) - static_cast<int>(std::max(umi[a].size(), umi[b].size())) < sim_thresh)
            continue;
        }
        const int score = cache.score(ia, ib, unique);
        if (score >= sim_thresh) {
          if (edges.size() >= max_edges)
            stop("UMI graph exceeds longcell.max_graph_edges. No group was silently dropped; increase the limit only with sufficient memory.");
          const double weight = a == b ? static_cast<double>(count[a]) * (count[a]-1) / 2 :
                                        static_cast<double>(count[a]) * count[b];
          edges.push_back({static_cast<double>(a), static_cast<double>(b),
                           static_cast<double>(score), weight});
        }
      }
    }
  }
  return edges;
}
