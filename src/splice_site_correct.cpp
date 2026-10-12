#include "splice_site_correct.h"

std::vector<std::string> splice_site_cpp(const std::string& input, const std::string& delimiters) {
  std::vector<std::string> tokens;
  istringstream stream(input);
  std::string token;

  while (getline(stream, token)) {
    size_t startPos = 0;
    size_t endPos;

    while ((endPos = token.find_first_of(delimiters, startPos)) != std::string::npos) {
      if (endPos != startPos) {
        tokens.push_back(token.substr(startPos, endPos - startPos));
      }

      startPos = endPos + 1;
    }

    if (startPos < token.length()) {
      tokens.push_back(token.substr(startPos));
    }
  }

  return tokens;
}

std::string getSubString(
    const std::string& strValue,
    const std::string& startChar,
    const std::string& endChar)
{
  std::string subString = "";
  // Get index position of first character
  size_t startPos = strValue.find(startChar);
  // Get index position of second character
  size_t endPos = strValue.rfind(endChar);
  // Check if both the characters exists in string
  if( startPos != std::string::npos &&
      endPos != std::string::npos && startPos != endPos)
  {
    // Get substring from start to end character
    subString = strValue.substr(startPos + 1, endPos - startPos - 1);
  }
  return subString;
}

std::map<std::string, int> isoform_count(std::vector<std::string> isoform,const std::string sep){
  std::map<std::string, int> isoCounts;

  int len = isoform.size();
  std::string temp = "";
  for(int i = 0;i < len;i++){
    temp = getSubString(isoform[i],sep,sep);
    isoCounts[temp]++;
  }
  return(isoCounts);
}


std::unordered_map<std::string, int> splice_site_count_cpp(std::vector<std::string> isoform,
                                            const std::string split,const std::string sep){

  std::map<std::string, int> isoCounts = isoform_count(isoform,sep);
  std::unordered_map<std::string, int> siteCounts;

  for (auto pair : isoCounts) {
    std::vector<std::string> sites = splice_site_cpp(pair.first,split+sep);

    int len = sites.size();
    for (int i = 0;i < len;i++) {
      siteCounts[sites[i]] += pair.second;
    }
  }

  return(siteCounts);
}


bool isin(const std::string element,const StringVector vec){
  int len = vec.size();
  for (int i = 0; i < len; i++) {
    if (vec[i] == element)
      return true;
  }
  return false;
}

// [[Rcpp::export]]
List splice_site_table_cpp(std::vector<std::string> isoform,
                                    const std::string split,const std::string sep,
                                    const int splice_site_thresh){

  std::unordered_map<std::string, int> ss_count = splice_site_count_cpp(isoform,split,sep);

  StringVector ss;
  for(auto pair : ss_count){
    if(pair.second >= splice_site_thresh){
      ss.push_back(pair.first);
    }
  }

  int len = ss.size();

  std::vector<int> id_vec;
  std::vector<std::string> start_vec;
  std::vector<std::string> end_vec;

  int iso_size = isoform.size();
  if(len > 0){
    NumericVector temp_ss (len);
    temp_ss.names()=ss;

    NumericMatrix mid_splice_sites(iso_size,len);

    for(int i = 0;i < iso_size;i++){
      std::string iso = isoform[i];
      std::vector<std::string> site = splice_site_cpp(iso,split+sep);
      int iso_len = site.size();
      std::string start = site[0];
      std::string end = site[iso_len-1];

      temp_ss.fill(0);
      bool flag = false;
      for(int i = 1;i < iso_len-1;i++){
        if(isin(site[i],ss)){
          temp_ss[site[i]] = 1;
        }
        else{
          flag = true;
          break;
        }
      }

      if(!flag){
        id_vec.push_back(i+1);
        start_vec.push_back(start);
        end_vec.push_back(end);
        mid_splice_sites(start_vec.size()-1,_) = temp_ss;
      }
    }
    NumericMatrix out_splice_sites;
    if(start_vec.size() > 0){
      out_splice_sites = mid_splice_sites(Range(0,start_vec.size()-1),_);
      colnames(out_splice_sites) = ss;
    }
    List out = List::create(Named("id") = id_vec,
                          Named("start") = start_vec,
                          Named("mid") = out_splice_sites,
                          Named("end") = end_vec);
    return(out);
    }
  else{
    for(int i = 0;i < iso_size;i++){
      std::string iso = isoform[i];
      std::vector<std::string> site = splice_site_cpp(iso,split+sep);
      int iso_len = site.size();
      if(iso_len == 2){
        std::string start = site[0];
        std::string end = site[iso_len-1];
        id_vec.push_back(i+1);
        start_vec.push_back(start);
        end_vec.push_back(end);
      }
    }
    List out = List::create(Named("id") = id_vec,
                          Named("start") = start_vec,
                          Named("end") = end_vec);
    return(out);
  }
}

// [[Rcpp::export]]
LogicalMatrix matrix_xor(IntegerMatrix mat) {
  int n = mat.nrow(), m = mat.ncol();
  LogicalMatrix result(n, n);

  for (int i = 0; i < n; ++i) {
    result(i, i) = true; // self is always non-conflicting
    for (int j = i + 1; j < n; ++j) {
      bool conflict = false;
      for (int k = 0; k < m; ++k) {
        int vi = mat(i, k);
        int vj = mat(j, k);
        if (vi != NA_INTEGER && vj != NA_INTEGER && vi != vj) {
          conflict = true;
          break;
        }
      }
      result(i, j) = result(j, i) = !conflict;
    }
  }

  return result;
}


