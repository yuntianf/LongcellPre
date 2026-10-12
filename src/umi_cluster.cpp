#include "umi_cluster.h"

int selectSum(std::map<std::string,int> count,std::vector<std::string> id){
  int len = id.size(),sum = 0;
  for(int i = 0;i < len;i++){
    sum+=count[id[i]];
  }
  return(sum);
}

// [[Rcpp::export]]
DataFrame shareNeighbor(std::vector<std::string> index, List neighbor,NumericVector count){
  int len = index.size();
  std::map<std::string,int> mapping;
  for(int i = 0;i < len;i++){
    mapping[index[i]] = count[i];
  }

  std::vector<std::string> id1;
  std::vector<std::string> id2;
  std::vector<int> share;

  List neighbor_sort;
  for(int i = 0;i < len;i++){
    std::vector<std::string> node = neighbor[i];
    std::sort(node.begin(), node.end());
    neighbor_sort[index[i]] = node;
  }


  int n = 0;
  std::vector<std::string> intersection;
  for(int i = 0;i < len;i++){
    std::vector<std::string> node1 = neighbor_sort[i];

    for(auto j:node1){
      std::vector<std::string> node2 = neighbor_sort[j];
      intersection.clear();

      std::set_intersection(node1.begin(), node1.end(),
                       node2.begin(), node2.end(),
                       std::back_inserter(intersection));

      n = selectSum(mapping,intersection);

      id1.push_back(index[i]);
      id2.push_back(j);
      share.push_back(n);
    }
  }

  DataFrame df = DataFrame::create(Named("node1") = id1, Named("node2") = id2, Named("share") = share);
  return(df);
}

