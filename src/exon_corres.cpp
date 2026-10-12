#include "splice_site_correct.h"
#include "exon_corres.h"

pair<int,int> exonstr2bin(string exon,const string& delimiters){
  std::vector<string> sites = splice_site_cpp(exon,delimiters);

  pair<int,int> bin;
  bin.first = stoi(sites[0]);
  bin.second = stoi(sites[1]);
  return(bin);
}

std::vector<pair<int,int> >isostr2bins(string isoform,const string& delimiters){
  std::vector<string> sites = splice_site_cpp(isoform,delimiters);

  int isolen = sites.size();
  pair<int,int> bin;
  std::vector<pair<int,int> > out;
  for(int i = 0; i < isolen/2;i++){
    bin.first = stoi(sites[i*2]);
    bin.second = stoi(sites[i*2+1]);
    out.push_back(bin);
  }
  return(out);
}

int bin_sum(std::vector<pair<int,int> > bins){
  int sum = 0;
  for(auto i:bins){
    if(i.first > i.second){
      stop("The end position of each bin should be larger than its start position!");
    }
    sum += i.second-i.first+1;
  }
  return(sum);
}

string paste(std::vector<string> s,string sep){
  string out = "";
  int len = s.size();
  for(int i =0;i < len;i++){
    out+=s[i];
    if(i!=len-1){
      out+=sep;
    }
  }
  return(out);
}

std::vector<string> exon_status(std::vector<string> exons,string split){
  int len = exons.size();
  if(len == 1){
    exons[0] = exons[0]+split+"2";
  }
  for(int i = 0;i < len;i++){
    if(i == 0){
      exons[i] = exons[i]+split+"-1";
    }
    else if(i == len-1){
      exons[i] = exons[i]+split+"1";
    }
    else{
      exons[i] = exons[i]+split+"0";
    }
  }
  return(exons);
}

string bin2exonid(const string bin_str,int status,
                std::vector<int> start, std::vector<int> end,
                std::vector<string> exon_id,
                int mid_bias,int end_bias, int end_overlap,
                string nonsense_label,string split,string sep){
  if(status < -1 || status > 2){
    stop("There are only four kinds of status: -1=start exon; 0=middle exon; 1=end exon; 2=only one exon. Please specify status within them four");
  }
  pair<int,int> bin = exonstr2bin(bin_str,sep);

  int exon_size = exon_id.size();

  std::vector<string> exon_seq;

  std::vector<pair<int,int> > left;

  bool flag = false;
  for(int i = 0;i < exon_size;i++){
    if(start[i] <= bin.second && end[i] >= bin.first){
      flag = true;
      if(status == 0){
        exon_seq.push_back(exon_id[i]);
      }
      else{
        if(min(bin.second,end[i])-max(bin.first,start[i])+1 >= end_overlap){
          exon_seq.push_back(exon_id[i]);
        }
      }
      if(start[i] != bin.first){
        if(status == 0 || status == 1){
          left.push_back({min(start[i],bin.first),max(start[i],bin.first)-1});
        }
        else if(status == -1 || status == 2){
          if(bin.first < start[i]-end_bias){
            left.push_back({bin.first,start[i]-1});
          }
        }
      }

      bin.first = end[i]+1;
      if(end[i] > bin.second){
        if(status == 0 || status == -1){
          left.push_back({bin.second+1,end[i]});
        }
        break;
      }
    }
    else{
      if(flag){
        break;
      }
    }
  }

  int left_len =bin_sum(left);
  int end_len = max(bin.second-bin.first+1,0);
  if(status == 1 || status == 2){
    if(end_len > end_bias){
      left_len += end_len;
    }
  }
  else{
    left_len += end_len;
  }

  if(left_len > mid_bias){
    return(nonsense_label);
  }

  return(paste(exon_seq,split));
}

map<string, string> bins2exonids(std::vector<string> bins,std::vector<int> status,
                             std::vector<int> start, std::vector<int> end,
                             std::vector<string> exon_id,
                             int mid_bias,int end_bias,
                             int end_overlap,
                             string nonsense_label,
                             string split,string sep){
  int bins_size = bins.size();
  if(bins.size()!= status.size()){
    stop("Each bin should have its corresponding status indication!");
  }

  map<string, string> exonid_vec;
  string temp;
  for(int i = 0;i < bins_size;i++){
    temp = bin2exonid(bins[i],status[i],start,end,exon_id,
                      mid_bias,end_bias,end_overlap,nonsense_label,split);
    exonid_vec[bins[i]+to_string(status[i])] = temp;
  }
  return(exonid_vec);
}



// [[Rcpp::export]]
List isos2exonids_index(std::vector<std::string> isoform,
                                 std::vector<int> start, std::vector<int> end,
                                 std::vector<std::string> exon_id,
                                 int mid_bias,int end_bias,
                                 int end_overlap,
                                 std::string nonsense_label,
                                 std::string split,std::string sep){
  std::set<std::string> exonstr_uniq_status_set;
  std::map<std::string, std::vector<std::string> > iso_exonstr;

  for(auto iso:isoform){
    std::vector<std::string> exonstr = splice_site_cpp(iso,split);

    std::vector<std::string> exonstr_status = exon_status(exonstr);
    iso_exonstr[iso] = exonstr_status;

    for(auto temp:exonstr_status){
      exonstr_uniq_status_set.insert(temp);
    }
  }

  std::vector<std::string> exonstr_uniq_status_vec(exonstr_uniq_status_set.begin(),
                                         exonstr_uniq_status_set.end());

  std::map<std::string,std::string> exonstr_status2exon_id;
  for(auto exonstr:exonstr_uniq_status_vec){
    std::vector<std::string> temp = splice_site_cpp(exonstr,split);
    exonstr_status2exon_id[exonstr] = bin2exonid(temp[0],stoi(temp[1]),start,end,
                                                 exon_id,mid_bias,end_bias,end_overlap,
                                                 nonsense_label,split,sep);
  }

  List iso_exonid;
  for(auto iso:isoform){
    std::string exonid = "";
    for(auto exonstr:iso_exonstr[iso]){
      exonid += exonstr_status2exon_id[exonstr]+split;
    }
    exonid = exonid.substr(0,exonid.size()-1);
    iso_exonid.push_back(exonid,iso);
  }

  return(iso_exonid);
}
