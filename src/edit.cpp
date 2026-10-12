#include "edit.h"

//' reverse
 //'
 //' This function reverses a string.
 //'
 //' @param s The string to be inversed.
std::vector<std::string> reverse(std::vector<std::string> s) {
    int n = s.size();
    for (int i = 0; i < n; i++) {
        std::reverse(s[i].begin(), s[i].end());
    }
    return(s);
}

//' kmer
 //'
 //' This function kmerizes a string
 //'
 //' @param s The string to be inversed.
 //' @param k The size of the kmer.
 //' @param step The distance between two kmers.
 //' @return A string set including all unique kmers.
set<std::string> kmer(std::vector<std::string> s, int k,int step) {
    int n = s.size();
    set<std::string> out;

    for (int i = 0; i < n; i++) {
        int l = s[i].length();
        for (int j = 0; j < l - k + 1; j+=step) {
            out.insert(s[i].substr(j, k));
        }
    }

    return(out);
}

//' kmer
 //'
 //' This function kmerizes a string
 //'
 //' @param s The string to be inversed.
 //' @param k The size of the kmer.
 //' @param step The distance between two kmers.
 //' @return A string vector including all kmers.
std::vector<std::string> kmer(std::string s, int k,int step) {
    int s_len = s.length();
    std::vector<std::string> out;

    for (int i = 0; i < s_len - k + 1; i+=step) {
        out.push_back(s.substr(i, k));
    }
    return(out);
}

//' editDist
 //'
 //' This function calculates the edit distance between two strings.
 //'
 //' @param word1 The input string 1.
 //' @param word2 The input string 2.
 //' @return An int as the edit distance.
int editDist(std::string word1, std::string word2) {

    std::vector<std::vector<int>> dp = std::vector<std::vector<int>>(word1.size() + 1, std::vector<int>(word2.size() + 1, 0));

    for (int i = 0; i <= word1.size(); i++) {
        dp[i][0] = i;
    }

    for (int j = 1; j <= word2.size(); j++) {
        dp[0][j] = j;
    }

    for (int i = 0; i < word1.size(); i++) {
        for (int j = 0; j < word2.size(); j++) {
            if (word1[i] == word2[j]) {
                dp[i + 1][j + 1] = dp[i][j];
            }
            else {
                dp[i + 1][j + 1] = std::min(std::min(dp[i][j + 1], dp[i + 1][j]), dp[i][j]) + 1;
            }
        }
    }

    return dp[word1.size()][word2.size()];
}

//' traceback
 //'
 //' This function to find the position at which two strings have the minimal edit distance.
 //'
 //' @param dp The matrix of atring alignment.
 //' @return A pair, the fist element is the minimal edit distance, the second element
 //' is the position.
pair<int, int> traceback(std::vector<std::vector<int>> dp) {
    int nrow = dp.size();
    int ncol = dp[0].size();

    int minEdit = nrow;
    int pos = 0;

    for (int i = 0; i < ncol; i++) {
        if (dp[nrow - 1][i] < minEdit) {
            minEdit = dp[nrow - 1][i];
            pos = i;
        }
    }
    //return(pair<int, int>(minEdit, pos));

    int loc = nrow - 1;
    while (loc > 1) {
      int temp = std::min(std::min(dp[loc][pos-1], dp[loc-1][pos]), dp[loc-1][pos-1]);
      if (dp[loc - 1][pos - 1] == temp) {
        loc = loc - 1;
        pos = pos - 1;
      }
      else if (dp[loc - 1][pos] == temp) {
        loc = loc - 1;
      }
      else {
        pos = pos - 1;
      }
    }
    return(pair<int, int>(minEdit, pos-1));
}

//' minEditDist
 //'
 //' This function to find best alignemnt between the sequence and barcode with the minimal edit distance.
 //'
 //' @param seq The input tag sequence.
 //' @param barcode The barcode to be matched.
 //' @return A pair, the fist element is the minimal edit distance, the second element
 //' is the position.
pair<int, int> minEditDist(std::string seq, std::string barcode) {
    int bs = barcode.size();
    int ss = seq.size();

    std::vector<std::vector<int>> dp = std::vector<std::vector<int>>(bs + 1, std::vector<int>(ss + 1, 0));

    for (int i = 0; i <= bs; i++) {
        dp[i][0] = i;
    }

    for (int i = 0; i < bs; i++) {
        for (int j = 0; j < ss; j++) {
            if (barcode[i] == seq[j]) {
                dp[i + 1][j + 1] = dp[i][j];
            }
            else {
                dp[i + 1][j + 1] = std::min(std::min(dp[i][j + 1], dp[i + 1][j]), dp[i][j]) + 1;
            }
        }
    }
    /*
    for (int i = 0; i < bs; i++) {
        for (int j = 0; j < ss; j++) {
            cout << dp[i][j] << " ";
        }
        cout << endl;
    }
    */
    return traceback(dp);
}
