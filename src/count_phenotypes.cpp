// [[Rcpp::plugins("cpp11")]]
#include <Rcpp.h>
#include <iostream>


using namespace Rcpp;
using namespace std;

// Integer vector hashing
struct HASH {
  int operator()(const vector<int> &V) const {
    int hash=0;
    for(int i=0;i<V.size();i++) {
      hash^=V[i];
    }
    return hash;
  }
};

// [[Rcpp::export]]
IntegerMatrix group_and_reduce_phenotypes(IntegerMatrix phenotypes, IntegerVector markers_to_drop) {
  
  int nrows = phenotypes.rows();
  int ncols = phenotypes.cols();
  int n_markers = markers_to_drop.length();

  CharacterVector ch = colnames(phenotypes);  
  
  //Remove Columns
  int cols_removed = 0;
  for(int j = 0; j<n_markers; j++){
    if(markers_to_drop(j)){
      cols_removed++;
      for(int i = 0; i<nrows; i++){
        phenotypes(i,j) = -1;
      }
    }
  }
  
  
  //If no columns were removed, return phenotyps
  if(!cols_removed) return(phenotypes);
  
  
  
  unordered_map<vector<int>,int,HASH> phen_map;
  
  int unique_counter = 0;
  for(int i = 0; i < nrows; i++){
    
    vector<int> phen_vec(n_markers);
    for(int j = 0; j < n_markers;j++){
      phen_vec[j] = phenotypes(i,j);
    }
    
    
    if (phen_map.find(phen_vec) == phen_map.end()){
      phen_map[phen_vec] = unique_counter;
      unique_counter++;
    }
  }
  
  
  // Create new matrix
  IntegerMatrix collapsed_phenotypes(unique_counter, ncols);
  colnames(collapsed_phenotypes) = ch;
  
  // Paste phenotypes
  for (auto& it: phen_map) {
    for(int i = 0; i < n_markers; i++) collapsed_phenotypes(it.second, i) = it.first[i];
  }
  
  // Sum counts
  for(int i = 0; i < nrows; i++){
    
    vector<int> phen_vec(n_markers);
    
    for(int j = 0; j < n_markers;j++){
      phen_vec[j] = phenotypes(i,j);
    }
    
    unordered_map<vector<int>,int,HASH>::const_iterator u_phen = phen_map.find(phen_vec);
    for(int j = n_markers; j < ncols; j++){
      collapsed_phenotypes(u_phen->second,j) = collapsed_phenotypes(u_phen->second,j) + phenotypes(i,j);
    }
    
  }
  
  return(collapsed_phenotypes);
}
