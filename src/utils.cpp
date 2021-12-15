#include <Rcpp.h>
#include<iostream>


using namespace Rcpp;
using namespace std;

//' Generates a human-readable phenotype name
//'
//' \code{make_phenotype_name} generates a human-readable phenotype name from a numeric representation
//' with the following encoding:
//' -1 : Drop marker
//' 0  : - (negative)
//' 1  : + (positive)
//' 2  : ++ (double positive)
//' More positive states can also be represented (e.g. 3, 4... : +++, ++++...)
//'
//' @param phenotype Numeric vector with phenotype encoding.
//' @param markers Names for markers. Must match the length of \code{phenotype}.
//' @return Phenotype name.
//'
//' @examples
//' make_phenotype_name(c(0,1,2,0),c("Marker1","Marker2","Marker3","Marker4"))
//' make_phenotype_name(c(0,1,-1),c("Marker1","Marker2","Marker3"))
//' 
//' @export
// [[Rcpp::export]]
String make_phenotype_name(IntegerVector phenotype, CharacterVector markers){
  
  vector<string> phen;
  
  for(int i = 0; i < phenotype.length(); i++){
    
    if(phenotype[i] == 0){
      
      string marker = Rcpp::as<string>(markers[i]); 
      phen.push_back(marker);
      phen.push_back("-");
      
    }else if(phenotype[i] > 0){
      
      string marker = Rcpp::as<string>(markers[i]); 
      phen.push_back(marker);
      for(int j = 0; j < phenotype[i]; j++) phen.push_back("+");
      
    }
  }
  
  stringstream ss;
  for_each(phen.begin(), phen.end(), [&ss] (const string& s) { ss << s; });
  
  return( ss.str() );
  
}


//' Count markers in phenotype
//'
//' \code{count_markers} counts the number of markers a phenotype has
//' given its numerical representation. Basically it returns the number of
//' values greater than -1 in an integer vector.
//'
//' @param phenotype Numeric vector with phenotype encoding.
//' @return Number of markers.
//'
//' @examples
//' count_markers(c(0,1,2,0))
//' count_markers(c(0,1,-1))
//'
//' @export
// [[Rcpp::export]]
int count_markers(IntegerVector phenotype){
  
  int n_markers = 0;  
  
  for(int i = 0; i<phenotype.length(); i++){
    if(phenotype[i] > -1) n_markers++;
  }
  
  return(n_markers);
}