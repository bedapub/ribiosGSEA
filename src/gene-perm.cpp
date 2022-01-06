#include "random"
#include "algorithm"
#include <Rcpp.h>

#ifdef _OPENMP
  #include <omp.h>  // This line won't add the library if you don't compile with -fopenmp option.
  #ifdef _MSC_VER
    // For Microsoft compiler
    #define OMP_FOR(n) __pragma(omp parallel for if(n>10))
  #else  // assuming "__GNUC__" is defined
    // For GCC compiler
    #define OMP_FOR(n) _Pragma("omp parallel for if(n>10)")
  #endif
#else
  #define omp_get_thread_num() 0
  #define OMP_FOR(n)
#endif

using namespace Rcpp;

typedef std::vector<int> IndexList;
typedef std::vector< std::vector<int> > IndexListSet;

int op_decrease (int i) { return --i;}

double subsetSum(const Rcpp::NumericVector &x,
		 const IndexList &index) {
  unsigned int i;
  double stat=0.0;
  for(i=0; i<index.size(); i++) {
    stat += x[ index[i] ]; 
  }
  return(stat);
}

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
RcppExport SEXP cpp_geneSetPerm(SEXP stats,
				SEXP rinds, // list
				SEXP Nsim) {
BEGIN_RCPP

  int i,j,ng;
  double psum;
  
  Rcpp::NumericVector xx(stats);
  Rcpp::NumericVector xp(Rcpp::clone(stats));
  Rcpp::List rindList(rinds);
  ng = rindList.size();
  std::vector<int> geneSetSizes(ng,0);
  
  IndexListSet indset;
  for(i=0; i<ng;i++) {
    indset.push_back(Rcpp::as< std::vector<int> >(rindList[i]));
    geneSetSizes[i]=indset[i].size();
    transform(indset[i].begin(), 
	      indset[i].end(), 
	      indset[i].begin(), op_decrease);
  }
  int nsim=Rcpp::as<int>(Nsim);

  std::vector<int> scounts(ng,0);
  Rcpp::NumericVector ps(ng);
  Rcpp::NumericVector gsStats(ng);
  
  for(i=0;i<ng;i++) {
    gsStats[i]=subsetSum(xx,indset[i]);
  }


  // https://meetingcpp.com/blog/items/stdrandom_shuffle-is-deprecated.html
  std::random_device rng;
  std::mt19937 urng(rng());

  for(i=0;i<nsim;i++) {
    std::shuffle(xp.begin(), xp.end(), urng);
    #pragma omp parallel for
    for(j=0; j<ng; j++) {
      psum=subsetSum(xp,indset[j]);
      if(psum>=gsStats[j])
        scounts[j]++;
    }
  }

  for(i=0;i<ng;i++) {
    if(geneSetSizes[i]==0) {
      ps[i]=1.0;
      gsStats[i] = NA_REAL;
    } else {
      ps[i]=(scounts[i]+1.0)/(nsim+1.0);
      gsStats[i] = gsStats[i]/geneSetSizes[i];
    }

  }
  
  Rcpp::DataFrame ret =
    Rcpp::DataFrame::create(Rcpp::Named("mean")=gsStats,
			    Rcpp::Named("geneSetSize")=geneSetSizes,
			    Rcpp::Named("p")=ps);
  
  
  return(ret);
END_RCPP
}
