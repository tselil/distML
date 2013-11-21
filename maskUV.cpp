#include <Rcpp.h>
using namespace Rcpp;

// Given matrix factors U and V, and a list of entries (is,js) returns
// a list of entries in maskUV so that maskUV(k) = UV'(is,js)
// requires maskUV to be pre-allocated in the R code that calls this
//
// U is an m-by-r matrix
// V is an n-by-r matrix

// [Rcpp::export]]
NumericVector maskUV(NumericMatrix U, NumericMatrix V, IntegerVector is, IntegerVector js)
{
	// Get the length of the entries list and the rank of UV'
	int l = is.size();
	int r = U.ncol();
	
	// Initialize the output vector to all zeros
	NumericVector maskUV(l,0.0);
	
	// Loop over non-zero entries and compute output vector
	int i = is(0);
	int j = js(0);
	for(int n = 0; n < l; n++)
	{
		i = is(n);
		j = js(n);
		maskUV(n) = 0;
		for(int k = 0; k < r; k++)
		{
			maskUV(n) += U(i,k)*V(j,k);
		}
	}
	return maskUV;
}
