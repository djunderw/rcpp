#include <cstdlib>
#include <iostream>
#include <vector>
#include <algorithm>
#include <Rcpp.h>

using namespace std;

/*!
 *  Efficient version of the "Mden_i_derv()" function 
 * 	from Jung-In's R code.
 *
 * 	@param i integer index
 * 	@param j integer index
 * 	@param length_ni number of events for each i
 * 	@param dat the matrix of data
 * 	@param x unknown vector
 * 	@param zz unknown matrix
 */
RcppExport SEXP Mden_i_derv_cpp(  	SEXP i,
									SEXP j,
									SEXP length_ni,
									SEXP dat,
									SEXP x,
									SEXP zz			) {

	Rcpp::NumericVector vec1(x);
	Rcpp::NumericVector vec2(zz);
	Rcpp::NumericVector product(vec1.size());
	Rcpp::NumericVector dotproduct(1);

	if(vec1.size()==vec2.size()) {
		for(int idx=0; idx<vec1.size(); ++idx) {
			product[idx]=vec1[idx]*vec2[idx];
		}
	} else {
		cout << "x and zz have unequal dimensions!\n";
	}

	dotproduct[0]=0.0;
	
	for(int idx=0; idx<product.size(); ++idx)
		dotproduct[0] += product[idx];

	return dotproduct;
}
