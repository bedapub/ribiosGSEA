/* Common definitions of ribiosGSEA */
#ifndef RIBIOS_GSEA_H
#define RIBIOS_GSEA_H

#include "Rcpp.h"

RcppExport SEXP cpp_geneSetPerm(SEXP, SEXP, SEXP);
RcppExport SEXP list2mat(SEXP);

#endif
