#ifndef RBCFLIB_H
#define RBCFLIB_H
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <htslib/hts.h>
#include <bcftools.h>
/* Version function declarations */
SEXP RC_HTSLibVersion(void);
SEXP RC_BCFToolsVersion(void);

#endif /* RBCFLIB_H */