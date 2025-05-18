#include "RBCFLib.h"
extern char* hts_version(void);
extern char* bcftools_version(void);
/* Function declarations */
SEXP RC_HTSLibVersion(void);
SEXP RC_BCFToolsVersion(void);

/* Function implementations */
/* 
 * Function to return htslib version as an R character string
 */
SEXP RC_HTSLibVersion(void) {
    SEXP result;
    const char *version = hts_version();
    
    PROTECT(result = mkString(version));
    UNPROTECT(1);
    
    return result;
}

/*
 * Function to return bcftools version as an R character string
 */
SEXP RC_BCFToolsVersion(void) {
    SEXP result;
    char *version = bcftools_version();
    
    PROTECT(result = mkString(version));
    UNPROTECT(1);
    
    return result;
}