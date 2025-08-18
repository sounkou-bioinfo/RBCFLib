#include "RBCFLib.h"
#include "htslib/hts.h"
#include "htslib/faidx.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

/* Function declarations */
SEXP RC_HTSLibVersion(void);
SEXP RC_BCFToolsVersion(void);
SEXP RC_FaidxIndexFasta(SEXP fasta_path);
SEXP RC_FaidxFetchRegion(SEXP fasta_path, SEXP seqname, SEXP start, SEXP end);
extern char *bcftools_version(void);

/* Binary path storage */
static char *cached_bcftools_path = NULL;

/* 
 * Function to get bcftools binary path
 */
const char* BCFToolsBinaryPath(void) {
    if (cached_bcftools_path == NULL) {
        // Use system.file to get the bin directory  
        SEXP call = PROTECT(lang3(install("system.file"),
                                 mkString("bin"),
                                 ScalarString(mkChar("package")), 
                                 mkString("RBCFLib")));
        SEXP pkg_dir = PROTECT(eval(call, R_GlobalEnv));
        
        if (TYPEOF(pkg_dir) == STRSXP && length(pkg_dir) > 0) {
            const char *bin_dir = CHAR(STRING_ELT(pkg_dir, 0));
            // Construct path: <bin_dir>/bcftools
            size_t path_len = strlen(bin_dir) + strlen("/bcftools") + 1;
            cached_bcftools_path = (char*)malloc(path_len);
            if (cached_bcftools_path != NULL) {
                snprintf(cached_bcftools_path, path_len, "%s/bcftools", bin_dir);
            }
        }
        UNPROTECT(2);
        
        if (cached_bcftools_path == NULL) {
            error("Failed to get bcftools binary path");
        }
    }
    
    return cached_bcftools_path;
}

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

/*
 * Function to index a FASTA file using htslib faidx and return the index path
 */
SEXP RC_FaidxIndexFasta(SEXP fasta_path) {
    const char *path = CHAR(STRING_ELT(fasta_path, 0));
    SEXP result;
    
    int ret = fai_build(path);
    if (ret != 0) {
        error("Failed to index FASTA file: %s", path);
    }
    
    // Create the fai path by appending .fai
    char *fai_path = R_alloc(strlen(path) + 5, sizeof(char));
    strcpy(fai_path, path);
    strcat(fai_path, ".fai");
    
    PROTECT(result = mkString(fai_path));
    UNPROTECT(1);
    
    return result;
}

/*
 * Function to fetch a sequence region from a FASTA file using htslib faidx
 */
SEXP RC_FaidxFetchRegion(SEXP fasta_path, SEXP seqname, SEXP start, SEXP end) {
    const char *path = CHAR(STRING_ELT(fasta_path, 0));
    const char *seq_name = CHAR(STRING_ELT(seqname, 0));
    int start_pos = INTEGER(start)[0];
    int end_pos = INTEGER(end)[0];
    
    SEXP result;
    
    faidx_t *fai = fai_load(path);
    if (!fai) {
        error("Failed to load FASTA index for %s", path);
    }
    
    int seq_len = 0;
    // htslib uses 0-based, end-inclusive coordinates, R is 1-based
    char *seq = faidx_fetch_seq(fai, seq_name, start_pos - 1, end_pos - 1, &seq_len);
    
    if (!seq) {
        fai_destroy(fai);
        error("Failed to fetch sequence for region %s:%d-%d", seq_name, start_pos, end_pos);
    }
    
    PROTECT(result = allocVector(STRSXP, 1));
    SET_STRING_ELT(result, 0, mkChar(seq));
    
    free(seq);
    fai_destroy(fai);
    UNPROTECT(1);
    
    return result;
}
