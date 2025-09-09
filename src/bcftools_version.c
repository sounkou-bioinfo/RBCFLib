#include "./bcftools-1.22/version.h"
#define SCORE_VERSION "2025-08-19"


/* Minimal bcftools version function for R package */
char *bcftools_version(void)
{
    return BCFTOOLS_VERSION;
}
char *bcftools_score_version(void)
{
    return SCORE_VERSION; // Placeholder version
}