#include "./bcftools-1.22/version.h"

/* Minimal bcftools version function for R package */
char *bcftools_version(void)
{
    return BCFTOOLS_VERSION;
}