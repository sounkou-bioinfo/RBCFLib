// this adappted from https://github.com/pysam-developers/pysam/blob/master/import/pysam.h
#ifndef BCFTOOLS_RBCFLIB_H
#define BCFTOOLS_RBCFLIB_H

#include <stdio.h>

#ifndef __has_attribute
#define __has_attribute(attribute) 0
#endif
#ifndef RBCFLIB_NORETURN
#if __has_attribute(__noreturn__) || __GNUC__ >= 3
#define RBCFLIB_NORETURN __attribute__((__noreturn__))
#else
#define RBCFLIB_NORETURN
#endif
#endif

/* Standard output and error streams for bcftools */
extern FILE * bcftools_stderr;
extern FILE * bcftools_stdout;
extern const char * bcftools_stdout_fn;

/* Set bcftools standard error to point to file descriptor
   Setting the stderr will close the previous stderr. */
FILE * bcftools_set_stderr(int fd);

/* Set bcftools standard output to point to file descriptor
   Setting the stdout will close the previous stdout. */
FILE * bcftools_set_stdout(int fd);

/* Set bcftools standard output to point to filename */
void bcftools_set_stdout_fn(const char * fn);

/* Close bcftools standard error and set to NULL */
void bcftools_close_stderr(void);

/* Close bcftools standard output and set to NULL */
void bcftools_close_stdout(void);

/* Write string to bcftools_stdout */
int bcftools_puts(const char *s);

/* Dispatch to appropriate bcftools command */
int bcftools_dispatch(int argc, char *argv[]);

/* Handle exit calls from bcftools */
void RBCFLIB_NORETURN bcftools_exit(int status);

/* Main entry point for bcftools */
extern int bcftools_main(int argc, char *argv[]);

/* These are only needed in C source, not R interface code */
#if !(defined CYTHON_ABI || defined CYTHON_HEX_VERSION || defined R_INTERFACE)

/* Several non-static function names are used in both samtools and bcftools.
   Both libraries may be loaded simultaneously, leading to collisions.
   #define these names so the actual symbol names include distinct prefixes. */
#define main_consensus bcftools_main_consensus
#define main_reheader bcftools_main_reheader
#define bam_smpl_init bcftools_bam_smpl_init
#define bam_smpl_destroy bcftools_bam_smpl_destroy
#define read_file_list bcftools_read_file_list


#endif

#endif /* BCFTOOLS_RBCFLIB_H */