#include <R.h>
#include <Rinternals.h>
#include "bcftools.RBCFLIB.h"
#include <fcntl.h>
#include <stdlib.h>
#include <unistd.h>          // for close()

SEXP RC_bcftools_run(SEXP args, SEXP capture_stdout, SEXP capture_stderr, 
                    SEXP stdout_file, SEXP stderr_file) {


  int i, status, nargs;
  const char *std_out, *std_err;
  int fd_stdout = -1, fd_stderr = -1;
  
  // Parse R arguments
  nargs = length(args);
  char *user_argv[nargs];
  
  for (i = 0; i < nargs; i++) {
    user_argv[i] = (char*) CHAR(STRING_ELT(args, i));
  }
  
  // Setup stdout redirection if requested
  if (asLogical(capture_stdout)) {
    std_out = CHAR(STRING_ELT(stdout_file, 0));
    fd_stdout = open(std_out, O_WRONLY|O_CREAT|O_TRUNC, 0666);
    if (fd_stdout == -1) {
      error("Could not open stdout file for writing");
    }
    bcftools_set_stdout(fd_stdout);
  } else {
    bcftools_set_stdout(1); // Standard stdout
  }
  
  // Setup stderr redirection if requested
  if (asLogical(capture_stderr)) {
    std_err = CHAR(STRING_ELT(stderr_file, 0));
    fd_stderr = open(std_err, O_WRONLY|O_CREAT|O_TRUNC, 0666);
    if (fd_stderr == -1) {
      if (fd_stdout != -1) close(fd_stdout);
      error("Could not open stderr file for writing");
    }
    bcftools_set_stderr(fd_stderr);
  } else {
    bcftools_set_stderr(2); // Standard stderr
  }
  
  // Call the bcftools dispatcher
  Rprintf("Running bcftools with %d arguments\n", nargs);
  for (i = 0; i < nargs; i++) {
    Rprintf("argv[%d]: %s\n", i, user_argv[i]);
  }

  // build an argv[] with an extra slot for the program name
  int user_n = length(args);
  int bc_n   = user_n + 1;            // +1 for "bcftools"
  char *bc_argv[bc_n];

  // argv[0] must be the program name
  bc_argv[0] = "bcftools";

  // shift user args into argv[1..]
  for (i = 0; i < user_n; i++) {
    bc_argv[i+1] = user_argv[i];
  }

  // now dispatch
  status = bcftools_dispatch(bc_n, bc_argv);
  
  // Close file handles
  bcftools_close_stdout();
  bcftools_close_stderr();
  
  // Return status code
  if (fd_stdout != -1) {
    close(fd_stdout);
  }
  if (fd_stderr != -1) {
    close(fd_stderr);
  }
  if (status == -1) {
    error("bcftools_dispatch failed");
  }
  if (status == 0) {
    Rprintf("bcftools completed successfully\n");
  } else {
    Rprintf("bcftools failed with status %d\n", status);
  }

  return ScalarInteger(status);
}