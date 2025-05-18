#include <Rinternals.h>
#include <fcntl.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "munge.RBCFLIB.c"

SEXP RC_bcftools_munge(
                  SEXP args, 
                  SEXP capture_stdout,
                  SEXP capture_stderr, 
                  SEXP stdout_file,
                  SEXP stderr_file,
                  SEXP is_usage) {


  int i, status, nargs;
  const char *std_out, *std_err;
  int fd_stdout = -1, fd_stderr = -1;
  
  // Parse R arguments
  nargs = length(args);
  char *user_argv[nargs];
  
  for (i = 0; i < nargs; i++) {
    user_argv[i] = (char*) CHAR(STRING_ELT(args, i));
  }
  
  // Get the command string
  char* cmd_str = "bcftools+munge";
  
  // Setup stdout redirection if requested
  if (asLogical(capture_stdout)) {
    std_out = CHAR(STRING_ELT(stdout_file, 0));
    fd_stdout = open(std_out, O_WRONLY|O_CREAT|O_TRUNC, 0666);
    if (fd_stdout == -1) {
      error("Could not open stdout file for writing");
    }
    bcftools_set_stdout(fd_stdout);
  } else {
    // Redirect to /dev/null if we don't want to capture stdout
    fd_stdout = open("/dev/null", O_WRONLY);
    if (fd_stdout == -1) {
      error("Could not open /dev/null for writing");
    }
    bcftools_set_stdout(fd_stdout);
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
    // Redirect stderr to /dev/null if we don't want to capture it
    fd_stderr = open("/dev/null", O_WRONLY);
    if (fd_stderr == -1) {
      if (fd_stdout != -1) close(fd_stdout);
      error("Could not open /dev/null for writing");
    }
    bcftools_set_stderr(fd_stderr);
  }
  
  // Build an argv[] with slots for "bcftools" and the command
  int user_n = length(args);
  int bc_n   = user_n + 1;            // +1 for "bcftools", +1 for command
  char *bc_argv[bc_n];

  // argv[0] must be the program name
  bc_argv[0] = cmd_str;

  // shift user args into argv[2..]
  for (i = 0; i < user_n; i++) {
    bc_argv[i+1] = user_argv[i];
  }
  
  // For debugging, can be enabled by an environment variable or verbose flag
  if (getenv("RBCFLIB_DEBUG") != NULL) {
    Rprintf("Running bcftools %s with %d arguments\n", cmd_str, nargs);
    for (i = 0; i < bc_n; i++) {
      Rprintf("argv[%d]: %s\n", i, bc_argv[i]);
    }
  }

  // now dispatch
  status = run_munge(bc_n, bc_argv);
  
  // Close file handles
  bcftools_close_stdout();
  bcftools_close_stderr();
  
  // Return status code with command attribute
  if (fd_stdout != -1 && fd_stdout != 1) close(fd_stdout); // Don't close stdout (fd 1)
  if (fd_stderr != -1 && fd_stderr != 2) close(fd_stderr); // Don't close stderr (fd 2)
  if (status == -1) error("bcftools_dispatch failed");

  // Build command character vector to return
  SEXP cmd = PROTECT(allocVector(STRSXP, bc_n));
  for (i = 0; i < bc_n; i++) {
    SET_STRING_ELT(cmd, i, mkChar(bc_argv[i]));
  }
  // Wrap status as integer with 'command' attribute
  SEXP res = PROTECT(ScalarInteger(status));
  setAttrib(res, install("command"), cmd);
  UNPROTECT(2);
  return res;
}