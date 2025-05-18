#include <Rinternals.h>
#include <fcntl.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
// Define bcftools_munge as the name for the run function in munge.RBCFLIB.c
#define run bcftools_munge
#include "BCFToolsScore/munge.RBCFLIB.c"

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
  
  // Get number of arguments
  int user_n = length(args);
  nargs = user_n + 1; // +1 for the command name
  // Create array for user arguments plus command name
  char *c_args[nargs];
  // Set the first argument to the command name
  c_args[0] = "+munge"; // Prefix with + for bcftools plugin
  for (i = 0; i < user_n; i++) {
    c_args[i+1] = user_argv[i];
  }
  
  // Always print debugging information to diagnose the issue
  Rprintf("Running bcftools munge with %d arguments\n", nargs);
  for (i = 0; i < nargs; i++) {
    Rprintf("arg[%d]: %s\n", i, c_args[i]);
  }

  // Set environment variable for debugging if needed
  // putenv("BCFTOOLS_DEBUG=1");
  
  // now run bcftools_munge directly
  status = bcftools_munge(nargs, c_args);
  
  // Close file handles
  bcftools_close_stdout();
  bcftools_close_stderr();
  
  // Return status code with command attribute
  if (fd_stdout != -1 && fd_stdout != 1) close(fd_stdout); // Don't close stdout (fd 1)
  if (fd_stderr != -1 && fd_stderr != 2) close(fd_stderr); // Don't close stderr (fd 2)
  if (status == -1) error("bcftools_munge failed");

  // Build command character vector to return
  SEXP cmd = PROTECT(allocVector(STRSXP, user_n));
  for (i = 0; i < user_n; i++) {
    SET_STRING_ELT(cmd, i, mkChar(user_argv[i]));
  }
  // Wrap status as integer with 'command' attribute
  SEXP res = PROTECT(ScalarInteger(status));
  setAttrib(res, install("command"), cmd);
  UNPROTECT(2);
  return res;
}