#include <R.h>
#include <Rinternals.h>
#include "bcftools.RBCFLIB.h"
#include <fcntl.h>
extern void bcftools_close_stdout(void);
extern void bcftools_close_stderr(void);
SEXP RC_bcftools_run(SEXP args, SEXP capture_stdout, SEXP capture_stderr, 
                    SEXP stdout_file, SEXP stderr_file) {
  int i, status, nargs;
  const char *std_out, *std_err;
  int fd_stdout = -1, fd_stderr = -1;
  
  // Parse R arguments
  nargs = length(args);
  char *argv[nargs];
  
  for (i = 0; i < nargs; i++) {
    argv[i] = (char*) CHAR(STRING_ELT(args, i));
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
  status = bcftools_dispatch(nargs, argv);
  
  // Close file handles
  bcftools_close_stdout();
  bcftools_close_stderr();
  
  // Return status code
  return ScalarInteger(status);
}