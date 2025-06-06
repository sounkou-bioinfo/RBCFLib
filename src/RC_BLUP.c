#include <Rinternals.h>
#include <fcntl.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "blup.RBCFLIB.c"


SEXP RC_bcftools_blup(
                  SEXP args, 
                  SEXP capture_stdout,
                  SEXP capture_stderr, 
                  SEXP stdout_file,
                  SEXP stderr_file) {

  /* Reset getopt()/getopt_long() processing. */
#if defined __GLIBC__
  optind = 0;
#elif defined _OPTRESET || defined _OPTRESET_DECLARED
  optreset = optind = 1;
#else
  optind = 1;
#endif

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
  char* cmd_str = "blup";
  
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
      error("Could not open stderr file for writing");
    }
    bcftools_set_stderr(fd_stderr);
  } else {
    // Redirect to /dev/null if we don't want to capture stderr
    fd_stderr = open("/dev/null", O_WRONLY);
    if (fd_stderr == -1) {
      error("Could not open /dev/null for writing");
    }
    bcftools_set_stderr(fd_stderr);
  }
  
  // Call the main function with our arguments
  status = run_blup(nargs, user_argv);
  
  // Close file handles
  bcftools_close_stdout();
  bcftools_close_stderr();
  
  // Return status code with command attribute
  if (fd_stdout != -1 && fd_stdout != 1) close(fd_stdout); // Don't close stdout (fd 1)
  if (fd_stderr != -1 && fd_stderr != 2) close(fd_stderr); // Don't close stderr (fd 2)
  if (status == -1) error("bcftools blup failed");
  
  // Return status and command string as attribute
  SEXP result = PROTECT(ScalarInteger(status));
  char full_cmd[1000] = "bcftools blup";
  for (i = 0; i < nargs; i++) {
    strcat(full_cmd, " ");
    strcat(full_cmd, user_argv[i]);
  }
  setAttrib(result, install("command"), mkString(full_cmd));
  UNPROTECT(1);
  
  return result;
}
