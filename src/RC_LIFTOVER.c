#include <Rinternals.h>
#include <fcntl.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "liftover.RBCFLIB.c"

SEXP RC_bcftools_liftover(
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
  char* cmd_str = "liftover";
  
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
  
  
  // Create empty BCF headers for input and output
  bcf_hdr_t *in_hdr = bcf_hdr_init("r");
  bcf_hdr_t *out_hdr = bcf_hdr_init("w");
  
  if (!in_hdr || !out_hdr) {
    error("Failed to allocate BCF headers");
  }
  
  // Run the command
  status = run_liftover(nargs, user_argv, in_hdr, out_hdr);
  
  // Clean up
  bcf_hdr_destroy(in_hdr);
  bcf_hdr_destroy(out_hdr);
  
  if (fd_stdout != -1) close(fd_stdout);
  if (fd_stderr != -1) close(fd_stderr);
  
  // Return status and command string as attribute
  SEXP result = PROTECT(ScalarInteger(status));
  char full_cmd[1000] = "bcftools liftover";
  for (i = 0; i < nargs; i++) {
    strcat(full_cmd, " ");
    strcat(full_cmd, user_argv[i]);
  }
  setAttrib(result, install("command"), mkString(full_cmd));
  UNPROTECT(1);
  
  return result;
}
