#include <Rinternals.h>
#include <unistd.h>
#include <sys/wait.h>
#include <string.h>
#include <stdio.h>
#include <fcntl.h>
#include <stdlib.h>
#include <signal.h>
#include "bcftools.RBCFLIB.h"

// Helper: convert SEXP of strings to char** (NULL-terminated)
static char **sexp_to_argv(SEXP args, SEXP command, int *argc_out) {
    int nargs = length(args);
    int argc = nargs + 2;
    char **argv = (char **)malloc((argc + 1) * sizeof(char *));
    if (argv == NULL) {
        error("Memory allocation failed");
    }
    
    argv[0] = strdup("bcftools");
    argv[1] = strdup(CHAR(STRING_ELT(command, 0)));
    for (int i = 0; i < nargs; i++)
        argv[i + 2] = strdup(CHAR(STRING_ELT(args, i)));
    argv[argc] = NULL;
    
    if (argc_out) *argc_out = argc;
    return argv;
}

// Helper: free memory allocated for argv
static void free_argv(char **argv) {
    if (!argv) return;
    for (int i = 0; argv[i] != NULL; i++) {
        free(argv[i]);
    }
    free(argv);
}

/**
 * Pipe the output of one bcftools command into another
 * 
 * @param command1 First bcftools command
 * @param args1 Arguments for first command
 * @param command2 Second bcftools command
 * @param args2 Arguments for second command
 * @param capture_stdout Whether to capture stdout from the second command
 * @param capture_stderr Whether to capture stderr from both commands
 * @param stdout_file File to capture stdout (for second command)
 * @param stderr_file File to capture stderr (for both commands)
 * 
 * @return Integer vector of exit statuses with 'command' attribute containing combined command
 */
SEXP RC_bcftools_pipe(
    SEXP command1, SEXP args1,
    SEXP command2, SEXP args2,
    SEXP capture_stdout,
    SEXP capture_stderr, 
    SEXP stdout_file,
    SEXP stderr_file
) {
    int pipefd[2];
    const char *std_out, *std_err;
    int fd_stdout = -1, fd_stderr = -1;
    int argc1, argc2;
    char **argv1 = NULL, **argv2 = NULL;
    SEXP res = R_NilValue, cmd = R_NilValue;
    int do_capture_stdout, do_capture_stderr;
    
    // Extract capture flags before the fork
    do_capture_stdout = asLogical(capture_stdout);
    do_capture_stderr = asLogical(capture_stderr);
    
    // Setup pipe
    if (pipe(pipefd) == -1) {
        error("pipe() failed");
    }
    
    // Setup stdout redirection if requested (for second command only)
    if (do_capture_stdout) {
        std_out = CHAR(STRING_ELT(stdout_file, 0));
        fd_stdout = open(std_out, O_WRONLY|O_CREAT|O_TRUNC, 0666);
        if (fd_stdout == -1) {
            close(pipefd[0]);
            close(pipefd[1]);
            error("Could not open stdout file for writing: %s", std_out);
        }
    } else {
        // Redirect to /dev/null if we don't want to capture stdout
        fd_stdout = open("/dev/null", O_WRONLY);
        if (fd_stdout == -1) {
            close(pipefd[0]);
            close(pipefd[1]);
            error("Could not open /dev/null for writing");
        }
    }
    
    // Setup stderr redirection if requested (for both commands)
    if (do_capture_stderr) {
        std_err = CHAR(STRING_ELT(stderr_file, 0));
        fd_stderr = open(std_err, O_WRONLY|O_CREAT|O_TRUNC, 0666);
        if (fd_stderr == -1) {
            close(pipefd[0]);
            close(pipefd[1]);
            if (fd_stdout != -1) close(fd_stdout);
            error("Could not open stderr file for writing");
        }
    } else {
        // Redirect stderr to /dev/null if we don't want to capture it
        fd_stderr = open("/dev/null", O_WRONLY);
        if (fd_stderr == -1) {
            close(pipefd[0]);
            close(pipefd[1]);
            if (fd_stdout != -1) close(fd_stdout);
            error("Could not open /dev/null for writing");
        }
    }
    
    // Convert R args to C argv arrays - allocating new memory
    argv1 = sexp_to_argv(args1, command1, &argc1);
    argv2 = sexp_to_argv(args2, command2, &argc2);
    
    // Debug output if requested
    if (getenv("RBCFLIB_DEBUG") != NULL) {
        Rprintf("Piping: bcftools %s | bcftools %s\n", 
                CHAR(STRING_ELT(command1, 0)), CHAR(STRING_ELT(command2, 0)));
        
        Rprintf("Command 1 arguments:\n");
        for (int i = 0; i < argc1; i++) {
            Rprintf("  argv[%d]: %s\n", i, argv1[i]);
        }
        
        Rprintf("Command 2 arguments:\n");
        for (int i = 0; i < argc2; i++) {
            Rprintf("  argv[%d]: %s\n", i, argv2[i]);
        }
    }
    
    // Create first child process
    pid_t pid1 = fork();
    if (pid1 < 0) {
        free_argv(argv1);
        free_argv(argv2);
        close(pipefd[0]);
        close(pipefd[1]);
        if (fd_stdout != -1) close(fd_stdout);
        if (fd_stderr != -1) close(fd_stderr);
        error("fork() failed for first command");
    }

    if (pid1 == 0) {
        // First child: run first command, output to pipe
        close(pipefd[0]);  // Close read end
        
        // Redirect stdout to pipe
        if (dup2(pipefd[1], STDOUT_FILENO) == -1) {
            perror("dup2 stdout");
            exit(1);
        }
        
        // Redirect stderr if requested
        if (do_capture_stderr) {
            if (dup2(fd_stderr, STDERR_FILENO) == -1) {
                perror("dup2 stderr");
                exit(1);
            }
        }
        
        // Set bcftools stdout/stderr
        bcftools_set_stdout(pipefd[1]);
        if (do_capture_stderr) {
            bcftools_set_stderr(fd_stderr);
        }
        
        // Run bcftools
        int status = bcftools_dispatch(argc1, argv1);
        
        // Clean up and exit
        bcftools_close_stdout();
        bcftools_close_stderr();
        free_argv(argv1);
        free_argv(argv2);
        close(pipefd[1]); // Close after running
        exit(status);
    }

    // Create second child process
    pid_t pid2 = fork();
    if (pid2 < 0) {
        free_argv(argv1);
        free_argv(argv2);
        close(pipefd[0]);
        close(pipefd[1]);
        if (fd_stdout != -1) close(fd_stdout);
        if (fd_stderr != -1) close(fd_stderr);
        kill(pid1, SIGTERM); // Kill first child
        error("fork() failed for second command");
    }

    if (pid2 == 0) {
        // Second child: run second command, input from pipe
        close(pipefd[1]);  // Close write end
        
        // Redirect stdin to pipe
        if (dup2(pipefd[0], STDIN_FILENO) == -1) {
            perror("dup2 stdin");
            exit(1);
        }
        
        // Redirect stdout if requested
        if (dup2(fd_stdout, STDOUT_FILENO) == -1) {
            perror("dup2 stdout");
            exit(1);
        }
        
        // Redirect stderr if requested
        if (do_capture_stderr) {
            if (dup2(fd_stderr, STDERR_FILENO) == -1) {
                perror("dup2 stderr");
                exit(1);
            }
        }
        
        // Set bcftools stdout/stderr
        bcftools_set_stdout(fd_stdout);
        if (do_capture_stderr) {
            bcftools_set_stderr(fd_stderr);
        }
        
        // Run bcftools
        int status = bcftools_dispatch(argc2, argv2);
        
        // Clean up and exit
        bcftools_close_stdout();
        bcftools_close_stderr();
        free_argv(argv1);
        free_argv(argv2);
        close(pipefd[0]); // Close after running
        exit(status);
    }

    // Parent: close pipe ends immediately
    close(pipefd[0]);
    close(pipefd[1]);
    
    // Wait for children to finish
    int status1, status2;
    waitpid(pid1, &status1, 0);
    waitpid(pid2, &status2, 0);
    
    // Now build the command attribute for the result
    int total_args = argc1 + argc2 + 1;  // +1 for the pipe symbol
    PROTECT(cmd = allocVector(STRSXP, total_args));
    
    int j = 0;
    // First command
    for (int i = 0; i < argc1; i++) {
        SET_STRING_ELT(cmd, j++, mkChar(argv1[i]));
    }
    // Add pipe symbol
    SET_STRING_ELT(cmd, j++, mkChar("|"));
    // Second command
    for (int i = 0; i < argc2; i++) {
        SET_STRING_ELT(cmd, j++, mkChar(argv2[i]));
    }
    
    // Clean up parent resources
    free_argv(argv1);
    free_argv(argv2);
    if (fd_stdout != -1 && fd_stdout != STDOUT_FILENO) close(fd_stdout);
    if (fd_stderr != -1 && fd_stderr != STDERR_FILENO) close(fd_stderr);
    
    // Create result
    PROTECT(res = allocVector(INTSXP, 2));
    INTEGER(res)[0] = WIFEXITED(status1) ? WEXITSTATUS(status1) : -1;
    INTEGER(res)[1] = WIFEXITED(status2) ? WEXITSTATUS(status2) : -1;
    
    // Set attribute with combined command
    setAttrib(res, install("command"), cmd);
    
    UNPROTECT(2);
    return res;
}
