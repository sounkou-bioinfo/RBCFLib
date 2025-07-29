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
 * Execute a pipeline of bcftools commands
 * 
 * @param commands List of bcftools commands
 * @param args List of arguments for each command
 * @param n_commands Number of commands in the pipeline
 * @param capture_stdout Whether to capture stdout from the last command
 * @param capture_stderr Whether to capture stderr from all commands
 * @param stdout_file File to capture stdout (for last command)
 * @param stderr_file File to capture stderr (for all commands)
 * 
 * @return Integer vector of exit statuses with 'command' attribute containing combined command
 */
SEXP RC_bcftools_pipeline(
    SEXP commands, SEXP args,
    SEXP n_commands,
    SEXP capture_stdout,
    SEXP capture_stderr, 
    SEXP stdout_file,
    SEXP stderr_file
) {
    int num_commands = asInteger(n_commands);
    int pipes[num_commands - 1][2]; // pipes[i] connects command i and i+1
    pid_t pids[num_commands];       // process ids for each command
    int statuses[num_commands];      // exit statuses for each command
    int argc_values[num_commands];   // argc values for each command
    char ***argv_values = NULL;      // argv values for each command
    const char *std_out, *std_err;
    int fd_stdout = -1, fd_stderr = -1;
    SEXP res = R_NilValue, cmd = R_NilValue;
    int do_capture_stdout, do_capture_stderr;
    
    // Check for valid number of commands
    if (num_commands < 1) {
        error("At least one command is required");
    }
    
    // Extract capture flags before the fork
    do_capture_stdout = asLogical(capture_stdout);
    do_capture_stderr = asLogical(capture_stderr);
    
    // Setup stderr redirection if requested (for all commands)
    if (do_capture_stderr) {
        std_err = CHAR(STRING_ELT(stderr_file, 0));
        fd_stderr = open(std_err, O_WRONLY|O_CREAT|O_TRUNC, 0666);
        if (fd_stderr == -1) {
            error("Could not open stderr file for writing");
        }
    } else {
        // Redirect stderr to /dev/null if we don't want to capture it
        fd_stderr = open("/dev/null", O_WRONLY);
        if (fd_stderr == -1) {
            error("Could not open /dev/null for writing");
        }
    }
    
    // Setup stdout redirection if requested (for last command only)
    if (do_capture_stdout) {
        std_out = CHAR(STRING_ELT(stdout_file, 0));
        fd_stdout = open(std_out, O_WRONLY|O_CREAT|O_TRUNC, 0666);
        if (fd_stdout == -1) {
            if (fd_stderr != -1) close(fd_stderr);
            error("Could not open stdout file for writing: %s", std_out);
        }
    } else {
        // Redirect to /dev/null if we don't want to capture stdout
        fd_stdout = open("/dev/null", O_WRONLY);
        if (fd_stdout == -1) {
            if (fd_stderr != -1) close(fd_stderr);
            error("Could not open /dev/null for writing");
        }
    }
    
    // Allocate memory for argv arrays
    argv_values = (char ***)malloc(num_commands * sizeof(char **));
    if (argv_values == NULL) {
        if (fd_stdout != -1) close(fd_stdout);
        if (fd_stderr != -1) close(fd_stderr);
        error("Memory allocation failed for argv_values");
    }
    
    // Convert R args to C argv arrays for each command
    for (int i = 0; i < num_commands; i++) {
        SEXP cmd_i = VECTOR_ELT(commands, i);
        SEXP args_i = VECTOR_ELT(args, i);
        argv_values[i] = sexp_to_argv(args_i, cmd_i, &argc_values[i]);
    }
    
    // Debug output if requested
    if (getenv("RBCFLIB_DEBUG") != NULL) {
        Rprintf("Piping %d commands:\n", num_commands);
        for (int i = 0; i < num_commands; i++) {
            Rprintf("Command %d: bcftools %s\n", i+1, CHAR(STRING_ELT(VECTOR_ELT(commands, i), 0)));
            Rprintf("  Arguments:\n");
            for (int j = 0; j < argc_values[i]; j++) {
                Rprintf("    argv[%d]: %s\n", j, argv_values[i][j]);
            }
        }
    }
    
    // Create pipes between commands
    for (int i = 0; i < num_commands - 1; i++) {
        if (pipe(pipes[i]) == -1) {
            // Clean up already created pipes
            for (int j = 0; j < i; j++) {
                close(pipes[j][0]);
                close(pipes[j][1]);
            }
            // Free argv arrays
            for (int j = 0; j < num_commands; j++) {
                free_argv(argv_values[j]);
            }
            free(argv_values);
            if (fd_stdout != -1) close(fd_stdout);
            if (fd_stderr != -1) close(fd_stderr);
            error("pipe() creation failed");
        }
    }
    
    // Create child processes for each command
    for (int i = 0; i < num_commands; i++) {
        pids[i] = fork();
        if (pids[i] < 0) {
            // Error on fork, clean up
            for (int j = 0; j < num_commands - 1; j++) {
                close(pipes[j][0]);
                close(pipes[j][1]);
            }
            // Kill already created children
            for (int j = 0; j < i; j++) {
                kill(pids[j], SIGTERM);
            }
            // Free argv arrays
            for (int j = 0; j < num_commands; j++) {
                free_argv(argv_values[j]);
            }
            free(argv_values);
            if (fd_stdout != -1) close(fd_stdout);
            if (fd_stderr != -1) close(fd_stderr);
            error("fork() failed for command %d", i+1);
        }
        
        if (pids[i] == 0) {
            // Child process for command i
            
            // Setup stdin (for all but first command)
            if (i > 0) {
                // Take input from previous pipe
                if (dup2(pipes[i-1][0], STDIN_FILENO) == -1) {
                    perror("dup2 stdin");
                    exit(1);
                }
            }
            
            // Setup stdout
            if (i < num_commands - 1) {
                // Pipe output to next command
                if (dup2(pipes[i][1], STDOUT_FILENO) == -1) {
                    perror("dup2 stdout");
                    exit(1);
                }
                bcftools_set_stdout(pipes[i][1]);
            } else {
                // Last command - write to file or stdout
                if (dup2(fd_stdout, STDOUT_FILENO) == -1) {
                    perror("dup2 stdout");
                    exit(1);
                }
                bcftools_set_stdout(fd_stdout);
            }
            
            // Setup stderr (common to all commands)
            if (do_capture_stderr) {
                if (dup2(fd_stderr, STDERR_FILENO) == -1) {
                    perror("dup2 stderr");
                    exit(1);
                }
                bcftools_set_stderr(fd_stderr);
            }
            
            // Close all pipe ends not needed by this process
            for (int j = 0; j < num_commands - 1; j++) {
                // Close all pipe ends except:
                // - stdin for this process (pipes[i-1][0] if i > 0) 
                // - stdout for this process (pipes[i][1] if i < num_commands-1)
                if (j == i - 1 && i > 0) {
                    close(pipes[j][1]); // Keep read end, close write end
                } else if (j == i && i < num_commands - 1) {
                    close(pipes[j][0]); // Keep write end, close read end
                } else {
                    // Close both ends for pipes not directly connected to this process
                    close(pipes[j][0]);
                    close(pipes[j][1]);
                }
            }
            
            // Run bcftools
            int status = bcftools_dispatch(argc_values[i], argv_values[i]);
            
            // Clean up bcftools streams
            bcftools_close_stdout();
            bcftools_close_stderr();
              
            // Free all argv arrays
            for (int j = 0; j < num_commands; j++) {
                free_argv(argv_values[j]);
            }
            free(argv_values);
            
            exit(status);
        }
    }
    
    // Parent: close all pipe ends
    for (int i = 0; i < num_commands - 1; i++) {
        close(pipes[i][0]);
        close(pipes[i][1]);
    }
    
    // Wait for all children to finish
    for (int i = 0; i < num_commands; i++) {
        int status;
        waitpid(pids[i], &status, 0);
        statuses[i] = WIFEXITED(status) ? WEXITSTATUS(status) : -1;
    }
    
    // Now build the command attribute for the result
    // Count total arguments across all commands
    int total_args = 0;
    for (int i = 0; i < num_commands; i++) {
        total_args += argc_values[i];
        if (i < num_commands - 1) total_args++; // Add for pipe symbol
    }
    
    PROTECT(cmd = allocVector(STRSXP, total_args));
    
    int cmd_idx = 0;
    // Add each command and its arguments
    for (int i = 0; i < num_commands; i++) {
        // Add command and its arguments
        for (int j = 0; j < argc_values[i]; j++) {
            SET_STRING_ELT(cmd, cmd_idx++, mkChar(argv_values[i][j]));
        }
        
        // Add pipe symbol between commands
        if (i < num_commands - 1) {
            SET_STRING_ELT(cmd, cmd_idx++, mkChar("|"));
        }
    }
    
    // Clean up parent resources
    for (int i = 0; i < num_commands; i++) {
        free_argv(argv_values[i]);
    }
    free(argv_values);
    if (fd_stdout != -1 && fd_stdout != STDOUT_FILENO) close(fd_stdout);
    if (fd_stderr != -1 && fd_stderr != STDERR_FILENO) close(fd_stderr);
    
    // Create result
    PROTECT(res = allocVector(INTSXP, num_commands));
    for (int i = 0; i < num_commands; i++) {
        INTEGER(res)[i] = statuses[i];
    }
    //Rprintf("BCFtools pipeline executed with %d commands\n", num_commands);
    // Set attribute with combined command
    setAttrib(res, install("command"), cmd);
    
    UNPROTECT(2);
    return res;
}
