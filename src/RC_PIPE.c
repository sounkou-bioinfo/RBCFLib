#include <Rinternals.h>
#include <R.h>

#ifdef _WIN32
#include <windows.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "RBCFLib.h"
#endif

#ifndef _WIN32
#include <unistd.h>
#include <sys/wait.h>
#include <string.h>
#include <stdio.h>
#include <fcntl.h>
#include <stdlib.h>
#include <signal.h>
#include <errno.h>
#include "RBCFLib.h"

// Global flag to track if SIGPIPE has been handled
// FIX ME: this is a global state, should be managed
// should probably be removed
static int sigpipe_handled = 0;

// Helper: Setup SIGPIPE handling for pipeline operations
// this is needed to prevent subprocesses getting terminated by SIGPIPE
// which is may be captured by R
static void setup_sigpipe_handling(void) {
    if (!sigpipe_handled) {
        // Ignore SIGPIPE globally - we'll handle EPIPE errors locally
        // This prevents the process from being terminated by SIGPIPE
        // More precisely, we avoid R's default SIGPIPE handler
        signal(SIGPIPE, SIG_IGN);
        sigpipe_handled = 1;
        
        if (getenv("RBCFLIB_DEBUG") != NULL) {
            Rprintf("SIGPIPE handling set to SIG_IGN\n");
        }
    }
}

// Helper: Safe close with SIGPIPE protection
static void safe_close_fd(int fd) {
    if (fd >= 0 && fd != STDIN_FILENO && fd != STDOUT_FILENO && fd != STDERR_FILENO) {
        // Close might trigger SIGPIPE if there's buffered data and the other end is closed
        // But we've set SIGPIPE to SIG_IGN, so this should just return an error
        if (close(fd) == -1 && errno == EPIPE) {
            if (getenv("RBCFLIB_DEBUG") != NULL) {
                Rprintf("EPIPE error on close(fd=%d) - expected when other end closed early\n", fd);
            }
        }
    }
}

// Helper: convert SEXP of strings to char** (NULL-terminated)
static char **sexp_to_argv(SEXP args, SEXP command, int *argc_out) {
    int nargs = length(args);
    int argc = nargs + 2;
    char **argv = (char **)malloc((argc + 1) * sizeof(char *));
    if (argv == NULL) {
        error("Memory allocation failed");
    }
    
    // Get the bcftools binary path
    const char *bcftools_path = BCFToolsBinaryPath();
    if (bcftools_path == NULL) {
        free(argv);
        error("BCFTools binary path not found");
    }
    
    argv[0] = strdup(bcftools_path);
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
    
    // Setup SIGPIPE handling before any pipe operations
    setup_sigpipe_handling();
    
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
        const char * bcftools_plugins_path = BCFToolsPluginsPath();
        Rprintf("Using BCFTOOLS_PLUGINS: %s\n", bcftools_plugins_path);
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
            
            // Isolate child process to prevent R state corruption
            // Create new process group to isolate from parent
            setpgid(0, 0);
            
            // Signal protection - ensure child doesn't interfere with parent
            signal(SIGINT, SIG_DFL);
            signal(SIGTERM, SIG_DFL);
            signal(SIGPIPE, SIG_IGN);
            
            // Set up BCFTOOLS_PLUGINS environment variable
            // Get the bcftools binary path and derive plugins directory
            // FIX: ME actually get plugin from function or environment variable
            const char *bcftools_path = BCFToolsBinaryPath();
            const char *bcftools_plugins_path = BCFToolsPluginsPath();
            if (bcftools_plugins_path != NULL && strlen(bcftools_plugins_path) > 0) {
                setenv("BCFTOOLS_PLUGINS", bcftools_plugins_path, 1);
                if(getenv("RBCFLIB_DEBUG") != NULL) {
                    Rprintf("Set BCFTOOLS_PLUGINS to %s\n", bcftools_plugins_path);
                }
            }
            
            // Setup stdin (for all but first command)
            if (i > 0) {
                // Take input from previous pipe
                if (dup2(pipes[i-1][0], STDIN_FILENO) == -1) {
                    perror("dup2 stdin");
                    _exit(1);  // Use _exit instead of exit
                }
            }
            
            // Setup stdout
            if (i < num_commands - 1) {
                // Pipe output to next command
                if (dup2(pipes[i][1], STDOUT_FILENO) == -1) {
                    perror("dup2 stdout");
                    _exit(1);  // Use _exit instead of exit
                }
            } else {
                // Last command - write to file or stdout
                if (dup2(fd_stdout, STDOUT_FILENO) == -1) {
                    perror("dup2 stdout");
                    _exit(1);  // Use _exit instead of exit
                }
            }
            
            // Setup stderr (common to all commands)
            if (do_capture_stderr) {
                if (dup2(fd_stderr, STDERR_FILENO) == -1) {
                    perror("dup2 stderr");
                    _exit(1);  // Use _exit instead of exit
                }
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
            
            // Additional SIGPIPE protection in child process
            signal(SIGPIPE, SIG_IGN);
            
            // Run bcftools using execv
            execv(argv_values[i][0], argv_values[i]);
            
            // If execv returns, it means it failed
            perror("execv failed");
            
            // Free all argv arrays before exit
            for (int j = 0; j < num_commands; j++) {
                free_argv(argv_values[j]);
            }
            free(argv_values);
            // use _exit to avoid R finalization code
 
            _exit(1);
        }
    }
    
    // Parent: close all pipe ends
    for (int i = 0; i < num_commands - 1; i++) {
        safe_close_fd(pipes[i][0]);
        safe_close_fd(pipes[i][1]);
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

    if (fd_stdout != -1 && fd_stdout != STDOUT_FILENO) safe_close_fd(fd_stdout);
    if (fd_stderr != -1 && fd_stderr != STDERR_FILENO) safe_close_fd(fd_stderr);
    // Create result
    PROTECT(res = allocVector(INTSXP, num_commands));
    for (int i = 0; i < num_commands; i++) {
        INTEGER(res)[i] = statuses[i];
    }
    
    // Set attribute with combined command
    setAttrib(res, Rf_install("command"), cmd);
    
    UNPROTECT(2);
    return res;
}
#endif

#ifdef _WIN32
// Windows implementation using CreateProcess and pipes

// Helper: Safe close handle
static void safe_close_handle(HANDLE handle) {
    if (handle != NULL && handle != INVALID_HANDLE_VALUE) {
        CloseHandle(handle);
    }
}

// Helper: Create a pipe for Windows
static int create_pipe_win32(HANDLE *read_pipe, HANDLE *write_pipe) {
    SECURITY_ATTRIBUTES sa;
    sa.nLength = sizeof(SECURITY_ATTRIBUTES);
    sa.bInheritHandle = TRUE;
    sa.lpSecurityDescriptor = NULL;
    
    return CreatePipe(read_pipe, write_pipe, &sa, 0) ? 0 : -1;
}

// Helper: Convert R string vector to Windows command line
static char* build_command_line(SEXP command, SEXP args) {
    const char *bcftools_path = BCFToolsBinaryPath();
    if (bcftools_path == NULL) {
        return NULL;
    }
    
    int nargs = length(args);
    const char *cmd = CHAR(STRING_ELT(command, 0));
    
    // Calculate total length needed
    size_t total_len = strlen(bcftools_path) + strlen(cmd) + 10; // extra space
    for (int i = 0; i < nargs; i++) {
        total_len += strlen(CHAR(STRING_ELT(args, i))) + 3; // quotes + space
    }
    
    char *cmdline = (char*)malloc(total_len);
    if (cmdline == NULL) return NULL;
    
    // Build command line: "bcftools" "command" "arg1" "arg2" ...
    snprintf(cmdline, total_len, "\"%s\" \"%s\"", bcftools_path, cmd);
    
    for (int i = 0; i < nargs; i++) {
        strcat(cmdline, " \"");
        strcat(cmdline, CHAR(STRING_ELT(args, i)));
        strcat(cmdline, "\"");
    }
    
    return cmdline;
}

/**
 * Windows implementation of pipeline execution
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
    HANDLE pipes[num_commands - 1][2]; // pipes[i][0] = read, pipes[i][1] = write
    PROCESS_INFORMATION pis[num_commands];
    STARTUPINFO sis[num_commands];
    int statuses[num_commands];
    char **cmdlines = NULL;
    const char *std_out, *std_err;
    HANDLE h_stdout = INVALID_HANDLE_VALUE, h_stderr = INVALID_HANDLE_VALUE;
    SEXP res = R_NilValue, cmd = R_NilValue;
    int do_capture_stdout, do_capture_stderr;
    
    // Check for valid number of commands
    if (num_commands < 1) {
        error("At least one command is required");
    }
    
    // Extract capture flags
    do_capture_stdout = asLogical(capture_stdout);
    do_capture_stderr = asLogical(capture_stderr);
    
    // Setup stderr redirection if requested
    if (do_capture_stderr) {
        std_err = CHAR(STRING_ELT(stderr_file, 0));
        h_stderr = CreateFile(std_err, GENERIC_WRITE, 0, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
        if (h_stderr == INVALID_HANDLE_VALUE) {
            error("Could not open stderr file for writing");
        }
    } else {
        h_stderr = CreateFile("NUL", GENERIC_WRITE, 0, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
        if (h_stderr == INVALID_HANDLE_VALUE) {
            error("Could not open NUL for writing");
        }
    }
    
    // Setup stdout redirection if requested
    if (do_capture_stdout) {
        std_out = CHAR(STRING_ELT(stdout_file, 0));
        h_stdout = CreateFile(std_out, GENERIC_WRITE, 0, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
        if (h_stdout == INVALID_HANDLE_VALUE) {
            safe_close_handle(h_stderr);
            error("Could not open stdout file for writing: %s", std_out);
        }
    } else {
        h_stdout = CreateFile("NUL", GENERIC_WRITE, 0, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
        if (h_stdout == INVALID_HANDLE_VALUE) {
            safe_close_handle(h_stderr);
            error("Could not open NUL for writing");
        }
    }
    
    // Allocate memory for command lines
    cmdlines = (char**)malloc(num_commands * sizeof(char*));
    if (cmdlines == NULL) {
        safe_close_handle(h_stdout);
        safe_close_handle(h_stderr);
        error("Memory allocation failed for cmdlines");
    }
    
    // Convert R args to Windows command lines
    for (int i = 0; i < num_commands; i++) {
        SEXP cmd_i = VECTOR_ELT(commands, i);
        SEXP args_i = VECTOR_ELT(args, i);
        cmdlines[i] = build_command_line(cmd_i, args_i);
        if (cmdlines[i] == NULL) {
            for (int j = 0; j < i; j++) free(cmdlines[j]);
            free(cmdlines);
            safe_close_handle(h_stdout);
            safe_close_handle(h_stderr);
            error("Failed to build command line for command %d", i+1);
        }
    }
    
    // Debug output if requested
    if (getenv("RBCFLIB_DEBUG") != NULL) {
        Rprintf("Piping %d commands on Windows:\n", num_commands);
        for (int i = 0; i < num_commands; i++) {
            Rprintf("Command %d: %s\n", i+1, cmdlines[i]);
        }
    }
    
    // Create pipes between commands
    for (int i = 0; i < num_commands - 1; i++) {
        if (create_pipe_win32(&pipes[i][0], &pipes[i][1]) == -1) {
            // Clean up already created pipes
            for (int j = 0; j < i; j++) {
                safe_close_handle(pipes[j][0]);
                safe_close_handle(pipes[j][1]);
            }
            for (int j = 0; j < num_commands; j++) free(cmdlines[j]);
            free(cmdlines);
            safe_close_handle(h_stdout);
            safe_close_handle(h_stderr);
            error("CreatePipe failed");
        }
    }
    
    // Create child processes for each command
    for (int i = 0; i < num_commands; i++) {
        ZeroMemory(&sis[i], sizeof(STARTUPINFO));
        sis[i].cb = sizeof(STARTUPINFO);
        sis[i].dwFlags = STARTF_USESTDHANDLES;
        
        // Setup stdin
        if (i > 0) {
            sis[i].hStdInput = pipes[i-1][0];
        } else {
            sis[i].hStdInput = GetStdHandle(STD_INPUT_HANDLE);
        }
        
        // Setup stdout
        if (i < num_commands - 1) {
            sis[i].hStdOutput = pipes[i][1];
        } else {
            sis[i].hStdOutput = h_stdout;
        }
        
        // Setup stderr
        sis[i].hStdError = h_stderr;
        
        ZeroMemory(&pis[i], sizeof(PROCESS_INFORMATION));
        
        // Create the process
        if (!CreateProcess(NULL, cmdlines[i], NULL, NULL, TRUE, 0, NULL, NULL, &sis[i], &pis[i])) {
            // Error creating process, clean up
            for (int j = 0; j < num_commands - 1; j++) {
                safe_close_handle(pipes[j][0]);
                safe_close_handle(pipes[j][1]);
            }
            // Terminate already created processes
            for (int j = 0; j < i; j++) {
                TerminateProcess(pis[j].hProcess, 1);
                safe_close_handle(pis[j].hProcess);
                safe_close_handle(pis[j].hThread);
            }
            for (int j = 0; j < num_commands; j++) free(cmdlines[j]);
            free(cmdlines);
            safe_close_handle(h_stdout);
            safe_close_handle(h_stderr);
            error("CreateProcess failed for command %d", i+1);
        }
    }
    
    // Close all pipe handles in parent
    for (int i = 0; i < num_commands - 1; i++) {
        safe_close_handle(pipes[i][0]);
        safe_close_handle(pipes[i][1]);
    }
    
    // Wait for all processes to finish
    for (int i = 0; i < num_commands; i++) {
        WaitForSingleObject(pis[i].hProcess, INFINITE);
        DWORD exit_code;
        GetExitCodeProcess(pis[i].hProcess, &exit_code);
        statuses[i] = (int)exit_code;
        safe_close_handle(pis[i].hProcess);
        safe_close_handle(pis[i].hThread);
    }
    
    // Build the command attribute for the result
    int total_args = 0;
    for (int i = 0; i < num_commands; i++) {
        total_args++; // For command line string
        if (i < num_commands - 1) total_args++; // Add for pipe symbol
    }
    
    PROTECT(cmd = allocVector(STRSXP, total_args));
    
    int cmd_idx = 0;
    for (int i = 0; i < num_commands; i++) {
        SET_STRING_ELT(cmd, cmd_idx++, mkChar(cmdlines[i]));
        if (i < num_commands - 1) {
            SET_STRING_ELT(cmd, cmd_idx++, mkChar("|"));
        }
    }
    
    // Clean up resources
    for (int i = 0; i < num_commands; i++) {
        free(cmdlines[i]);
    }
    free(cmdlines);
    safe_close_handle(h_stdout);
    safe_close_handle(h_stderr);
    
    // Create result
    PROTECT(res = allocVector(INTSXP, num_commands));
    for (int i = 0; i < num_commands; i++) {
        INTEGER(res)[i] = statuses[i];
    }
    
    // Set attribute with combined command
    setAttrib(res, Rf_install("command"), cmd);
    
    UNPROTECT(2);
    return res;
}
#endif