# RBCFLib Command Line Interfaces

This directory contains command-line interfaces (CLIs) for the RBCFLib package that allow users to access the bcftools functionality directly from the command line.

## Available CLI Tools

### 1. RBcftools

`RBcftools` is a direct command-line wrapper for bcftools commands using the RBCFLib package's `BCFToolsRun` function.

Usage:
```
RBcftools [command] [options]
```

Example:
```
RBcftools view -h sample.vcf.gz
```

### 2. BCFToolsCli

`BCFToolsCli` provides the same functionality as `RBcftools` but uses PascalCase command names for a more R-like style. The commands work exactly the same way, but are capitalized.

Usage:
```
BCFToolsCli <Command> [options]
```

Example:
```
BCFToolsCli View -h sample.vcf.gz
```

## Installation

When the RBCFLib package is installed, these scripts are placed in the package's installation directory. To make them accessible from your command line, you may need to:

1. Find where they are installed:
   ```r
   system.file("bin", package = "RBCFLib")
   ```

2. Add that directory to your PATH, or create symlinks to the scripts in a directory that's already in your PATH.

## Commands

Both CLIs support all standard bcftools commands:

- version/Version
- view/View
- index/Index
- query/Query
- call/Call
- mpileup/Mpileup
- concat/Concat
- merge/Merge
- norm/Norm
- stats/Stats
- annotate/Annotate
- cnv/Cnv
- consensus/Consensus
- convert/Convert
- csq/Csq
- filter/Filter
- gtcheck/Gtcheck
- plugin/Plugin
- roh/Roh
- isec/Isec
- reheader/Reheader
- sort/Sort
- head/Head
- help/Help

For more information on using a specific command, use:
```
RBcftools help [command]
```
or
```
BCFToolsCli Help [Command]
```
