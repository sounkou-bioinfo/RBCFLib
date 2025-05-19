# RBCFLib Command Line Interfaces a.k.a bcftools bundler in a R package :D

This directory contains command-line interfaces (CLIs) for the RBCFLib package that allow users to access the bcftools functionality directly from the command line but distributed in the R package. THIS DOES NOT SUPPORT PLUGINS except the ones provided as separed commnds.

## Available CLI Tools

### 1. DownloadGenomes

`DownloadGenomes` is a command-line interface for downloading human reference genomes using the RBCFLib package's `DownloadHumanReferenceGenomes` function.

Usage:
```
DownloadGenomes [options]
```

Options:
```
  --grch37-dir=DIR      Directory to store GRCh37 reference (default: ~/GRCh37)
  --grch38-dir=DIR      Directory to store GRCh38 reference (default: ~/GRCh38)
  --cytoband            Also download cytoband files
  --chain               Also download chain files for liftover
  --genomes             Download genome FASTA files (default behavior)
  --grch37-fasta=URL    URL for GRCh37 FASTA (optional)
  --grch38-fasta=URL    URL for GRCh38 FASTA (optional)
  --method=METHOD       Download method: 'wget', 'curl', 'auto', 'internal', 'libcurl', or 'lynx' (default: 'wget')
  --extra=ARG           Extra command-line arguments for 'wget' or 'curl' methods. For multiple arguments, use
                        comma-separated values, e.g., '--extra=-C,on' for wget's continue feature
  --help                Display this help message
```

Example:
```
# Download genomes with wget and enable continue feature for interrupted downloads
DownloadGenomes --method=wget --extra=-C,on

# Download genomes and also download cytoband and chain files
DownloadGenomes --cytoband --chain
```

### 2. RBcftools

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

`BCFToolsCli` provides the same functionality as `RBcftools` but uses PascalCase command names. The commands work exactly the same way, but are capitalized.

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
