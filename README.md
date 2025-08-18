# RBCFLIB

**RBCFLib** provides a minimalist and GWAS-VCF centric R wrapper for `htslib` and `bcftools`, enabling efficient manipulation of BCF/VCF genomic data files directly within R. The package includes bundled bcftools source code that is compiled during installation, then executes bcftools commands as optimized subprocesses through a C interface that handles piping, error capture, and process management.

The goal that motivated **RBCFLib** is to provide a minimalist and performant R-based GWAS-VCF conversion and analysis tool as an alternative to existing sumstats manipulation tools, such as the very good [{MungeSumstats}](https://github.com/Al-Murphy/MungeSumstats), by leveraging the speed of [`bcftools munge`](https://github.com/freeseek/score). Unlike library-based approaches, RBCFLib executes the full bcftools CLI as optimized subprocesses, ensuring compatibility with all bcftools features and plugins while maintaining high performance through efficient process management and piping. For more comprehensive and performant BCF/VCF I/O options in R, consider [{vcfppR}](https://github.com/Zilong-Li/vcfppR) or [{VariantAnnotation}](https://github.com/Bioconductor/VariantAnnotation).

**Currently in the box:**

-   **BCFTools Executable Wrapper:** Execute any `bcftools` command directly from R using `BCFToolsRun()` - the package bundles the complete bcftools executable and runs it as optimized subprocesses with full error handling and output capture.
-   **Piping and Pipeline Support:** Create complex data processing workflows with idiomatic R syntax
    -   **`BCFToolsPipeline`**: Chain multiple bcftools commands in efficient pipelines, with automatic process management and piping handled in optimized C code for maximum performance
-   **Download Reference Genomes:** download human reference genomes and liftover files files with `DownloadHumanReferenceGenomes` .
-   **BCFTools Score Plugin wrappers:**
    -   **`BCFToolsMunge`**: Convert summary statistics to GWAS-VCF format
    -   **`BCFToolsScore`**: Compute polygenic scores using genotype data and weights
    -   **`BCFToolsLiftover`**: Convert variants between genome assemblies
    -   **`BCFToolsMetal`**: Perform meta-analysis on GWAS-VCF files
    -   **`BCFToolsPGS`**: Compute polygenic score loadings (requires CHOLMOD, this is checked at install time)
    -   **`BCFToolsBLUP`** : compute blup

For more details on these methods, refer to the BCFTools Score documentation: https://github.com/freeseek/score

-   **Version Information:** Retrieve `htslib` and `bcftools` library versions using `HTSLibVersion()` and `BCFToolsVersion()`.
-   **Utility Functions:** Access the path to the bundled `bcftools` executable via `BCFToolsBinaryPath()`
-   **Command-Line Interfaces:** Use the package's functions as command-line tools via the scripts in `inst/bin/`.

**Installation**

The package includes bundled bcftools source code and automatically compiles the complete bcftools suite during installation, ensuring compatibility and performance across different systems without requiring external dependencies.

``` r

# to install from r-universe

install.packages('RBCFLib', repos = c('https://sounkou-bioinfo.r-universe.dev'))

# install the developement version from github

devtools::install_github("sounkou-bioinfo/RBCFLib@develop")
```

**Example Usage**

``` r
# Load the library
library(RBCFLib)

# Get bcftools version (from the bundled executable)
BCFToolsVersion()

# Run a bcftools command - executes as subprocess with full error handling
vcfFile <- system.file("exdata", "imputed.gt.vcf.gz", package = "RBCFLib")
results <- BCFToolsRun("view", c("-i" , ' ALT == "A" ' , vcfFile))
print(results$stdout)

# Create a pipeline of multiple bcftools commands - efficient C-managed piping
# The commands are executed as connected subprocesses with automatic pipe management
results <- BCFToolsPipeline(
    "view", c("-r", "chr21", "-Ob", vcfFile),
    "view", c("--no-version"),
    "query", c( "-f", "%CHROM\\t%POS\\t%REF\\t%ALT\\n")
)
print(results$stdout)
```

**TODO:**

The following features are planned for future development:

-   **Enhanced `bcftools` Wrapping:**
    -   [ ] more plugins for basic statistics
-   **Direct BCF/VCF Data Manipulation in R:**
    -   [ ] Develop functions for reading/scanning BCF/VCF files into R data structures ? (e.g., data frames or similar).
    -   [ ] Minimal streaming like rbcf and vcfppR
-   **Some Interval Operations (kfunctions):**
    -   [ ] Provide some utils for working with genomic intervals.
-   **Tabix Support:**
    -   [ ] A tabix wrapper