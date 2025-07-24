# RBCFLIB

**RBCFLib** provides a minimamlist and GWAS-VCF centric R bindings to `htslib` and `bcftools`, enabling efficient manipulation of BCF/VCF genomic data files directly within R. It includes a wrapper for the `bcftools` command-line interface, with an implementation strategy adapted from  [pysam](https://github.com/pysam-developers/pysam)'s approach to wrapping samtools/bcftools command line libraries.

The goal that motivated **RBCFLib** is to provide a minimalist and performant R based GWAS-VCF based conversion and analysis tool as an alternative to existing sumstats manipulation tools, such as the very good [{MungeSumstats}](https://github.com/Al-Murphy/MungeSumstats), by leveraging the speed of [`bcftools munge`](https://github.com/freeseek/score). For more comprehensive and performant BCF/VCF I/O options in R, [{vcffppR}](https://github.com/Zilong-Li/vcfppR) or [{VariantAnnotation}](https://github.com/Bioconductor/VariantAnnotation).

**Currently in the box:**

*   **BCFTools Command-Line Wrapper:** Execute a subset of `bcftools` commands directly from R using `BCFToolsRun()`.
*   **Piping and Pipeline Support:** to be made more idiomatic R with NSE
    *   **`BCFToolsPipeline`**: Create a pipeline of multiple bcftools commands, executing them in sequence. **Note:** this should be valid bcftools pipeline
*   **Download Reference Genomes:** download human reference genomes and liftover files files with `DownloadHumanReferenceGenomes` .
*   **BCFTools Score Plugin wrappers:**
    *   **`BCFToolsMunge`**: Convert summary statistics to GWAS-VCF format
    *   **`BCFToolsScore`**: Compute polygenic scores using genotype data and weights 
    *   **`BCFToolsLiftover`**: Convert variants between genome assemblies
    *   **`BCFToolsMetal`**: Perform meta-analysis on GWAS-VCF files
    *   **`BCFToolsPGS`**: Compute polygenic score loadings (requires CHOLMOD, this is checked at install time)
    *   **`BCFToolsBLUP`** : compute blup 

 For more details on this methods, refer to the documentations of BCFTools Score https://github.com/freeseek/score

*   **Version Information:** Retrieve `htslib` and `bcftools` library versions using `HTSLibVersion()` and `BCFToolsVersion()`.
*   **Utility Functions:** Access the path to the bundled `bcftools` executable via `BCFToolsCLIPath()`
*   **Command-Line Interfaces:** Use the package's functions as command-line tools via the scripts in `inst/bin/`.




**Installation**

```r

# to install from r-universe

install.packages('RBCFLib', repos = c('https://sounkou-bioinfo.r-universe.dev'))

# install the developement version from github

devtools::install_github("sounkou-bioinfo/RBCFLib@develop")

```

**Example Usage**

```r
# Load the library
library(RBCFLib)

# Get bcftools version
BCFToolsVersion()

# Run a bcftools command (example: view a VCF file)
vcfFile <- system.file("exdata", "imputed.gt.vcf.gz", package = "RBCFLib")
BCFToolsRun("view", c("-i" , ' ALT == "A" ' , vcfFile))

# Create a pipeline of multiple bcftools commands
# Note: The input file can be included in the first command's arguments
result <- BCFToolsPipeline(
  "view", c("-r", "chr1:1000-2000", vcfFile), 
  "view", c("-i", "QUAL>20"),
  "view", c("-H")
)
print(result$stdout)
```

**TODO:**

The following features are planned for future development:

*   **Enhanced `bcftools` Wrapping:**
    *   [ ] Improve handling of stdout/stderr for interactive R sessions. Probably change the BcftoolsRun function signature
    *   [ ] Implement robust user interrupt handling (Ctrl+C)  This might involve running `bcftools` commands in a separate thread or background R process.
    *   [ ] more plugins for basic statistics

*   **Direct BCF/VCF Data Manipulation in R:**
    *   [ ] Develop functions for reading/scanning BCF/VCF files into R data structures ? (e.g., data frames or similar).
    *   [ ] Minimal streaming like rbcf and vcfppR

*   **Some Interval Operations (kfunctions):**
    *   [ ] Provide some utils for working with genomic intervals.

*   **Tabix Support:**
    *   [ ] A tabix wrapper
