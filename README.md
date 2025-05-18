# RBCLIB

**RBCFLib** provides R bindings to `htslib` and `bcftools`, enabling efficient manipulation of BCF/VCF genomic data files directly within R. It includes a wrapper for the `bcftools` command-line interface, with an implementation strategy inspired by [pysam](https://github.com/pysam-developers/pysam)'s approach to wrapping samtools/bcftools command line libraries tools.

A primary goal that motivated **RBCFLib** is to performance R based GWAS-VCF centric alternative to existing sumstats manipulation tools, such as [{MungeSumstats}](https://github.com/Al-Murphy/MungeSumstats), by leveraging the speed of `bcftools munge`.

**Current Features:**

*   **BCFTools Command-Line Wrapper:** Execute `bcftools` commands directly from R using `BCFToolsRun()`.
*   **Version Information:** Retrieve `htslib` and `bcftools` library versions using `HTSLibVersion()` and `BCFToolsVersion()`.
*   **Utility Functions:** Access the path to the bundled `bcftools` executable via `BCFToolsCLIPath()`.

**TODO:**

The following features are planned for future development:

*   **Enhanced `bcftools` Wrapping:**
    *   [ ] Improve handling of stdout/stderr for interactive R sessions.
    *   [ ] Implement robust user interrupt handling (Ctrl+C). This might involve running `bcftools` commands in a separate thread or background R process.

*   **Direct BCF/VCF Data Manipulation in R:**
    *   [ ] Develop functions for reading/scanning BCF/VCF files into R data structures (e.g., data frames or similar).
    *   [ ] Explore streaming capabilities for large VCF/BCF files to manage memory efficiently.

*   **`faidx` Integration:**
    *   [ ] Implement functions for fast FASTA file indexing and sequence retrieval.

*   **Interval Operations (kfunctions):**
    *   [ ] Provide tools for working with genomic intervals.

*   **Tabix Support:**
    *   [ ] Enable importing data from Tabix-indexed tabular files into BCF/VCF format or R objects.

**Installation**

```r
# Instructions for installing from GitHub (once available)
devtools::install_github("sounkou-bioinfo/RBCFLib")
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
# Use bcftools munge help
tsvFile <- system.file("exdata", "test_plink.tsv", package = "RBCFLib")
refFile <- system.file("exdata", "Test.fa", package = "RBCFLib")
outFile <- tempfile(fileext = ".vcf")

out <- BCFToolsMunge(
    input_file = tsvFile,
    fasta_ref = refFile,
    output = outFile,
    columns = "PLINK" # Adding required columns parameter
)
```
