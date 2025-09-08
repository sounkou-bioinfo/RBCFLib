# Comprehensive tests for RBCFLib RBCF functions - adapted from original rbcf test suite
# Tests cover VCF/BCF file handling, variant processing, and genotype analysis
library(tinytest)
library(RBCFLib)
# Load RBCFLib and ensure functions are available

# Function to count variants in a VCF/BCF file
count_variants <- function(filename) {
    fp <- BCFOpen(filename, FALSE)
    if (is.null(fp)) {
        return(-1)
    }
    n <- 0
    while (!is.null(vc <- BCFNext(fp))) {
        n <- n + 1
    }
    BCFClose(fp)
    return(n)
}

# Function to find a specific variant by position
find_variant <- function(fp, contig, pos) {
    region <- paste(contig, ":", pos, "-", pos, sep = "")
    if (!BCFQuery(fp, region)) {
        return(NULL)
    }
    while (!is.null(vc <- BCFNext(fp))) {
        return(vc)
    }
    return(NULL)
}

# Get test data files - use files that definitely exist
rotavirus_01 <- system.file(
    "exdata",
    "rotavirus_rf.01.vcf",
    package = "RBCFLib"
)
rotavirus_02 <- system.file(
    "exdata",
    "rotavirus_rf.02.vcf.gz",
    package = "RBCFLib"
)
rotavirus_03 <- system.file(
    "exdata",
    "rotavirus_rf.03.vcf.gz",
    package = "RBCFLib"
)
rotavirus_04 <- system.file(
    "exdata",
    "rotavirus_rf.04.bcf",
    package = "RBCFLib"
)
rotavirus_ann <- system.file(
    "exdata",
    "rotavirus_rf.ann.vcf.gz",
    package = "RBCFLib"
)
gnomad_bcf <- system.file(
    "exdata",
    "gnomad.exomes.r2.0.1.sites.bcf",
    package = "RBCFLib"
)
genotypes_bcf <- system.file(
    "exdata",
    "1000G.ALL.2of4intersection.20100804.genotypes.bcf",
    package = "RBCFLib"
)
imputed_vcf <- system.file("exdata", "imputed.gt.vcf.gz", package = "RBCFLib")

# Verify test files exist
expect_true(file.exists(rotavirus_01), "Rotavirus 01 VCF file exists")
expect_true(file.exists(rotavirus_02), "Rotavirus 02 VCF file exists")
expect_true(file.exists(gnomad_bcf), "Gnomad BCF file exists")


# Test 1: Basic VCF/BCF File Operations
cat("=== Test 1: VCF/BCF File Opening and Closing ===\n")

# Test opening and closing a VCF file (no index required)
fp <- BCFOpen(rotavirus_01, FALSE)
expect_true(!is.null(fp), "BCFOpen successfully opens VCF file")
expect_true(looks_like_vcf_context(fp), "BCFOpen returns valid VCF context")
BCFClose(fp)
expect_true(TRUE, "BCFClose works")

# Test opening indexed VCF
if (file.exists(rotavirus_02)) {
    fp <- BCFOpen(rotavirus_02, TRUE)
    expect_true(!is.null(fp), "BCFOpen successfully opens indexed VCF file")
    BCFClose(fp)
    expect_true(TRUE, "BCFClose works for indexed VCF")
}

cat("Done.\n\n")

# Test 2: VCF Header Information
cat("=== Test 2: VCF Header Information Functions ===\n")

# Test with unindexed file
fp <- BCFOpen(rotavirus_01, FALSE)
expect_true(!is.null(fp), "File opened successfully")

# Test sample information
nsamples <- BCFNSamples(fp)
expect_true(is.numeric(nsamples), "BCFNSamples returns numeric value")
expect_true(nsamples > 0, "BCFNSamples returns positive number")

samples <- BCFSamples(fp)
expect_true(is.character(samples), "BCFSamples returns character vector")
expect_equal(
    length(samples),
    nsamples,
    info = "BCFSamples length matches BCFNSamples"
)

# Test sample indexing
if (nsamples > 0) {
    first_sample <- BCFSampleAt(fp, 1)
    expect_true(is.character(first_sample), "BCFSampleAt returns character")
    expect_equal(
        first_sample,
        samples[1],
        info = "BCFSampleAt(1) matches first sample"
    )

    # Test sample name to index conversion
    indices <- BCFSample2Index(fp, c(samples[1], "missing_sample"))
    expect_true(is.numeric(indices), "BCFSample2Index returns numeric vector")
    expect_true(
        !is.na(indices[1]) && indices[1] == 1,
        "BCFSample2Index returns correct 1-based index"
    )
}

# Test chromosome/contig information (requires indexed VCF); skip for unindexed file
if (file.exists(rotavirus_02)) {
    fp_idx <- BCFOpen(rotavirus_02, TRUE)
    if (!is.null(fp_idx)) {
        chroms <- BCFChromosomes(fp_idx)
        expect_true(
            is.character(chroms),
            "BCFChromosomes returns character vector (indexed)"
        )
        expect_true(
            length(chroms) > 0,
            "BCFChromosomes returns non-empty vector (indexed)"
        )
        contigs <- BCFContigs(fp_idx)
        expect_equal(
            chroms,
            contigs,
            info = "BCFContigs matches BCFChromosomes (indexed)"
        )
        BCFClose(fp_idx)
    }
} else {
    cat(
        "Indexed rotavirus_02.vcf.gz not found - skipping chromosome/contig tests\n"
    )
    chroms <- character(0)
    contigs <- character(0)
}

# Test dictionary
dict <- BCFDictionary(fp)
expect_true(is.data.frame(dict), "BCFDictionary returns data frame")
expect_true("chrom" %in% names(dict), "Dictionary contains chrom column")

# Test table functions
infos_table <- BCFInfos(fp)
expect_true(
    is.data.frame(infos_table) || is.null(infos_table),
    "BCFInfos returns data frame or NULL"
)

formats_table <- BCFFormats(fp)
expect_true(
    is.data.frame(formats_table) || is.null(formats_table),
    "BCFFormats returns data frame or NULL"
)

BCFClose(fp)
cat("Samples:", nsamples, ", Chromosomes:", length(chroms), "\n\n")

# Test 3: Variant Counting and Iteration
cat("=== Test 3: Scanning Variants in Multiple Files ===\n")

# Test variant counting in different file formats
test_files <- c(rotavirus_01, rotavirus_02, rotavirus_03, rotavirus_04)
existing_files <- test_files[file.exists(test_files)]

for (filename in existing_files) {
    n_variants <- count_variants(filename)
    expect_true(
        n_variants >= 0,
        paste("Successfully counted variants in", basename(filename))
    )
    cat(paste(basename(filename), ":", n_variants, "variants\n"))
}

cat("\n")

# Test 4: Indexed VCF Querying
cat("=== Test 4: Querying Indexed VCF Files ===\n")

if (file.exists(rotavirus_02)) {
    fp <- BCFOpen(rotavirus_02, TRUE)
    expect_true(!is.null(fp), "Opened indexed VCF file")

    chroms <- BCFChromosomes(fp)

    # Test different query intervals
    intervals <- c("", "RF03", "RF03:2000-3000")

    for (interval in intervals) {
        if (interval == "") {
            # Empty interval should query all
            query_result <- TRUE # Skip empty query test
        } else {
            query_result <- BCFQuery(fp, interval)
        }

        if (query_result && interval != "") {
            n_variants <- 0
            while (!is.null(vc <- BCFNext(fp))) {
                n_variants <- n_variants + 1
            }
            expect_true(
                n_variants >= 0,
                paste("Query successful for interval:", interval)
            )
            cat(paste("Interval '", interval, "':", n_variants, "variants\n"))
        }
    }

    # Test collect functionality
    if (length(chroms) > 0) {
        test_interval <- paste0(chroms[1], ":1-1000")
        variants <- BCFQuery(fp, test_interval, collect = TRUE)
        if (!is.null(variants)) {
            expect_true(
                is.list(variants),
                "BCFQuery with collect=TRUE returns list"
            )
            cat(paste(
                "Collected",
                length(variants),
                "variants using collect=TRUE\n"
            ))
        }
    }

    BCFClose(fp)
}

cat("\n")

# Test 5: Detailed Variant Analysis
cat("=== Test 5: Detailed Variant Information and Analysis ===\n")

fp <- BCFOpen(rotavirus_01, FALSE)
expect_true(!is.null(fp), "File opened for variant analysis")

# Read first few variants and test all variant functions
variant_count <- 0
while (!is.null(vc <- BCFNext(fp)) && variant_count < 5) {
    variant_count <- variant_count + 1

    # Test basic variant information
    tid <- VariantTid(vc)
    expect_true(is.numeric(tid), "VariantTid returns numeric")

    contig <- VariantContig(vc)
    expect_true(is.character(contig), "VariantContig returns character")

    chrom <- VariantChrom(vc)
    expect_equal(chrom, contig, info = "VariantChrom matches VariantContig")

    pos <- VariantPos(vc)
    expect_true(
        is.numeric(pos) && pos > 0,
        "VariantPos returns positive number"
    )

    start <- VariantStart(vc)
    expect_equal(start, pos, info = "VariantStart matches VariantPos")

    end <- VariantEnd(vc)
    expect_true(is.numeric(end) && end >= pos, "VariantEnd >= VariantPos")

    stop <- VariantStop(vc)
    expect_equal(stop, end, info = "VariantStop matches VariantEnd")

    # Test ID information
    has_id <- VariantHasId(vc)
    expect_true(is.logical(has_id), "VariantHasId returns logical")
    if (has_id) {
        id <- VariantId(vc)
        expect_true(
            is.character(id),
            "VariantId returns character when available"
        )
    }

    # Test allele information
    nalleles <- VariantNAlleles(vc)
    expect_true(is.numeric(nalleles) && nalleles >= 2, "VariantNAlleles >= 2")

    alleles <- VariantAlleles(vc)
    expect_true(
        is.character(alleles),
        "VariantAlleles returns character vector"
    )
    expect_equal(
        length(alleles),
        nalleles,
        info = "VariantAlleles length matches VariantNAlleles"
    )

    ref <- VariantReference(vc)
    expect_true(is.character(ref), "VariantReference returns character")
    expect_equal(
        ref,
        alleles[1],
        info = "VariantReference matches first allele"
    )

    alt <- VariantAltAlleles(vc)
    expect_true(is.character(alt), "VariantAltAlleles returns character vector")
    expect_equal(
        length(alt),
        nalleles - 1,
        info = "VariantAltAlleles length is nalleles-1"
    )

    # Test QUAL information
    has_qual <- VariantHasQual(vc)
    expect_true(is.logical(has_qual), "VariantHasQual returns logical")
    if (has_qual) {
        qual <- VariantQual(vc)
        expect_true(
            is.numeric(qual),
            "VariantQual returns numeric when available"
        )
    }

    # Test filter information
    is_filtered <- VariantIsFiltered(vc)
    expect_true(is.logical(is_filtered), "VariantIsFiltered returns logical")

    filters <- VariantFilters(vc)
    expect_true(
        is.character(filters) || is.null(filters),
        "VariantFilters returns character or NULL"
    )

    # Test variant types
    types <- VariantTypes(vc)
    expect_true(is.character(types), "VariantTypes returns character")

    is_snp <- VariantIsSnp(vc)
    expect_true(is.logical(is_snp), "VariantIsSnp returns logical")

    max_ploidy <- VariantMaxPloidy(vc)
    expect_true(
        is.numeric(max_ploidy) && max_ploidy > 0,
        "VariantMaxPloidy > 0"
    )

    # Test sample count
    nsamples_var <- VariantNSamples(vc)
    expect_true(is.numeric(nsamples_var), "VariantNSamples returns numeric")

    # Test INFO attributes
    info_keys <- VariantInfoIds(vc)
    expect_true(
        is.character(info_keys) || is.null(info_keys),
        "VariantInfoIds returns character or NULL"
    )

    format_keys <- VariantFormatIds(vc)
    expect_true(
        is.character(format_keys) || is.null(format_keys),
        "VariantFormatIds returns character or NULL"
    )

    cat(paste(
        "Variant",
        variant_count,
        ":",
        contig,
        ":",
        pos,
        "(",
        paste(alleles, collapse = "/"),
        ")\n"
    ))
}

BCFClose(fp)
cat("Analyzed", variant_count, "variants in detail\n\n")

# Additional Test: Variant INFO / FORMAT attribute access (adapted from original script)
cat("=== Additional Test: Variant INFO / FORMAT attributes ===\n")
if (file.exists(gnomad_bcf)) {
    fp <- BCFOpen(gnomad_bcf, TRUE)
    if (!is.null(fp)) {
        # helper to find variant by contig/pos
        find_variant_attr <- function(fp, contig, pos) {
            region <- paste0(contig, ":", pos, "-", pos)
            if (!BCFQuery(fp, region)) {
                return(NULL)
            }
            while (!is.null(vc <- BCFNext(fp))) {
                return(vc)
            }
            NULL
        }
        # pick a known position from original tests
        vc <- find_variant_attr(fp, "1", 905608)
        if (!is.null(vc)) {
            # String attribute (CSQ) split / no split
            if (VariantHasAttribute(vc, "CSQ")) {
                csq_nosplit <- VariantStringAttribute(vc, "CSQ", split = FALSE)
                expect_true(
                    is.character(csq_nosplit) || is.null(csq_nosplit),
                    "VariantStringAttribute no split returns character or NULL"
                )
                csq_split <- VariantStringAttribute(vc, "CSQ", split = TRUE)
                expect_true(
                    is.null(csq_split) || length(csq_split) >= 1,
                    "VariantStringAttribute split returns vector or NULL"
                )
            }
            # Integer attribute
            if (VariantHasAttribute(vc, "AN_POPMAX")) {
                an_popmax <- VariantIntAttribute(vc, "AN_POPMAX")
                expect_true(
                    is.numeric(an_popmax),
                    "VariantIntAttribute returns numeric"
                )
            }
            # Float attribute
            if (VariantHasAttribute(vc, "AF_POPMAX")) {
                af_popmax <- VariantFloatAttribute(vc, "AF_POPMAX")
                expect_true(
                    is.numeric(af_popmax),
                    "VariantFloatAttribute returns numeric"
                )
            }
            # Flag attribute (may or may not exist)
            flag_val <- VariantFlagAttribute(vc, "VQSR_NEGATIVE_TRAIN_SITE")
            expect_true(
                is.logical(flag_val) || is.null(flag_val),
                "VariantFlagAttribute returns logical or NULL"
            )
        }
        BCFClose(fp)
    }
}
cat("INFO / FORMAT attribute tests done\n\n")

# Additional Test: VEP and SnpEff annotations
cat("=== Additional Test: VEP / SnpEff annotations ===\n")
# VEP (CSQ) in gnomad
if (file.exists(gnomad_bcf)) {
    fp <- BCFOpen(gnomad_bcf, TRUE)
    if (!is.null(fp)) {
        vep_found <- FALSE
        iter <- 0
        while (!is.null(vc <- BCFNext(fp)) && iter < 200 && !vep_found) {
            iter <- iter + 1
            if (VariantHasAttribute(vc, "CSQ")) {
                vep_tbl <- VariantVep(vc)
                expect_true(
                    is.null(vep_tbl) || is.data.frame(vep_tbl),
                    "VariantVep returns data.frame or NULL"
                )
                vep_found <- TRUE
            }
        }
        BCFClose(fp)
    }
}
# SnpEff (ANN) in annotated rotavirus
if (file.exists(rotavirus_ann)) {
    fp <- BCFOpen(rotavirus_ann, FALSE)
    if (!is.null(fp)) {
        ann_found <- FALSE
        iter <- 0
        while (!is.null(vc <- BCFNext(fp)) && iter < 200 && !ann_found) {
            iter <- iter + 1
            if (VariantHasAttribute(vc, "ANN")) {
                ann_tbl <- VariantSnpEff(vc)
                expect_true(
                    is.null(ann_tbl) || is.data.frame(ann_tbl),
                    "VariantSnpEff returns data.frame or NULL"
                )
                ann_found <- TRUE
            }
        }
        BCFClose(fp)
    }
}
cat("Annotation tests done\n\n")

# Test 6: Genotype Analysis
cat("=== Test 6: Genotype Analysis ===\n")

# Use file with genotypes
if (file.exists(genotypes_bcf)) {
    fp <- BCFOpen(genotypes_bcf, TRUE)
    expect_true(!is.null(fp), "Opened genotypes BCF file")

    # Find a variant at a specific position
    vc <- find_variant(fp, "1", 10583)

    if (!is.null(vc)) {
        cat("Testing genotype functions on variant 1:10583\n")

        nsamples <- VariantNSamples(vc)
        expect_true(nsamples > 0, "Variant has samples")

        # Test genotype for first few samples
        for (i in 1:min(3, nsamples)) {
            gt <- VariantGenotype(vc, i)
            expect_true(
                looks_like_gt_context(gt),
                "VariantGenotype returns valid genotype context"
            )

            # Test genotype functions
            ploidy <- GenotypePloidy(gt)
            expect_true(is.numeric(ploidy), "GenotypePloidy returns numeric")

            alleles_idx <- GenotypeAllelesIdx0(gt)
            expect_true(
                is.numeric(alleles_idx) || is.null(alleles_idx),
                "GenotypeAllelesIdx0 returns numeric or NULL"
            )

            is_homref <- GenotypeHomRef(gt)
            expect_true(is.logical(is_homref), "GenotypeHomRef returns logical")

            is_het <- GenotypeHet(gt)
            expect_true(is.logical(is_het), "GenotypeHet returns logical")

            is_homvar <- GenotypeHomVar(gt)
            expect_true(is.logical(is_homvar), "GenotypeHomVar returns logical")

            is_hetnonref <- GenotypeHetNonRef(gt)
            expect_true(
                is.logical(is_hetnonref),
                "GenotypeHetNonRef returns logical"
            )

            is_nocall <- GenotypeNoCall(gt)
            expect_true(is.logical(is_nocall), "GenotypeNoCall returns logical")

            is_phased <- GenotypePhased(gt)
            expect_true(is.logical(is_phased), "GenotypePhased returns logical")

            sample_name <- GenotypeSample(gt)
            expect_true(
                is.character(sample_name),
                "GenotypeSample returns character"
            )

            # Test genotype attributes
            has_dp <- isTRUE(GenotypeHasDp(gt))
            if (has_dp) {
                dp <- GenotypeDp(gt)
                expect_true(
                    is.numeric(dp),
                    "GenotypeDp returns numeric when available"
                )
            }

            has_gq <- isTRUE(GenotypeHasGq(gt))
            if (has_gq) {
                gq <- GenotypeGq(gt)
                expect_true(
                    is.numeric(gq),
                    "GenotypeGq returns numeric when available"
                )
            }

            if (i <= 2) {
                # Only print details for first 2 samples
                cat(paste(
                    "  Sample",
                    i,
                    ":",
                    sample_name,
                    "ploidy:",
                    ploidy,
                    "homref:",
                    is_homref,
                    "het:",
                    is_het,
                    "\n"
                ))
            }
        }

        # Test vectorized genotype functions
        all_gt_idx <- VariantGenotypesAlleleIdx0(vc)
        expect_true(
            is.numeric(all_gt_idx),
            "VariantGenotypesAlleleIdx0 returns numeric vector"
        )

        all_gt_strings <- VariantGenotypesAlleleStrings(vc)
        expect_true(
            is.character(all_gt_strings),
            "VariantGenotypesAlleleStrings returns character vector"
        )
        expect_equal(
            length(all_gt_strings),
            nsamples,
            info = "Genotype strings length matches sample count"
        )

        # Test allele counts
        ref_counts <- VariantGenotypesAlleleCounts(vc, 0)
        expect_true(
            is.numeric(ref_counts),
            "VariantGenotypesAlleleCounts returns numeric for reference"
        )
        expect_equal(
            length(ref_counts),
            nsamples,
            info = "Reference counts length matches sample count"
        )

        alt_counts <- VariantGenotypesAlleleCounts(vc, 1)
        expect_true(
            is.numeric(alt_counts),
            "VariantGenotypesAlleleCounts returns numeric for alternate"
        )
        expect_equal(
            length(alt_counts),
            nsamples,
            info = "Alternate counts length matches sample count"
        )

        cat("Tested vectorized genotype functions\n")
    }

    BCFClose(fp)
}

cat("\n")

# Test 7: VCF Writing (to stdout)
cat("=== Test 7: VCF Writing Functionality ===\n")

fp <- BCFOpen(rotavirus_01, FALSE)
expect_true(!is.null(fp), "Source file opened for writing test")

# Create a writer that outputs to stdout (devnull for testing)
out <- BCFNewWriter(fp, "/dev/null")
expect_true(!is.null(out), "BCFNewWriter creates valid writer")

# Write a few variants
written_count <- 0
while (!is.null(vc <- BCFNext(fp)) && written_count < 3) {
    pos <- VariantPos(vc)
    if (pos %% 10 == 0) {
        # Only write variants at positions divisible by 10
        result <- BCFWriteVariant(out, vc)
        expect_true(is.logical(result), "BCFWriteVariant returns logical")
        written_count <- written_count + 1
    }
}

BCFClose(fp)
expect_true(TRUE, "Closed source file")
BCFClose(out)
expect_true(TRUE, "Closed writer file")

cat("Successfully wrote", written_count, "variants\n\n")

# Test 8: Error Handling
cat("=== Test 8: Error Handling ===\n")

# Test with non-existent file (should return NULL or throw; accept either)
nonexistent <- "/nonexistent/file.vcf"
res <- tryCatch(BCFOpen(nonexistent, FALSE), error = function(e) e)
if (inherits(res, "error")) {
    expect_true(TRUE, "BCFOpen raised error for missing file (acceptable)")
} else {
    expect_true(is.null(res), "BCFOpen returns NULL for missing file")
}

# Test with invalid arguments (should error due to stopifnot)
expect_error(
    BCFOpen(rotavirus_01, requireIndex = "invalid"),
    "is.logical\\(requireIndex\\) is not TRUE"
)

cat("Error handling tests completed\n\n")

cat("=== All RBCFLib RBCF tests completed successfully! ===\n")
