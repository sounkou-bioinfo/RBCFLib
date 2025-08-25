#!/bin/bash
set -e

VCF="../../1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
VBI="../../1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.vbi"

# Index test
echo "Indexing $VCF -> $VBI"
time ./vbi_index index "$VCF" "$VBI"

# Query test (whole chr21)
echo "Querying chr21"
time ./vbi_index query --vcf "$VCF" --vbi "$VBI" 21 | head -20

echo "Querying chr21:10000000-10100000"
time ./vbi_index query --vcf "$VCF" --vbi "$VBI" 21:10000000-10100000 | head -20

echo "All tests completed."
