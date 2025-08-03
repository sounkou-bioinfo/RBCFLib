#!/usr/bin/bash
set -euo pipefail
# download bcftools from github
version=1.22
cd $(dirname $0)/src || exit 1
echo "$PWD"
ls -d ${PWD}/bcftools* | xargs -I {} -t rm -rf {} || true 
ls -d ${PWD}/htslib-* | xargs -I {} -t rm -rf {} || true
echo "Downloading bcftools version ${version} from GitHub"
wget https://github.com/samtools/bcftools/releases/download/${version}/bcftools-${version}.tar.bz2
echo "Unpacking bcftools version ${version}"
tar -xjf bcftools-${version}.tar.bz2
echo "Cleaning up"
rm -f bcftools-${version}.tar.bz2 || true
cd bcftools-${version} || exit 1
mv -f htslib-${version} ../ || exit 1
