#!/usr/bin/bash
set -eufo pipefail
# download bcftools from github
version=1.21
cd $(dirname $0)/src || exit 1
echo "Downloading bcftools version ${version} from GitHub"
wget https://github.com/samtools/bcftools/releases/download/${version}/bcftools-${version}.tar.bz2
echo "Unpacking bcftools version ${version}"
tar -xjf bcftools-${version}.tar.bz2
echo "Cleaning up"
rm -f bcftools-${version}.tar.bz2 || true
cd bcftools-${version} || exit 1
[[ -d htslib-${version} ]] && rm -rf htslib-${version} || true
mv -f htslib-{version} ../ || exit 1
