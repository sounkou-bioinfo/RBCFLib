#!/usr/bin/bash
set -euo pipefail
# download bcftools from github
version=1.22
THISDir=`dirname $0`
THISDir=`realpath ${THISDir}`
bcftoolsVersion="1.22"
# Download BCFToolsScore source as a tar.gz at the specified commit and unpack into plugins
commit=db4da22b6254636fe5e6865ea8ffe28133ca33eb
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
echo "Preparing to download BCFToolsScore at commit ${commit}"
cd ${THISDir}/src/bcftools-${bcftoolsVersion}/plugins || exit 1
# Download the tar.gz archive from GitHub (adjust owner/repo if necessary)
wget -O "BCFToolsScore.tar.gz" \
    "https://github.com/freeseek/score/archive/${commit}.tar.gz"
tar -xzf "BCFToolsScore.tar.gz" -C "${PWD}" --strip-components=1
rm -f "BCFToolsScore.tar.gz"
