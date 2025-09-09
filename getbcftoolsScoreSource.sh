#!/usr/bin/bash
set -euo pipefail
# download bcftools from github
version=1.22
THISDir=`dirname $0`
THISDir=`realpath ${THISDir}`
bcftoolsVersion="1.22"
# Download BCFToolsScore source as a tar.gz at the specified commit and unpack into plugins
commit=a9ffa435913101439974fce7c4812235b61e6df5
cd $(dirname $0)/src || exit 1
echo "$PWD"
echo "Preparing to download BCFToolsScore at commit ${commit}"
cd ${THISDir}/src/bcftools-${bcftoolsVersion}/plugins || exit 1
# Download the tar.gz archive from GitHub (adjust owner/repo if necessary)
wget -O "BCFToolsScore.tar.gz" \
    "https://github.com/freeseek/score/archive/${commit}.tar.gz"
tar -xzf "BCFToolsScore.tar.gz" -C "${PWD}" --strip-components=1
ls -lthr *
rm -f "BCFToolsScore.tar.gz"
