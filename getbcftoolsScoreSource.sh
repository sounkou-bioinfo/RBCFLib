#!/usr/bin/bash
set -eufo pipefail
THISDir=`dirname $0`
THISDir=`realpath ${THISDir}`
bcftoolsVersion="1.22"
# Download BCFToolsScore source as a tar.gz at the specified commit and unpack into src/BCFToolsScore
commit=db4da22b6254636fe5e6865ea8ffe28133ca33eb

destDir="$(dirname "$0")/src/BCFToolsScore"
echo "Preparing to download BCFToolsScore at commit ${commit}"
cd ${THISDir}/src/bcftools-${bcftoolsVersion}/plugins || exit 1
mkdir -p "${destDir}" || exit 1
# Download the tar.gz archive from GitHub (adjust owner/repo if necessary)
wget -O "BCFToolsScore.tar.gz" \
    "https://github.com/freeseek/score/archive/${commit}.tar.gz"
tar -xzf "BCFToolsScore.tar.gz" -C "${destDir}" --strip-components=1
rm -f "BCFToolsScore.tar.gz"
exit 0
