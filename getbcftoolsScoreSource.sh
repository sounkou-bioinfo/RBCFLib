#!/usr/bin/bash
set -eufo pipefail
# Download BCFToolsScore source as a tar.gz at the specified commit and unpack into src/BCFToolsScore
commit=db4da22b6254636fe5e6865ea8ffe28133ca33eb

destDir="$(dirname "$0")/src/BCFToolsScore"
echo "Preparing to download BCFToolsScore at commit ${commit}"

mkdir -p "${destDir}" || exit 1
# Download the tar.gz archive from GitHub (adjust owner/repo if necessary)
wget -O "${destDir}/BCFToolsScore.tar.gz" \
    "https://github.com/freeseek/score/archive/${commit}.tar.gz"

echo "Unpacking BCFToolsScore into ${destDir}"
tar -xzf "${destDir}/BCFToolsScore.tar.gz" -C "${destDir}" --strip-components=1
rm -f "${destDir}/BCFToolsScore.tar.gz"
echo "BCFToolsScore source is ready in ${destDir}"
