#!/bin/env bash
#copied and adapted from https://github.com/cran/Matrix/blob/4936a1e23b6753683b1ca84aefb1a7c6c67018b0/inst/scripts/ssget.sh
pkg=RBCFLib
ssdir=$(dirname $0)/src/SuiteSparse
ssver=7.6.0
sspfx=SuiteSparse-${ssver}
sstgz=${sspfx}.tar.gz
ssurl=https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/refs/tags/v${ssver}.tar.gz
[[ -d ${ssdir} ]] || mkdir -p ${ssdir}
rm -rf ${ssdir}/*
echo "Downloading ${ssurl}"
wget -O ${sstgz} ${ssurl} || exit 1
echo "Extracting only required SuiteSparse components into src/SuiteSparse"
tar -xzf ${sstgz} --strip-components=1 -C ${ssdir} \
    ${sspfx}/CHOLMOD \
    ${sspfx}/AMD \
    ${sspfx}/CAMD \
    ${sspfx}/CCOLAMD \
    ${sspfx}/COLAMD \
    ${sspfx}/CXSparse \
    ${sspfx}/SuiteSparse_config \
    ${sspfx}/LICENSE.txt \
    ${sspfx}/README.md \
    ${sspfx}/Makefile \
    ${sspfx}/CMakeLists.txt
rm -rf ${sstgz}
# Remove unwanted files (images, pdfs, etc.)
find ${ssdir} -iname "*.pdf" -exec rm -f {} +
find ${ssdir} -iname "*.ps" -exec rm -f {} +
find ${ssdir} -iname "*.png" -exec rm -f {} +
find ${ssdir} -iname "*.jpg" -exec rm -f {} +
find ${ssdir} -iname "*.jpeg" -exec rm -f {} +
find ${ssdir} -iname "*.tif" -exec rm -f {} +
find ${ssdir} -iname "*.tiff" -exec rm -f {} +
find ${ssdir} -iname "*.eps" -exec rm -f {} +
echo "Done"
exit 0