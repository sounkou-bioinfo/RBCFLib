#!/bin/env bash
#copied and adapted from https://github.com/cran/Matrix/blob/4936a1e23b6753683b1ca84aefb1a7c6c67018b0/inst/scripts/ssget.sh
pkg=RBCFLib
ssdir=$(dirname $0)/src/SuiteSparse
ssver=7.6.0
sspfx=SuiteSparse-${ssver}
sstgz=${sspfx}.tar.gz
ssurl=https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/refs/tags/v${ssver}.tar.gz
[[ -d ${ssdir} ]] || mkdir -p ${ssdir}
cd ${ssdir} || exit 1
rm -rf $PWD/*
echo "Downloading ${ssurl}"
wget -O ${sstgz} ${ssurl} || exit 1
echo "Unpacking ${sstgz}"
tar -xzf ${sstgz} || exit 1
rm -rf ${sstgz}
mv SuiteSparse-${ssver}/* . || exit 1
rm -rf SuiteSparse-${ssver} || exit 1
find SparseSuite -iname "*pdf" | xargs -I {} rm -f {}
echo "Done"
cd - || exit 1
exit 0