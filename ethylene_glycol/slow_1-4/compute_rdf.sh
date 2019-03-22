#!/usr/bin/env bash
dir=$1
name=dt${2}fs

# Compute radial distribution functions:
cd $dir/$name
mkdir -p rdf
cd rdf
root=../../../..
$root/travis/exe/travis -i $root/ethylene_glycol/travis_rdf.inp -p ../$name.xyz #> ${name}_rdf.output
gOO="rdf_C2H6O2_#2_C2H6O2_[Or_Oo].csv"
gOH="rdf_C2H6O2_#2_C2H6O2_[Or_H1-2o].csv"
gHH="rdf_C2H6O2_#2_C2H6O2_[H1-2r_H1-2o].csv"
paste -d"," <(cut -d";" -f1 $gOO | sed 's/# //g')  \
            <(cut -d";" -f2 $gOO | sed 's/g(r)/g(O-O)/g') \
            <(cut -d";" -f2 $gOH | sed 's/g(r)/g(O-H)/g') \
            <(cut -d";" -f2 $gHH | sed 's/g(r)/g(H-H)/g') > $root/ethylene_glycol/$dir/results/${name}_rdf.csv
