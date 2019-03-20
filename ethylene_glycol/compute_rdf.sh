#!/usr/bin/env bash
dir=$1
name=dt${2}fs

# Compute radial distribution functions:
cd $dir/$name
mkdir -p rdf
cd rdf
root=../../../..
$root/travis/exe/travis -i $root/water/travis_rdf.inp -p ../$name.xyz #> ${name}_rdf.output
gOO="rdf_H2O_#2_H2O_[Or_Oo].csv"
gOH="rdf_H2O_#2_H2O_[Or_Ho].csv"
gHH="rdf_H2O_#2_H2O_[Hr_Ho].csv"
paste -d"," <(cut -d";" -f1 $gOO | sed 's/# //g')  \
            <(cut -d";" -f2 $gOO | sed 's/g(r)/g(O-O)/g') \
            <(cut -d";" -f2 $gOH | sed 's/g(r)/g(O-H)/g') \
            <(cut -d";" -f2 $gHH | sed 's/g(r)/g(H-H)/g') > $root/water/$dir/results/${name}_rdf.csv
