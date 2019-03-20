#!/usr/bin/env bash
dir=$1
name=dt${2}fs

# Compute radial distribution functions:
cd $dir/$name
mkdir -p bdf
cd bdf
root=../../../..
$root/travis/exe/travis -i $root/ethylene_glycol/travis_bond.inp -p ../$name.xyz > ${name}_bond.output

file="rdf_C2H6O2_#2_[C1r_C2r].csv"
cut -d";" -f1,2 $file > $root/ethylene_glycol/$dir/results/${name}_cc_bond.csv

file1="rdf_C2H6O2_#2_[C1r_O1r].csv"
file2="rdf_C2H6O2_#2_[C2r_O2r].csv"
paste -d, <(cut -d";" -f1,2 $file1 | sed 's/;/,/g') <(cut -d";" -f2 $file2) > $root/ethylene_glycol/$dir/results/${name}_co_bond.csv

file1="rdf_C2H6O2_#2_[O1r_H1r].csv"
file2="rdf_C2H6O2_#2_[O2r_H2r].csv"
paste -d, <(cut -d";" -f1,2 $file1 | sed 's/;/,/g') <(cut -d";" -f2 $file2) > $root/ethylene_glycol/$dir/results/${name}_oh_bond.csv

file=${name}_bond.output
file=travis.log
values=$(grep "Mean value" $file | sed -e 's/:/,/g' -e 's/pm/,/g' | cut -d"," -f2,4 | sed 's/ //g')
values=$(echo $values)
dt=$(echo $2 | sed 's/p/./')
line=$(echo "$dt,$values" | sed -e 's/ /,/g' | sed 's/,$//')
eval "sed -i 's/^$dt.*/$line/' $root/ethylene_glycol/$dir/results/bond_stats.csv"
