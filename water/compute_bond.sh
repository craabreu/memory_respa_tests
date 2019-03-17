#!/usr/bin/env bash
dir=$1
name=dt${2}fs

# Compute radial distribution functions:
cd $dir/$name
mkdir -p bdf
cd bdf
root=../../../..
#$root/travis/exe/travis -i $root/water/travis_bond.inp -p ../$name.xyz > ${name}_bond.output
file="rdf_H2O_#2_[Or_Hr].csv"
sed 's/;/,/g' $file > $root/water/$dir/results/${name}_bond.csv

values=$(grep "Mean value" ${name}_bond.output | sed -e 's/:/,/g' -e 's/pm/,/g' | cut -d"," -f2,4 | sed 's/ //g')
dt=$(echo $2 | sed 's/p/./')
line=$(echo "$dt,$values")
echo $line
eval "sed -i 's/^$dt.*/$line/' $root/water/$dir/results/bond_stats.csv"
