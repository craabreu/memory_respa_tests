#!/usr/bin/env bash
dir=$1
name=dt${2}fs

# Compute radial distribution functions:
cd $dir/$name
mkdir -p adf
cd adf
root=../../../..
#$root/travis/exe/travis -i $root/water/travis_angle.inp -p ../$name.xyz > ${name}_angle.output
file="adf_H2O_#2_[O1r_H1r]-[O1r_H2r].csv"
sed 's/;/,/g' $file > $root/water/$dir/results/${name}_angle.csv

values=$(grep "Mean value" ${name}_angle.output | sed -e 's/:/,/g' -e 's/degree/,/g' | cut -d"," -f2,4 | sed 's/ //g')
dt=$(echo $2 | sed 's/p/./')
line=$(echo "$dt,$values")
eval "sed -i 's/^$dt.*/$line/' $root/water/$dir/results/angle_stats.csv"
