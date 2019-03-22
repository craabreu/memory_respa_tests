#!/usr/bin/env bash
dir=$1
name=dt${2}fs

# Compute radial distribution functions:
cd $dir/$name
mkdir -p ddf
cd ddf
root=../../../..
$root/travis/exe/travis -i $root/ethylene_glycol/travis_dihedral.inp -p ../$name.xyz > ${name}_dihedral.output

file="ddf_C2H6O2_#2_[C1r_O1r]-[C2r_O2r]-[C1r_C2r].csv"
cut -d";" -f1,2 $file > $root/ethylene_glycol/$dir/results/${name}_occo_dihedral.csv
sed -i 's/;/,/g' $root/ethylene_glycol/$dir/results/${name}_occo_dihedral.csv

file1="ddf_C2H6O2_#2_[C1r_C2r]-[O1r_H1r]-[C1r_O1r].csv"
file2="ddf_C2H6O2_#2_[C2r_C1r]-[O2r_H2r]-[C2r_O2r].csv"
paste -d, <(cut -d";" -f1,2 $file1 | sed 's/;/,/g') <(cut -d";" -f2 $file2) > $root/ethylene_glycol/$dir/results/${name}_ccoh_dihedral.csv
sed -i 's/Occurrence,  Occurrence/Occurrence1,  Occurrence2/g' $root/ethylene_glycol/$dir/results/${name}_ccoh_dihedral.csv

file=${name}_dihedral.output
file=travis.log
values=$(grep "Mean value" $file | sed -e 's/:/,/g' -e 's/pm/,/g' | cut -d"," -f2,4 | sed 's/ //g')
values=$(echo $values)
dt=$(echo $2 | sed 's/p/./')
line=$(echo "$dt,$values" | sed -e 's/ /,/g' | sed 's/,$//')
eval "sed -i 's/^$dt.*/$line/' $root/ethylene_glycol/$dir/results/dihedral_stats.csv"
