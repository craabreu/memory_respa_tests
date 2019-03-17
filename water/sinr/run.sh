#!/bin/bash

case+=("0p5")
nsteps+=("7200000")
device+=("0")
secdev+=("1")

case+=("01")
nsteps+=("3600000")
device+=("2")
secdev+=("3")

case+=("03")
nsteps+=("1200000")
device+=("0")
secdev+=("1")

case+=("06")
nsteps+=("600000")
device+=("2")
secdev+=("3")

case+=("09")
nsteps+=("400000")
device+=("0")
secdev+=("1")

case+=("15")
nsteps+=("240000")
device+=("2")
secdev+=("3")

case+=("30")
nsteps+=("120000")
device+=("0")
secdev+=("1")

case+=("45")
nsteps+=("80000")
device+=("2")
secdev+=("3")

case+=("90")
nsteps+=("40000")
device+=("0")
secdev+=("1")

name=sinr-${case[$1]}fs
cd $name
#python simulate.py --device ${device[$1]} --secdev ${secdev[$1]} --timestep ${case[$1]} --nsteps ${nsteps[$1]}
#babel -ipdb $name.pdb -oxyz $name.xyz
#sed -i -e 's/Ow/O/g' -e 's/Hw/H/g' $name.xyz

properties='Temp PotEng Press MolPress Virial NBVirial MolVirial'
# title='dt'
# for prop in $properties; do
#     title+=,$prop,rmse[$prop]
# done
# echo $title

dt=$(echo ${case[$1]} | sed 's/p/./')
values=$(postlammps -d comma -mm -c 83 -nt ineff $properties < ${name}.csv | \
         cut -d',' -f2,3 | \
         tr '\n' ',' | \
         sed 's/,$//')
line=$(echo "$dt,$values")
eval "sed -i 's/^$dt.*/$line/' ../results/properties.csv"
