#!/bin/bash

case+=("01")
nsteps+=("3600000")
device+=("0")
secdev+=("1")

case+=("03")
nsteps+=("1200000")
device+=("2")
secdev+=("3")

case+=("06")
nsteps+=("600000")
device+=("0")
secdev+=("1")

case+=("09")
nsteps+=("400000")
device+=("2")
secdev+=("3")

case+=("15")
nsteps+=("240000")
device+=("0")
secdev+=("1")

case+=("30")
nsteps+=("120000")
device+=("2")
secdev+=("3")

case+=("45")
nsteps+=("80000")
device+=("0")
secdev+=("1")

case+=("90")
nsteps+=("40000")
device+=("2")
secdev+=("3")

name=dt${case[$1]}fs
cd $name
python simulate.py --device ${device[$1]} --secdev ${secdev[$1]} --timestep ${case[$1]} --nsteps ${nsteps[$1]}

properties='Temp PotEng Press MolPress'

dt=$(echo ${case[$1]} | sed 's/p/./')
values=$(postlammps -d comma -mm -c 83 -nt ineff $properties < ${name}.csv | \
         cut -d',' -f2,3 | \
         tr '\n' ',' | \
         sed 's/,$//')
line=$(echo "$dt,$values")
eval "sed -i 's/^$dt.*/$line/' ../results/properties.csv"
