#!/bin/bash

name=dt0p5fs
python simulate.py --device 0 --secdev 1 --nsteps 7200000

properties='Temp PotEng Press MolPress'

dt="0.5"
values=$(postlammps -d comma -mm -c 83 -nt ineff $properties < ${name}.csv | \
         cut -d',' -f2,3 | \
         tr '\n' ',' | \
         sed 's/,$//')
line=$(echo "$dt,$values")
eval "sed -i 's/^$dt.*/$line/' ../results/properties.csv"
echo $line
