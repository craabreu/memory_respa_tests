#!/bin/bash

case+=("01")
nsteps+=("1800000")
device+=("0")

case+=("03")
nsteps+=("600000")
device+=("1")

case+=("06")
nsteps+=("300000")
device+=("2")

case+=("09")
nsteps+=("200000")
device+=("3")

case+=("15")
nsteps+=("1200000")
device+=("0")

case+=("30")
nsteps+=("60000")
device+=("1")

case+=("45")
nsteps+=("40000")
device+=("2")

case+=("90")
nsteps+=("20000")
device+=("3")

name=dt${case[$1]}fs
cd $name
python ../solvation.py --device ${device[$1]} --timestep ${case[$1]} --nsteps ${nsteps[$1]}
