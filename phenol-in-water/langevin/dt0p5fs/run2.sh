#!/bin/bash
device=$1
part=$(expr $1 + 5)
python solvation.py --device $device --part $part --nsteps 3600000
