#!/bin/bash
# Shear series, flow x / gradient y; two runs at a time. Output: stress_<tag>.dat
cd "$(dirname "$0")"
[ -f kg_equil.data ] || lmp -in in.equil -log log.equil > out.equil 2>&1
run(){ r=$1; tag=r$(echo $r | tr -d '.'); lmp -in in.shear -var rate $r -var tag $tag -var seed 11 -var nsteps 300000 -log log.$tag > out.$tag 2>&1; }
run 0.02 & run 0.05 & wait
run 0.1 & run 0.2 & wait
