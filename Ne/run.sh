#!/bin/bash

# numbers according to model and resource
num_procs=544 # number of processors
num_chains=27
num_atoms=54000
num_atoms_per_mol=2000
griddel=0.5
rvdw=0.95

# do not change the following lines
max_rows=10000000
c_pi=3.141592

rm var.*

echo -e '\t'parameter'('max_rows=$max_rows')' > var.default
echo -e '\t'parameter'('c_pi=$c_pi')' >> var.default
echo -e '\t'parameter'('num_atoms=$num_atoms')' > var.options
echo -e '\t'parameter'('num_atoms_per_mol=$num_atoms_per_mol')' >> var.options
echo -e '\t'parameter'('griddel=$griddel')' >> var.options
echo -e '\t'parameter'('rvdw=$rvdw')' >> var.options

num=1
while [ $num -le $num_chains ]
do
	echo $num > chain.txt
	ibrun -np $num_procs a.out > stdout.$num &&

	num=$(($num+1))
	rm chain.txt
done
