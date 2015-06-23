#! /bin/bash

clear
LT=/home/servando/Documents/Multiobjetivizacion/Lattices/2D_Square.lat
for j in 4
do
	echo '\n'
	for name in LOCALITY DEC1 DEC2
	do
		for k in 0 1 2 3 4 5
		do
			for i in 20
			do
	
			./nsga2r .5 $LT /home/servando/Documents/Multiobjetivizacion/Instances/2D_Square/2d$j.hp MO $name DYN30 $i 200 1000 2 2 .9 .5 $j $i $k
			done
		done
	done
done

