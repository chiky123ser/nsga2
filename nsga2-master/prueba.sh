#! /bin/bash

clear
LT=/home/servando/Documents/Multiobjetivizacion/Lattices/2D_Square.lat
PL=/home/servando/Documents/Multiobjetivizacion/Instances/2D_Square

for j in 4
do
	echo '\n'
	for name in DEC1
	do
		for k in 5 #0 1 2 3 4 5
		do
			for i in 20
			do
	
			./nsga2r .5 $LT $PL/2d$j.hp MO $name DYN30 $i 4 1 2 2 .9 .5 $j $i $k

			done
		done
	done
done

