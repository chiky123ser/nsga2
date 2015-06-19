#! /bin/bash

clear
LT=/home/servando/Documents/Multiobjetivizacion/Lattices/2D_Square.lat
for j in 4
do
	echo '\n'
	for i in 20
	do
	
	./nsga2r .5 $LT /home/servando/Documents/Multiobjetivizacion/Instances/2D_Square/2d$j.hp MO PARITY DYN30 $i 200 1000 2 2 .9 .5 $j $i
	done
done

