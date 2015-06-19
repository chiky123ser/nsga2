#! /bin/bash
clear
for j in 4
do
	for i in 17 25 65 92
	do
		
	../nsga2r .5  /home/servando/Documents/Multiobjetivizacion/Lattices/2D_Square.lat /home/servando/Documents/Multiobjetivizacion/Instances/2D_Square/2d$j.hp MO PARITY DET $i 200 1000 2 2 .9 .5 $j $i
	done
done
echo 'finish'

