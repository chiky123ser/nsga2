#! /bin/bash

clear
LT=/home/servando/Documents/Multiobjetivizacion/Lattices/2D_Square.lat
PL=/home/servando/Documents/Multiobjetivizacion/Instances/2D_Square
for j in 4 #3 5 6 7 8 9 10 11 12 13 14 15
do
	echo '\n'
	for name in LOCALITY #DEC1 DEC2
	do
		echo $name
		for k in 0 1 2 3 4 5
		do
			echo '\t'$k
			#
			for i in 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 #1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 # 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100
			do
			echo '\t'$i
			./nsga2r .5 $LT $PL/2d$j.hp MO $name DYN30 $i 200 1000 2 2 .8 .5 $j $i $k
			done
		done
	done
done
echo 'finish'

#debug
#run .5 /home/servando/Documents/Multiobjetivizacion/Lattices/2D_Square.lat /home/servando/Documents/Multiobjetivizacion/Instances/2D_Square/2d4.hp MO LOCALITY DYN30 45 20 1000 2 2 .9 .5 4 45 2
#shutdown -h now
