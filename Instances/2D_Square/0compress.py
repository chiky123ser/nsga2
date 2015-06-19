import os, sys




sequences = ["HHPPPPPHHPPPHPPPHP",
			"HPHPHHHPPPHHHHPPHH",
			"PHPPHPHHHPHHPHHHHH",
			"HPHPPHHPHPPHPHHPPHPH",

text=""

i=0
j=0
while (i<len(sequence)):
	frec=0
	j=i
	while(j<len(sequence)):
		
		if sequence[j] == sequence[i]: 
			frec+=1
		else: break
		j+=1		
	text += "%d%c" % (frec,sequence[i])
	i = j
	
print text