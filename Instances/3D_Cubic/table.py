import os, sys


instances = ['S20_11', 'S24_13', 'S25_09', 'S36_18', 'S46_32', 'S48_32','S50_31', 'S58_44', 'S60_52','S64_55', 'S67_56', 'S88_72', 'S103_56', 'S124_71', 'S136_80']

print "\\begin{table*}[h]\n\centering\n\small\n\\begin{tabular}{|c|p{4.8in}|c|c|}\n\hline"

print "\\rowcolor[rgb]{0,0,0} & \centering\\textcolor{white}{\\textbf{Sequence}}  & \\textcolor{white}{\\textbf{Len.}}  & \\textcolor{white}{$HHtc^*$}\\\\ \hline"

cont = 0
for instance in instances:		

	cont += 1
	
	fp = open("%s.hp" % instance, "r")
	lineas = fp.readlines()
	fp.close()
	
	l=int(lineas[0])
	sequence=lineas[1][:-1]
	e=int(lineas[2])
	
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
		if frec>1:
			text += "%c$_{%d}$" % (sequence[i], frec)
		else: 
			text += "%c" % (sequence[i])
		i = j
	
	#~ print "\\textcolor{white}{\\textbf{S%d}} & %s & %d & %d \\\\ \hline" % (cont, text, l, e)	
	print "\\textcolor{white}{\\textbf{%s}} & %s & %d & %d \\\\ \hline" % (instance.replace('_', '\_'), text, l, e)	
	
print "\end{tabular}\n\caption{HP benchmark sequences adopted for this study.}\label{tab:sequences}\n\end{table*}"


