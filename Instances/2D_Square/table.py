import os, sys


instances = ['S18_04', 'S18_08', 'S18_09', 'S20_09', 'S20_10', 'S24_09', 'S25_08', 'S36_14', 'S48_23', 'S50_21', 'S60_36', 'S64_42', 'S85_53', 'S100_48', 'S100_50']

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



#~ print "\\begin{table}[h]\n\centering\n\small\n\\begin{tabular}{|c|p{0.5\\textwidth}|c|c|c|}\n\hline"

#~ print "\\rowcolor[rgb]{0,0,0} \\textcolor{white}{\\textbf{Acronym}} & \centering\\textcolor{white}{\\textbf{Sequence}}  & \\textcolor{white}{\\textbf{Len.}}  & \\textcolor{white}{$E^*$} & \\textcolor{white}{\\textbf{References}}\\\\ \hline"

#~ for instance in instances:		

	#~ fp = open("%s.hp" % instance, "r")
	#~ lineas = fp.readlines()
	#~ fp.close()
	
	#~ ['18\n', 'HHPPPPPHHPPPHPPPHP\n', '-4\n', '\n', 'Krasnogor2002\n', 'Cutello2007\n', 'Islam2009\n', '\n', '\n']
	
	#~ l=int(lineas[0])
	#~ seq=lineas[1][:-1]
	#~ e=int(lineas[2])
	#~ ref=""
	#~ for i in range(3, len(lineas)):
		#~ if lineas[i]!='\n':
			#~ ref += "%s, " % lineas[i][:-1]
			
	#~ ref = "\cite{"+ref[:-2]+"}"
	
	#~ print "%s & %s & %d & %d & %s \\\\ \hline" % (instance.replace("_", "\_"), seq, l, e, ref)	
	
#~ print "\end{tabular}\n\caption{HP benchmark sequences adopted for this study.}\label{tab:sequences}\n\end{table}"
	
	
	
	
	
	
	
	
	


#~ S20\_09 & HPHPPHHPHPPHPHHPPHPH & 20 & -9 & \\\hline
#~ S20\_10 & HHHPPHPHPHPPHPHPHPPH & 20 & -10 & \\\hline
#~ S24\_09 & HHPPHPPHPPHPPHPPHPPHPPHH & 24 &  -9 &  \\\hline
#~ S25\_08 & PPHPPHHPPPPHHPPPPHHPPPPHH& 25 & -8 & \\\hline
#~ S36\_14 & PPPHHPPHHPPPPPHHHHHHHPPHHPPPP HHPPHPP & 36 & -14 & \\\hline
#~ S48\_23 & PPHPPHHPPHHPPPPPHHHHHHHHHHPPPPPPHHPP...& 48 & -23 & \\\hline
#~ S50\_21 & HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPH... & 50 & -21 & \\\hline
#~ S60\_36 & PPHHHPHHHHHHHHPPPHHHHHHHHHHPHPPPHHHH... & 60 & -36 & \\\hline
#~ S64\_42 & HHHHHHHHHHHHPHPHPPHHPPHHPPHPPHHPPHHP... & 64 & -42 & \\\hline
#~ S85\_53 & HHHHPPPPHHHHHHHHHHHHPPPPPPHHHHHHHHHH... & 85 & -53 & \\\hline
#~ S100\_48 & PPPPPPHPHHPPPPPHHHPHHHHHPHHPPPPHHPP... & 100 & -48 & \\\hline
#~ S100\_50 & PPPHHPPHHHHPPHHHPHHPHHPHHHHPPPPPPPP... & 100 & -50 & \\\hline
 #~ 