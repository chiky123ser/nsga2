data=c(NULL)
funcion=c("Parity","Locality")
cross=c("Cruce","Cycle","Non_wrap","One","Static","Two")
imp="/home/servando/Documents/Frankenstein/nsga2-master/output/R/Impresiones/"
#Indice para for de documentos
x = seq(from=1, to=100, by=1)
#Indice para funciones
indFun=seq(from=1,to=2,by=1)
#Indice para cross
indCross=seq(from=1,to=6,by=1)


for (i in indFun) {

	for (j in indCross) {
		file="/home/servando/Documents/Frankenstein/nsga2-master/output/"
		file=paste(file,paste(funcion[i],"/",sep=""),sep="")
		file=paste(file,cross[j],sep="")
		for (k in x) {
			
			fileX=paste(file,"/best/a1_best_pop.out4_",sep="")
			fileX=paste(fileX,k,sep="")

			datos= read.table(fileX)
			datos2=datos["V1"]+datos["V2"]
			minimo=min(datos2)
			data=c(data,minimo)
			#print(i)
			#print(minimo)
		}
		newmin=min(data)
		mean=mean(data)
		desviacion=sd(data)	
		#print(rle(sort(data)))
		frecuencia=sum(data==newmin)
		t=paste(imp,paste(funcion[i],cross[j],sep="_"),sep="")
		write(data,t)
		write(newmin,t, append = TRUE)
		write(mean,t, append = TRUE)
		write(desviacion,t, append = TRUE)
		write(frecuencia,t, append = TRUE)
		
		data=c(NULL)

		#print(data)
	}
}

Cruce_P=scan("/home/servando/Documents/Frankenstein/nsga2-master/output/R/Impresiones/Parity_Cruce", nmax = 100)
Cycle_P=scan("/home/servando/Documents/Frankenstein/nsga2-master/output/R/Impresiones/Parity_Cycle", nmax = 100)
Non_wrap_P=scan("/home/servando/Documents/Frankenstein/nsga2-master/output/R/Impresiones/Parity_Non_wrap", nmax = 100)
One_P=scan("/home/servando/Documents/Frankenstein/nsga2-master/output/R/Impresiones/Parity_One", nmax = 100)
Static_P=scan("/home/servando/Documents/Frankenstein/nsga2-master/output/R/Impresiones/Parity_Static", nmax = 100)
Two_P=scan("/home/servando/Documents/Frankenstein/nsga2-master/output/R/Impresiones/Parity_Two", nmax = 100)

Cruce_L=scan("/home/servando/Documents/Frankenstein/nsga2-master/output/R/Impresiones/Locality_Cruce", nmax = 100)
Cycle_L=scan("/home/servando/Documents/Frankenstein/nsga2-master/output/R/Impresiones/Locality_Cycle", nmax = 100)
Non_wrap_L=scan("/home/servando/Documents/Frankenstein/nsga2-master/output/R/Impresiones/Locality_Non_wrap", nmax = 100)
One_L=scan("/home/servando/Documents/Frankenstein/nsga2-master/output/R/Impresiones/Locality_One", nmax = 100)
Static_L=scan("/home/servando/Documents/Frankenstein/nsga2-master/output/R/Impresiones/Locality_Static", nmax = 100)
Two_L=scan("/home/servando/Documents/Frankenstein/nsga2-master/output/R/Impresiones/Locality_Two", nmax = 100)


#Cruce_D=scan("/home/servando/Documents/Frankenstein/nsga2-master/output/R/Impresiones/Dec_Cruce", nmax = 100)
#Cycle_D=scan("/home/servando/Documents/Frankenstein/nsga2-master/output/R/Impresiones/Dec_Cycle", nmax = 100)
#Non_wrap_D=scan("/home/servando/Documents/Frankenstein/nsga2-master/output/R/Impresiones/Dec_Non_wrap", nmax = 100)
#One_D=scan("/home/servando/Documents/Frankenstein/nsga2-master/output/R/Impresiones/Dec_One", nmax = 100)
#Static_D=scan("/home/servando/Documents/Frankenstein/nsga2-master/output/R/Impresiones/Dec_Static", nmax = 100)
#Two_D=scan("/home/servando/Documents/Frankenstein/nsga2-master/output/R/Impresiones/Dec_Two", nmax = 100)


frame1=data.frame(Cruce_P,Cycle_P,Non_wrap_P,One_P,Static_P,Two_P)
frame2=data.frame(Cruce_L,Cycle_L,Non_wrap_L,One_L,Static_L,Two_L)
#frame3=data.frame(Cruce_D,Cycle_D,Non_wrap_D,One_D,Static_D,Two_D)
frame4=data.frame(frame1,frame2)

#box1=boxplot(frame1,las = 2,at = c(1,2,3, 5,6,7),par(mar = c(12, 5, 4, 2)+ 0.1))
#boxplot(frame2,las = 2,at = c(1,2,3, 5,6,7),par(mar = c(12, 5, 4, 2)+ 0.1))
#boxplot(frame3,las = 2,at = c(1,2,3, 5,6,7),par(mar = c(12, 5, 4, 2)+ 0.1))
#boxplot(frame4,las = 2,at = c(1,2,3,4,5,6, 8,9,10,11,12,13, 15,16,17,18,19,20),par(mar = c(12, 5, 4, 2)+ 0.1))
