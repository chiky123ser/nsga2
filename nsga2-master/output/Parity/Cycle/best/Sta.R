data=c(NULL)
x = seq(from=1, to=100, by=1)
for (i in x) {
file="a1_best_pop.out4_"
file=paste(file,i,sep="")
datos= read.table(file)
datos2=datos["V1"]+datos["V2"]
minimo=min(datos2)
data=c(data,minimo)
print(i)
print(minimo)

}
print(data)
boxplot(data)