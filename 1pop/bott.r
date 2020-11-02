#wrapper script to simulate one population bottleneck
#Rscript --vanilla ./bott.r <nchr> <locuslength> <nlociXsim> <rec.rate=1e-8> <nsimulations>
args<-commandArgs(trailingOnly=TRUE)
ms<-"/opt/software/genetics/ms/ms"
cpd<-"./"
mod<-"bott"
nchr<-as.character(args[1])
tgen<-1
mu<-1e-8
recomb<-as.numeric(args[4])#recombination rate
ll<-as.numeric(args[2])#locus length
nsims<-as.numeric(args[5])#number of ABC simulations
nloci<-as.numeric(args[3])#loci to simulate in each sim
out<-paste(mod,"_ll",as.character(ll),"_nl",as.character(nloci),"_r",as.character(recomb),"_nc",nchr,sep="")

#Prior distributions parameters
na<-sample(500:50000,nsims,replace=T)
intensity<-runif(nsims,min=10,max=100)
n1<-na/intensity
tb<-sample(100:20000,nsims,replace=T)

#Parameters transformations for ms simulator
tn1<-4*n1*mu*ll
tbs<-(tb/tgen)/(4*n1)
srec<-4*n1*(recomb*(ll-1))

#Output: parameters files
partable<-cbind(n1,intensity,tb,na)
colnames(partable)<-c("n1","intensity","tb","na")
partablescaled<-cbind(tn1,tbs,srec)

write.table(partable,paste(out,".param",sep=""),row.names=F,quote=F,sep="\t")
write.table(partablescaled,paste(out,".paramscaled",sep=""),row.names=F,col.names=T,quote=F,sep="\t")

#Output summary statistics: FDSS
i<-1
write(paste(ms," ",as.character(nchr)," ",as.character(nloci)," -t ",as.character(tn1[i])," -r ",as.character(srec[i])," ",as.character(ll)," -eN ",as.character(tbs[i])," ",as.character(intensity[i])," | ",cpd,"compute_ss.py -np 1 -nc ",nchr," -w 100 -b 20 > ",out,".tab",sep=""),stderr())
for (i in 1:nsims){
	print(i)
	if (i==1){
		system(paste(ms," ",as.character(nchr)," ",as.character(nloci)," -t ",as.character(tn1[i])," -r ",as.character(srec[i])," ",as.character(ll)," -eN ",as.character(tbs[i])," ",as.character(intensity[i])," | ",cpd,"compute_ss.py -np 1 -nc ",nchr," -w 100 -b 20 > ",out,".tab",sep=""))
	}
	else{
		system(paste(ms," ",as.character(nchr)," ",as.character(nloci)," -t ",as.character(tn1[i])," -r ",as.character(srec[i])," ",as.character(ll)," -eN ",as.character(tbs[i])," ",as.character(intensity[i])," | ",cpd,"compute_ss.py -np 1 -nc ",nchr," -w 100 -b 20 >> ",out,".tab",sep=""))
	}
}
