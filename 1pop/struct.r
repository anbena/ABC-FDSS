#wrapper script to simulate one structured population
#Rscript --vanilla ./struct.r <nchr> <locuslength> <nlociXsim> <rec.rate=1e-8> <nsimulations>
args<-commandArgs(trailingOnly=TRUE)
ms<-"/opt/software/genetics/ms/ms"
cpd<-"./"
mod<-"struct"
nchr<-as.character(args[1])
tgen<-1
mu<-1e-8
recomb<-as.numeric(args[4])#recombination rate
ll<-as.numeric(args[2])#locus length
nsims<-as.numeric(args[5])#number of ABC simulations
nloci<-as.numeric(args[3])#loci to simulate in each sim
out<-paste(mod,"_ll",as.character(ll),"_nl",as.character(nloci),"_r",as.character(recomb),"_nc",nchr,sep="")

#Prior distributions parameters
n1<-sample(500:50000,nsims,replace=T)
deme<-sample(2:10,nsims,replace=T)
m<-rexp(nsims,0.1)*4

#Parameters transformations for ms simulator
tn1<-4*n1*mu*ll
srec<-4*n1*(recomb*(ll-1))

#Output: parameters files
partable<-cbind(n1,m,deme)
colnames(partable)<-c("n1","m","deme")
partablescaled<-cbind(tn1,m,srec)
write.table(partable,paste(out,".param",sep=""),row.names=F,quote=F,sep="\t")
write.table(partablescaled,paste(out,".paramscaled",sep=""),row.names=F,col.names=T,quote=F,sep="\t")

#Output summary statistics: FDSS
i<-1
write(paste(ms," ",as.character(nchr)," ",as.character(nloci)," -t ",as.character(tn1[i])," -r ",as.character(srec[i])," ",as.character(ll)," -I ",as.character(deme[i])," ",as.character(nchr)," ",paste(rep("0",(deme[i]-1)),collapse=" ")," ",as.character(m[i])," | ",cpd,"compute_ss.py -np 1 -nc ",nchr," -w 100 -b 20 > ",out,".tab",sep=""),stderr())
for (i in 1:nsims){
	print(i)
	if (i==1){
		system(paste(ms," ",as.character(nchr)," ",as.character(nloci)," -t ",as.character(tn1[i])," -r ",as.character(srec[i])," ",as.character(ll)," -I ",as.character(deme[i])," ",as.character(nchr)," ",paste(rep("0",(deme[i]-1)),collapse=" ")," ",as.character(m[i])," | ",cpd,"compute_ss.py -np 1 -nc ",nchr," -w 100 -b 20 > ",out,".tab",sep=""))
	}
	else{
		system(paste(ms," ",as.character(nchr)," ",as.character(nloci)," -t ",as.character(tn1[i])," -r ",as.character(srec[i])," ",as.character(ll)," -I ",as.character(deme[i])," ",as.character(nchr)," ",paste(rep("0",(deme[i]-1)),collapse=" ")," ",as.character(m[i])," | ",cpd,"compute_ss.py -np 1 -nc ",nchr," -w 100 -b 20 >> ",out,".tab",sep=""))
	}
}
