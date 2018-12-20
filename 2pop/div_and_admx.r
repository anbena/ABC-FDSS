#wrapper script to simulate two populations divergence with an asymmetrical single pulse of admixture
#Rscript --vanilla ./m1.r <nchr> <locuslength> <nlociXsim> <recrate=1.12e-8> <nsimulations>
args<-commandArgs(trailingOnly=TRUE)
ms<-"/opt/software/genetics/ms/ms"
cpd<-"./"
mod<-"m4"
nchr<-as.character(args[1])
tgen<-1
mu<-1.25e-8
recomb<-as.numeric(args[4])
ll<-as.numeric(args[2])#locus length
nsims<-as.numeric(args[5])#number of ABC simulations
nloci<-as.numeric(args[3])#loci to simulate in each sim
out<-paste(mod,"_ll",as.character(ll),"_nl",as.character(nloci),"_r",as.character(recomb),"_nc",nchr,sep="")
#main param
n1<-sample(500:50000,nsims,replace=T)
n2<-sample(500:50000,nsims,replace=T)
na<-sample(500:50000,nsims,replace=T)
td2<-sample(50:2500,nsims,replace=T)
td<-td2+sample(250:17500,nsims,replace=T)
#admixture rates follow an exponential distribution with mean 5% truncated in the interval [1%-25%]
adm12<-rexp(nsims*3,1/0.1)
adm12<-adm12[adm12>=0.01 & adm12<=0.25][1:nsims]
adm12<-1-adm12
adm21<-rexp(nsims*3,1/0.1)
adm21<-adm21[adm21>=0.01 & adm21<=0.25][1:nsims]
adm21<-1-adm21
#param transformations
tn1<-4*n1*mu*ll
rn2<-n2*4*mu*ll/tn1
rna<-na*4*mu*ll/tn1
td2s<-(td2/tgen)/(4*n1)
tds<-(td/tgen)/(4*n1)
srec<-4*n1*(recomb*(ll-1))

partable<-cbind(n1,n2,na,td2,td,1-adm12,1-adm21)
colnames(partable)<-c("n1","n2","na","td2","td","adm12","adm21")
partablescaled<-cbind(tn1,rn2,tds,td2s,rna,srec)
write.table(partable,paste(out,".param",sep=""),row.names=F,quote=F,sep="\t")
write.table(partablescaled,paste(out,".paramscaled",sep=""),row.names=F,col.names=T,quote=F,sep="\t")
i<-1
write(paste(ms," ",as.character(2*as.numeric(nchr))," ",as.character(nloci)," -t ",as.character(tn1[i])," -I 2 ",nchr," ",nchr," -n 2 ",as.character(rn2[i])," -es ",as.character(td2s[i])," 2 ",as.character(adm12[i])," -es ",as.character(td2s[i])," 1 ",as.character(adm21[i])," -ej ",as.character(td2s[i]+0.000001)," 3 1 -ej ",as.character(td2s[i]+0.000001)," 4 2 -ej ",as.character(tds[i])," 2 1 -en ",as.character(tds[i])," 1 ",as.character(rna[i]),"| ",cpd,"compute_pd.py -np 2 -nc ",nchr," -w 30 -b 50 -s > ",out,".tab",sep=""),stderr())
for (i in 1:nsims){
	print(i)
	if (i==1){
		system(paste(ms," ",as.character(2*as.numeric(nchr))," ",as.character(nloci)," -t ",as.character(tn1[i])," -I 2 ",nchr," ",nchr," -n 2 ",as.character(rn2[i])," -es ",as.character(td2s[i])," 2 ",as.character(adm12[i])," -es ",as.character(td2s[i])," 1 ",as.character(adm21[i])," -ej ",as.character(td2s[i]+0.000001)," 3 1 -ej ",as.character(td2s[i]+0.000001)," 4 2 -ej ",as.character(tds[i])," 2 1 -en ",as.character(tds[i])," 1 ",as.character(rna[i]),"| ",cpd,"compute_pd.py -np 2 -nc ",nchr," -w 30 -b 50 -s > ",out,".tab",sep=""))
	}
	else{
		system(paste(ms," ",as.character(2*as.numeric(nchr))," ",as.character(nloci)," -t ",as.character(tn1[i])," -I 2 ",nchr," ",nchr," -n 2 ",as.character(rn2[i])," -es ",as.character(td2s[i])," 2 ",as.character(adm12[i])," -es ",as.character(td2s[i])," 1 ",as.character(adm21[i])," -ej ",as.character(td2s[i]+0.000001)," 3 1 -ej ",as.character(td2s[i]+0.000001)," 4 2 -ej ",as.character(tds[i])," 2 1 -en ",as.character(tds[i])," 1 ",as.character(rna[i]),"| ",cpd,"compute_pd.py -np 2 -nc ",nchr," -w 30 -b 50 -s >> ",out,".tab",sep=""))
	}
}
