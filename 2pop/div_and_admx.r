#wrapper script to simulate two populations divergence with an asymmetrical single pulse of admixture
#Rscript --vanilla ./div_and_admx.r <nchr> <locuslength> <nlociXsim> <rec.rate=1e-8> <nsimulations>
args<-commandArgs(trailingOnly=TRUE)
ms<-"/opt/software/genetics/ms/ms"
cpd<-"./"
mod<-"m3"
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
n2<-sample(500:50000,nsims,replace=T)
na<-sample(500:50000,nsims,replace=T)
td2<-sample(50:2500,nsims*10,replace=T)#admixture time
td<-td2+sample(250:17500,nsims*10,replace=T)#populations divergence time
#Admixture time was sampled only if: (Tsep-Tadm)/Tsep was between 0.2 and 0.8
a<-(td-td2)/td
td<-td[a>=0.2 & a<=0.8][1:nsims]
td2<-td2[a>=0.2 & a<=0.8][1:nsims]

#admixture rates follow an uniform distribution with range [5%-20%]
adm12<-runif(nsims,0.05,0.20)
adm12<-1-adm12
adm21<-runif(nsims,0.05,0.20)
adm21<-1-adm21

#Parameters transformation for ms simulator
tn1<-4*n1*mu*ll
rn2<-n2*4*mu*ll/tn1
rna<-na*4*mu*ll/tn1
td2s<-(td2/tgen)/(4*n1)
tds<-(td/tgen)/(4*n1)
srec<-4*n1*(recomb*(ll-1))

#Output: parameters files
partable<-cbind(n1,n2,na,td2,td,1-adm12,1-adm21)
colnames(partable)<-c("n1","n2","na","td2","td","adm12","adm21")
partablescaled<-cbind(tn1,rn2,tds,td2s,rna,srec)
write.table(partable,paste(out,".param",sep=""),row.names=F,quote=F,sep="\t")
write.table(partablescaled,paste(out,".paramscaled",sep=""),row.names=F,col.names=T,quote=F,sep="\t")

#Output summary statistics: FDSS
i<-1
write(paste(ms," ",as.character(2*as.numeric(nchr))," ",as.character(nloci)," -t ",as.character(tn1[i])," -I 2 ",nchr," ",nchr," -n 2 ",as.character(rn2[i])," -es ",as.character(td2s[i])," 2 ",as.character(adm12[i])," -es ",as.character(td2s[i])," 1 ",as.character(adm21[i])," -ej ",as.character(td2s[i]+0.000001)," 3 1 -ej ",as.character(td2s[i]+0.000001)," 4 2 -ej ",as.character(tds[i])," 2 1 -en ",as.character(tds[i])," 1 ",as.character(rna[i]),"| ",cpd,"compute_ss.py -np 2 -nc ",nchr," -w 30 -b 50 -s > ",out,".tab",sep=""),stderr())
for (i in 1:nsims){
	print(i)
	if (i==1){
		system(paste(ms," ",as.character(2*as.numeric(nchr))," ",as.character(nloci)," -t ",as.character(tn1[i])," -I 2 ",nchr," ",nchr," -n 2 ",as.character(rn2[i])," -es ",as.character(td2s[i])," 2 ",as.character(adm12[i])," -es ",as.character(td2s[i])," 1 ",as.character(adm21[i])," -ej ",as.character(td2s[i]+0.000001)," 3 1 -ej ",as.character(td2s[i]+0.000001)," 4 2 -ej ",as.character(tds[i])," 2 1 -en ",as.character(tds[i])," 1 ",as.character(rna[i]),"| ",cpd,"compute_ss.py -np 2 -nc ",nchr," -w 30 -b 50 -s > ",out,".tab",sep=""))
	}
	else{
		system(paste(ms," ",as.character(2*as.numeric(nchr))," ",as.character(nloci)," -t ",as.character(tn1[i])," -I 2 ",nchr," ",nchr," -n 2 ",as.character(rn2[i])," -es ",as.character(td2s[i])," 2 ",as.character(adm12[i])," -es ",as.character(td2s[i])," 1 ",as.character(adm21[i])," -ej ",as.character(td2s[i]+0.000001)," 3 1 -ej ",as.character(td2s[i]+0.000001)," 4 2 -ej ",as.character(tds[i])," 2 1 -en ",as.character(tds[i])," 1 ",as.character(rna[i]),"| ",cpd,"compute_ss.py -np 2 -nc ",nchr," -w 30 -b 50 -s >> ",out,".tab",sep=""))
	}
}
