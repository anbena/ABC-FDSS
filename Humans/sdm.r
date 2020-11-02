#Single dispersal model
#Rscript --vanilla ./sdm.r <nchr> <locuslength> <nlociXsim> <rec.rate=1.12e-8> <nsimulations>

samp_int_vec<-function(x=1,y=1:10){
#x is an integer, y is a vector
out<-c()
for (i in 1:length(y)){
	if (x!=y[i]){
		out[i]<-sample(x:y[i],1,replace=T)
	}else{
		out[i]<-x
	}
}
return(out)
}

samp_vec_vec<-function(x=1:10,y=1:10){
#x is a vector, y is a vector
out<-c()
for (i in 1:length(y)){
	if (x[i]!=y[i]){
		out[i]<-sample(x[i]:y[i],1,replace=T)
	}else{
		out[i]<-x[i]
	}
}
return(out)
}

args<-commandArgs(trailingOnly=TRUE)
ms<-"/opt/software/genetics/msms/bin/msms"
cpd<-"./"
mod<-"sdm"
nchr<-as.character(args[1])
tgen<-29
mu<-1.25e-8
recomb<-as.numeric(args[4])#recombination rate
ll<-as.numeric(args[2])#locus length
nsims<-as.numeric(args[5])#number of ABC simulations
nloci<-as.numeric(args[3])#loci to simulate in each sim
out<-paste(mod,"_ll",as.character(ll),"_nl",as.character(nloci),"_r",as.character(recomb),"_nc",nchr,sep="")

#Prior distributions parameters as in Malaspinas et al. 2016
tbot<-2900
stN<-85735
stD<-67570
nAR<-sample(50:50000,nsims,replace=T)
nD<-sample(50:50000,nsims,replace=T)
nDR<-sample(50:50000,nsims,replace=T)
nN<-sample(50:50000,nsims,replace=T)
nNR<-sample(50:50000,nsims,replace=T)
nY<-sample(50:50000,nsims,replace=T)
nG<-sample(50:50000,nsims,replace=T)
nE<-sample(50:50000,nsims,replace=T)
nA<-sample(50:50000,nsims,replace=T)
nP<-sample(50:50000,nsims,replace=T)
nYG<-sample(50:50000,nsims,replace=T)
nNNR<-sample(50:50000,nsims,replace=T)
nDDR<-sample(50:50000,nsims,replace=T)
nDN<-sample(50:50000,nsims,replace=T)
nADN<-sample(50:50000,nsims,replace=T)
nAM<-sample(50:50000,nsims,replace=T)
rYG<-nYG/nY
rNNR<-nNNR/nN
rDDR<-nDDR/nD
rNDN<-nDN/nDDR
rADN<-nADN/nAR
rAM<-nAM/nADN

m67<-runif(nsims, 10^-6, 10^-3)*4*nAR
m76<-runif(nsims, 10^-6, 10^-3)*4*nAR
m78<-runif(nsims, 10^-6, 10^-3)*4*nAR
m87<-runif(nsims, 10^-6, 10^-3)*4*nAR
m89<-runif(nsims, 10^-6, 10^-3)*4*nAR
m98<-runif(nsims, 10^-6, 10^-3)*4*nAR
m910<-runif(nsims, 10^-6, 10^-3)*4*nAR
m109<-runif(nsims, 10^-6, 10^-3)*4*nAR

m1_78<-runif(nsims, 10^-6, 10^-3)*4*nAR
m1_87<-runif(nsims, 10^-6, 10^-3)*4*nAR
m1_810<-runif(nsims, 10^-6, 10^-3)*4*nAR
m1_108<-runif(nsims, 10^-6, 10^-3)*4*nAR


rP<-1/runif(nsims,min=2,max=100)
rEA<-1/runif(nsims,min=2,max=100)
rG<-1/runif(nsims,min=2,max=100)

tdYG<-sample(40000:145000,nsims,replace=T)
tGbot<-(tdYG-tbot)
tdOA<-samp_int_vec(35000,tGbot)
tOAbot<-(tdOA-tbot)
tdEA<-sample(20000:30000,nsims,replace=T)

taNA<-samp_int_vec(20000,tdEA)
paNA<-runif(nsims, 10^-3, 10^-1)
paNA<-1-paNA

taDP<-samp_int_vec(30000,tOAbot)
paDP<-runif(nsims, 10^-3, 10^-1)
paDP<-1-paDP

taARP<-samp_vec_vec(taDP,tOAbot)
paARP<-runif(nsims, 10^-3, 10^-1)
paARP<-1-paARP

taNEA<-samp_vec_vec(tdEA,tOAbot)
paNEA<-runif(nsims, 10^-3, 10^-1)
paNEA<-1-paNEA

taNG<-samp_vec_vec(tdOA,tGbot)
paNG<-runif(nsims, 10^-3, 10^-1)
paNG<-1-paNG

tdNNR<-110000
tdDDR<-393000
tdDN<-495000
tdADN<-580000
tdAM<-638000


#param transformations
tnAR<-4*nAR*mu*ll
snD<-nD*4*mu*ll/tnAR
snDR<-nDR*4*mu*ll/tnAR
snN<-nN*4*mu*ll/tnAR
snNR<-nNR*4*mu*ll/tnAR
snY<-nY*4*mu*ll/tnAR
snG<-nG*4*mu*ll/tnAR
snE<-nE*4*mu*ll/tnAR
snA<-nA*4*mu*ll/tnAR
snP<-nP*4*mu*ll/tnAR
snYG<-nYG*4*mu*ll/tnAR
snNNR<-nNNR*4*mu*ll/tnAR
snDDR<-nDDR*4*mu*ll/tnAR
snDN<-nDN*4*mu*ll/tnAR
snADN<-nADN*4*mu*ll/tnAR
snAM<-nAM*4*mu*ll/tnAR

stdYG<-(tdYG/tgen)/(4*nAR)
stGbot<-(tGbot/tgen)/(4*nAR)
stdOA<-(tdOA/tgen)/(4*nAR)
stOAbot<-(tOAbot/tgen)/(4*nAR)
stdEA<-(tdEA/tgen)/(4*nAR)

staNA<-(taNA/tgen)/(4*nAR)
staDP<-(taDP/tgen)/(4*nAR) 
staARP<-(taARP/tgen)/(4*nAR)
staNEA<-(taNEA/tgen)/(4*nAR)
staNG<-(taNG/tgen)/(4*nAR)

sstN<-(stN/tgen)/(4*nAR)
sstD<-(stD/tgen)/(4*nAR)

stdNNR<-(tdNNR/tgen)/(4*nAR)
stdDDR<-(tdDDR/tgen)/(4*nAR)
stdDN<-(tdDN/tgen)/(4*nAR)
stdADN<-(tdADN/tgen)/(4*nAR)
stdAM<-(tdAM/tgen)/(4*nAR)
tm1<-(tdEA/tgen)/(4*nAR)
srec<-4*nAR*(recomb*(ll-1))

#Output: parameters files
partable<-cbind(nAR,nD,nDR,nN,nNR,nY,nG,nE,nA,nP,nYG,nNNR,nDDR,nDN,nADN,nAM,rYG,rNNR,rDDR,rNDN,rADN,rAM,rP,rEA,rG,tdYG,tGbot,tdOA,tOAbot,tdEA,taNA,paNA,taDP,paDP,taARP,paARP,taNEA,paNEA,taNG,paNG,stN,stD,tdNNR,tdDDR,tdDN,tdADN,tdAM)
colnames(partable)<-c("nAR","nD","nDR","nN","nNR","nY","nG","nE","nA","nP","nYG","nNNR","nDDR","nDN","nADN","nAM","rYG","rNNR","rDDR","rNDN","rADN","rAM","rP","rEA","rG","tdYG","tGbot","tdOA","tOAbot","tdEA","taNA","paNA","taDP","paDP","taARP","paARP","taNEA","paNEA","taNG","paNG","stN","stD","tdNNR","tdDDR","tdDN","tdADN","tdAM")
partablescaled<-cbind(tnAR,snD,snDR,snN,snNR,snY,snG,snE,snA,snP,snYG,snNNR,snDDR,snDN,snADN,snAM,srec,stdYG,stGbot,stdOA,stOAbot,stdEA,staNA,staDP,staARP,staNEA,staNG,sstN,sstD,stdNNR,stdDDR,stdDN,stdADN,stdAM,m67,m76,m78,m87,m89,m98,m910,m109,tm1,m1_78,m1_87,m1_810,m1_108)
write.table(partable,paste(out,".param",sep=""),row.names=F,quote=F,sep="\t")
write.table(partablescaled,paste(out,".paramscaled",sep=""),row.names=F,col.names=T,quote=F,sep="\t")

#Output summary statistics: FDSS
i<-1
for (i in 1:nsims){
	#admixture neander > eurasia prima (forward) di admixture arcaica e admixture denisova
	li1<-paste(ms," ",as.character(6*as.numeric(nchr))," ",as.character(nloci)," -t ",as.character(tnAR[i])," -r ",as.character(srec[i])," ",as.character(ll)," -IT 3 10 ","0 0 0 0 0 0 ",nchr," 0 ",nchr," ",nchr," ",nchr," ",as.character(sstD[i])," 0 ",nchr," 0 0 0 0 0 0 0 0 ",as.character(sstN[i])," 0 0 0 ",nchr," 0 0 0 0 0 0 -n 6 ",as.character(snY[i])," -n 7 ",as.character(snG[i])," -n 8 ",as.character(snE[i])," -n 9 ",as.character(snA[i])," -n 10 ",as.character(snP[i])," -m 6 7 ",as.character(m67[i])," -m 7 6 ",as.character(m76[i])," -m 7 8 ",as.character(m78[i])," -m 8 7 ",as.character(m87[i])," -m 8 9 ",as.character(m89[i])," -m 9 8 ",as.character(m98[i])," -m 9 10 ",as.character(m910[i])," -m 10 9 ",as.character(m109[i])," -es ",as.character(staNA[i])," 5 ",as.character(paNA[i])," -ej ",as.character(staNA[i]+0.000001)," 11 9 -ej ",as.character(stdEA[i])," 9 8 -en ",as.character(stdEA[i])," 8 ",as.character(snE[i])," -em ",as.character(tm1[i])," 7 8 ",as.character(m1_78[i])," -em ",as.character(tm1[i])," 8 7 ",as.character(m1_87[i])," -em ",as.character(tm1[i])," 8 10 ",as.character(m1_810[i])," -em ",as.character(tm1[i])," 10 8 ",as.character(m1_108[i])," -es ",as.character(staDP[i])," 3 ",as.character(paDP[i])," -ej ",as.character(staDP[i]+0.000001)," 12 10 -es ",as.character(staARP[i])," 1 ",as.character(paARP[i])," -ej ",as.character(staARP[i]+0.000001)," 13 10 -es ",as.character(staNEA[i])," 5 ",as.character(paNEA[i])," -ej ",as.character(staNEA[i]+0.000001)," 14 8 -en ",as.character(stOAbot[i])," 8 ",as.character(rEA[i])," -en ",as.character(stOAbot[i])," 10 ",as.character(rP[i])," -ej ",as.character(stdOA[i])," 8 7 -ej ",as.character(stdOA[i])," 10 7 -es ",as.character(staNG[i])," 5 ",as.character(paNG[i])," -ej ",as.character(staNG[i]+0.000001)," 15 7 -ej ",as.character(stdNNR[i])," 5 4 -en ",as.character(stdNNR[i])," 4 ",as.character(snNNR[i])," -en ",as.character(stGbot[i])," 7 ",as.character(rG[i])," -ej ",as.character(stdYG[i])," 7 6 -en ",as.character(stdYG[i])," 6 ",as.character(snYG[i])," -ej ",as.character(stdDDR[i])," 3 2 -en ",as.character(stdDDR[i])," 2 ",as.character(snDDR[i])," -ej ",as.character(stdDN[i])," 4 2 -en ",as.character(stdDN[i])," 2 ",as.character(snDN[i])," -ej ",as.character(stdADN[i])," 2 1 -en ",as.character(stdADN[i])," 1 ",as.character(snDN[i])," -ej ",as.character(stdAM[i])," 6 1 -en ",as.character(stdAM[i])," 1 ",as.character(snAM[i]),sep="")
	#admixture neander > eurasia dopo (forward) di admixture arcaica e admixture denisova
	li2<-paste(ms," ",as.character(6*as.numeric(nchr))," ",as.character(nloci)," -t ",as.character(tnAR[i])," -r ",as.character(srec[i])," ",as.character(ll)," -IT 3 10 ","0 0 0 0 0 0 ",nchr," 0 ",nchr," ",nchr," ",nchr," ",as.character(sstD[i])," 0 ",nchr," 0 0 0 0 0 0 0 0 ",as.character(sstN[i])," 0 0 0 ",nchr," 0 0 0 0 0 0 -n 6 ",as.character(snY[i])," -n 7 ",as.character(snG[i])," -n 8 ",as.character(snE[i])," -n 9 ",as.character(snA[i])," -n 10 ",as.character(snP[i])," -m 6 7 ",as.character(m67[i])," -m 7 6 ",as.character(m76[i])," -m 7 8 ",as.character(m78[i])," -m 8 7 ",as.character(m87[i])," -m 8 9 ",as.character(m89[i])," -m 9 8 ",as.character(m98[i])," -m 9 10 ",as.character(m910[i])," -m 10 9 ",as.character(m109[i])," -es ",as.character(staNA[i])," 5 ",as.character(paNA[i])," -ej ",as.character(staNA[i]+0.000001)," 11 9 -ej ",as.character(stdEA[i])," 9 8 -en ",as.character(stdEA[i])," 8 ",as.character(snE[i])," -em ",as.character(tm1[i])," 7 8 ",as.character(m1_78[i])," -em ",as.character(tm1[i])," 8 7 ",as.character(m1_87[i])," -em ",as.character(tm1[i])," 8 10 ",as.character(m1_810[i])," -em ",as.character(tm1[i])," 10 8 ",as.character(m1_108[i])," -es ",as.character(staNEA[i])," 5 ",as.character(paNEA[i])," -ej ",as.character(staNEA[i]+0.000001)," 12 8 -es ",as.character(staDP[i])," 3 ",as.character(paDP[i])," -ej ",as.character(staDP[i]+0.000001)," 13 10 -es ",as.character(staARP[i])," 1 ",as.character(paARP[i])," -ej ",as.character(staARP[i]+0.000001)," 14 10 -en ",as.character(stOAbot[i])," 8 ",as.character(rEA[i])," -en ",as.character(stOAbot[i])," 10 ",as.character(rP[i])," -ej ",as.character(stdOA[i])," 8 7 -ej ",as.character(stdOA[i])," 10 7 -es ",as.character(staNG[i])," 5 ",as.character(paNG[i])," -ej ",as.character(staNG[i]+0.000001)," 15 7 -ej ",as.character(stdNNR[i])," 5 4 -en ",as.character(stdNNR[i])," 4 ",as.character(snNNR[i])," -en ",as.character(stGbot[i])," 7 ",as.character(rG[i])," -ej ",as.character(stdYG[i])," 7 6 -en ",as.character(stdYG[i])," 6 ",as.character(snYG[i])," -ej ",as.character(stdDDR[i])," 3 2 -en ",as.character(stdDDR[i])," 2 ",as.character(snDDR[i])," -ej ",as.character(stdDN[i])," 4 2 -en ",as.character(stdDN[i])," 2 ",as.character(snDN[i])," -ej ",as.character(stdADN[i])," 2 1 -en ",as.character(stdADN[i])," 1 ",as.character(snDN[i])," -ej ",as.character(stdAM[i])," 6 1 -en ",as.character(stdAM[i])," 1 ",as.character(snAM[i]),sep="")
	#admixture neander > eurasia dopo (forward) di admixture arcaica e prima admixture admixture denisova
	li3<-paste(ms," ",as.character(6*as.numeric(nchr))," ",as.character(nloci)," -t ",as.character(tnAR[i])," -r ",as.character(srec[i])," ",as.character(ll)," -IT 3 10 ","0 0 0 0 0 0 ",nchr," 0 ",nchr," ",nchr," ",nchr," ",as.character(sstD[i])," 0 ",nchr," 0 0 0 0 0 0 0 0 ",as.character(sstN[i])," 0 0 0 ",nchr," 0 0 0 0 0 0 -n 6 ",as.character(snY[i])," -n 7 ",as.character(snG[i])," -n 8 ",as.character(snE[i])," -n 9 ",as.character(snA[i])," -n 10 ",as.character(snP[i])," -m 6 7 ",as.character(m67[i])," -m 7 6 ",as.character(m76[i])," -m 7 8 ",as.character(m78[i])," -m 8 7 ",as.character(m87[i])," -m 8 9 ",as.character(m89[i])," -m 9 8 ",as.character(m98[i])," -m 9 10 ",as.character(m910[i])," -m 10 9 ",as.character(m109[i])," -es ",as.character(staNA[i])," 5 ",as.character(paNA[i])," -ej ",as.character(staNA[i]+0.000001)," 11 9 -ej ",as.character(stdEA[i])," 9 8 -en ",as.character(stdEA[i])," 8 ",as.character(snE[i])," -em ",as.character(tm1[i])," 7 8 ",as.character(m1_78[i])," -em ",as.character(tm1[i])," 8 7 ",as.character(m1_87[i])," -em ",as.character(tm1[i])," 8 10 ",as.character(m1_810[i])," -em ",as.character(tm1[i])," 10 8 ",as.character(m1_108[i])," -es ",as.character(staDP[i])," 3 ",as.character(paDP[i])," -ej ",as.character(staDP[i]+0.000001)," 12 10 -es ",as.character(staNEA[i])," 5 ",as.character(paNEA[i])," -ej ",as.character(staNEA[i]+0.000001)," 13 8 -es ",as.character(staARP[i])," 1 ",as.character(paARP[i])," -ej ",as.character(staARP[i]+0.000001)," 14 10 -en ",as.character(stOAbot[i])," 8 ",as.character(rEA[i])," -en ",as.character(stOAbot[i])," 10 ",as.character(rP[i])," -ej ",as.character(stdOA[i])," 8 7 -ej ",as.character(stdOA[i])," 10 7 -es ",as.character(staNG[i])," 5 ",as.character(paNG[i])," -ej ",as.character(staNG[i]+0.000001)," 15 7 -ej ",as.character(stdNNR[i])," 5 4 -en ",as.character(stdNNR[i])," 4 ",as.character(snNNR[i])," -en ",as.character(stGbot[i])," 7 ",as.character(rG[i])," -ej ",as.character(stdYG[i])," 7 6 -en ",as.character(stdYG[i])," 6 ",as.character(snYG[i])," -ej ",as.character(stdDDR[i])," 3 2 -en ",as.character(stdDDR[i])," 2 ",as.character(snDDR[i])," -ej ",as.character(stdDN[i])," 4 2 -en ",as.character(stdDN[i])," 2 ",as.character(snDN[i])," -ej ",as.character(stdADN[i])," 2 1 -en ",as.character(stdADN[i])," 1 ",as.character(snDN[i])," -ej ",as.character(stdAM[i])," 6 1 -en ",as.character(stdAM[i])," 1 ",as.character(snAM[i]),sep="")
	print(i)
	#print (c(staNEA[i],staDP[i],staARP[i]))
	if (i==1){
		if (staNEA[i] >= staARP[i] ){
			#print("uno")
			#write(li1,stderr())
			system(paste(li1," | ",cpd,"compute_ss.py -np 6 -nc ",nchr," -w 30 -b 50 -s > ",out,".tab",sep=""))
		}else if (staNEA[i] >= staDP[i] & staNEA[i] <= staARP[i]){
			#print("tre")
			#write(li3,stderr())
			system(paste(li3," | ",cpd,"compute_ss.py -np 6 -nc ",nchr," -w 30 -b 50 -s > ",out,".tab",sep=""))
		}else if (staNEA[i] < staDP[i]) {
			#print("due")
			#write(li2,stderr())
			system(paste(li2," | ",cpd,"compute_ss.py -np 6 -nc ",nchr," -w 30 -b 50 -s > ",out,".tab",sep=""))
		}
	}else{
		if (staNEA[i] >= staARP[i] ){
			#print("uno")
			#write(li1,stderr())
			system(paste(li1," | ",cpd,"compute_ss.py -np 6 -nc ",nchr," -w 30 -b 50 -s >> ",out,".tab",sep=""))
		}else if (staNEA[i] >= staDP[i] & staNEA[i] <= staARP[i]){
			#print("tre")
			#write(li3,stderr())
			system(paste(li3," | ",cpd,"compute_ss.py -np 6 -nc ",nchr," -w 30 -b 50 -s >> ",out,".tab",sep=""))
		}else if (staNEA[i] < staDP[i]) {
			#print("due")
			#write(li2,stderr())
			system(paste(li2," | ",cpd,"compute_ss.py -np 6 -nc ",nchr," -w 30 -b 50 -s >> ",out,".tab",sep=""))
		}
	}
}
