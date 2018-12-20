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


#wrapper script to simulate two population split with asymmetrical gene flow from divergence to the present
#Rscript --vanilla ./ooa_m2.r <nchr> <locuslength> <nlociXsim> <recrate=1.12e-8> <nsimulations>
args<-commandArgs(trailingOnly=TRUE)
ms<-"/opt/software/genetics/msms/bin/msms"
cpd<-"./"
mod<-"ooa_m2"
nchr<-as.character(args[1])
tgen<-29
mu<-1.25e-8
recomb<-as.numeric(args[4])
ll<-as.numeric(args[2])#locus length
nsims<-as.numeric(args[5])#number of ABC simulations
nloci<-as.numeric(args[3])#loci to simulate in each sim
out<-paste(mod,"_ll",as.character(ll),"_nl",as.character(nloci),"_r",as.character(recomb),"_nc",nchr,sep="")
#main param
tbot<-2900
stN<-85735
stD<-67570
nAR<-sample(50:50000,nsims,replace=T)
nD<-sample(50:50000,nsims,replace=T)
nDR<-sample(50:50000,nsims,replace=T)
nN<-sample(50:50000,nsims,replace=T)
nNR<-sample(50:50000,nsims,replace=T)
nY<-sample(50:50000,nsims,replace=T)
nG1<-sample(50:50000,nsims,replace=T)
nG2<-sample(50:50000,nsims,replace=T)
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
m68<-runif(nsims, 10^-6, 10^-3)*4*nAR
m86<-runif(nsims, 10^-6, 10^-3)*4*nAR
m89<-runif(nsims, 10^-6, 10^-3)*4*nAR
m98<-runif(nsims, 10^-6, 10^-3)*4*nAR
m910<-runif(nsims, 10^-6, 10^-3)*4*nAR
m109<-runif(nsims, 10^-6, 10^-3)*4*nAR
m1011<-runif(nsims, 10^-6, 10^-3)*4*nAR
m1110<-runif(nsims, 10^-6, 10^-3)*4*nAR

m1_89<-runif(nsims, 10^-6, 10^-3)*4*nAR
m1_98<-runif(nsims, 10^-6, 10^-3)*4*nAR
m1_911<-runif(nsims, 10^-6, 10^-3)*4*nAR
m1_119<-runif(nsims, 10^-6, 10^-3)*4*nAR


rP<-1/runif(nsims,min=2,max=100)
rEA<-1/runif(nsims,min=2,max=100)

tdYG1<-sample(40000:145000,nsims,replace=T)
tdYG2<-tdYG1
#tdOA1<-sample(35000:tdYG1,nsims,replace=T)
tdOA1<-samp_int_vec(40000,tdYG1)
tOAbot1<-(tdOA1-tbot)
#tdOA2<-sample(30000:tdOA1,nsims,replace=T)
tdOA2<-samp_int_vec(35000,tOAbot1)                      
tOAbot2<-(tdOA2-tbot)
#tdEA<-sample(20000:tOAbot2,nsims,replace=T)
tdEA<-samp_int_vec(20000,tOAbot2)

#taNA<-sample(20000:tdEA,nsims,replace=T)
taNA<-samp_int_vec(20000,tdEA)
paNA<-runif(nsims, 10^-3, 10^-1)
paNA<-1-paNA
#taDP<-sample(30000:tOAbot1,nsims,replace=T) 
taDP<-samp_int_vec(30000,tOAbot1)
paDP<-runif(nsims, 10^-3, 10^-1)
paDP<-1-paDP
#taARP<-sample(taDP:tOAbot1,nsims,replace=T)
taARP<-samp_vec_vec(taDP,tOAbot1)
paARP<-runif(nsims, 10^-3, 10^-1)
paARP<-1-paARP
#taNEA<-sample(tdEA:tOAbot2,nsims,replace=T)
taNEA<-samp_vec_vec(tdEA,tOAbot2)
paNEA<-runif(nsims, 10^-3, 10^-1)
paNEA<-1-paNEA
#taNG2<-sample(tdOA2:tdYG2,nsims,replace=T)
taNG2<-samp_vec_vec(tdOA2,tdYG2)
paNG2<-runif(nsims, 10^-3, 10^-1)
paNG2<-1-paNG2

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
snG1<-nG1*4*mu*ll/tnAR
snG2<-nG2*4*mu*ll/tnAR
snE<-nE*4*mu*ll/tnAR
snA<-nA*4*mu*ll/tnAR
snP<-nP*4*mu*ll/tnAR
snYG<-nYG*4*mu*ll/tnAR
snNNR<-nNNR*4*mu*ll/tnAR
snDDR<-nDDR*4*mu*ll/tnAR
snDN<-nDN*4*mu*ll/tnAR
snADN<-nADN*4*mu*ll/tnAR
snAM<-nAM*4*mu*ll/tnAR


stdYG1<-(tdYG1/tgen)/(4*nAR)
stdYG2<-(tdYG2/tgen)/(4*nAR)
stdOA1<-(tdOA1/tgen)/(4*nAR)
stdOA2<-(tdOA2/tgen)/(4*nAR)
stOAbot1<-(tOAbot1/tgen)/(4*nAR)
stOAbot2<-(tOAbot2/tgen)/(4*nAR)
stdEA<-(tdEA/tgen)/(4*nAR)

staNA<-(taNA/tgen)/(4*nAR)
staDP<-(taDP/tgen)/(4*nAR) 
staARP<-(taARP/tgen)/(4*nAR)
staNEA<-(taNEA/tgen)/(4*nAR)
staNG2<-(taNG2/tgen)/(4*nAR)

sstN<-(stN/tgen)/(4*nAR)
sstD<-(stD/tgen)/(4*nAR)

stdNNR<-(tdNNR/tgen)/(4*nAR)
stdDDR<-(tdDDR/tgen)/(4*nAR)
stdDN<-(tdDN/tgen)/(4*nAR)
stdADN<-(tdADN/tgen)/(4*nAR)
stdAM<-(tdAM/tgen)/(4*nAR)

tm1<-(tdEA/tgen)/(4*nAR)
srec<-4*nAR*(recomb*(ll-1))

partable<-cbind(nAR,nD,nDR,nN,nNR,nY,nG1,nG2,nE,nA,nP,nYG,nNNR,nDDR,nDN,nADN,nAM,rYG,rNNR,rDDR,rNDN,rADN,rAM,rP,rEA,tdYG1,tdYG2,tdOA1,tdOA2,tOAbot1,tOAbot2,tdEA,taNA,paNA,taDP,paDP,taARP,paARP,taNEA,paNEA,taNG2,paNG2,stN,stD,tdNNR,tdDDR,tdDN,tdADN,tdAM)
colnames(partable)<-c("nAR","nD","nDR","nN","nNR","nY","nG1","nG2","nE","nA","nP","nYG","nNNR","nDDR","nDN","nADN","nAM","rYG","rNNR","rDDR","rNDN","rADN","rAM","rP","rEA","tdYG1","tdYG2","tdOA1","tdOA2","tOAbot1","tOAbot2","tdEA","taNA","paNA","taDP","paDP","taARP","paARP","taNEA","paNEA","taNG2","paNG2","stN","stD","tdNNR","tdDDR","tdDN","tdADN","tdAM")
partablescaled<-cbind(tnAR,snD,snDR,snN,snNR,snY,snG1,snG2,snE,snA,snP,snYG,snNNR,snDDR,snDN,snADN,snAM,srec,stdYG1,stdYG2,stdOA1,stdOA2,stOAbot1,stOAbot2,stdEA,staNA,staDP,staARP,staNEA,staNG2,sstN,sstD,stdNNR,stdDDR,stdDN,stdADN,stdAM,m67,m76,m78,m87,m89,m98,m910,m109,m1011,m1110,tm1,m1_89,m1_98,m1_911,m1_119)
write.table(partable,paste(out,".param",sep=""),row.names=F,quote=F,sep="\t")
write.table(partablescaled,paste(out,".paramscaled",sep=""),row.names=F,col.names=T,quote=F,sep="\t")
i<-1
for (i in 1:nsims){	
	s<-c()
	s[1]<-paste(" -es ",as.character(staNA[i])," 5 ",as.character(paNA[i])," -ej ",as.character(staNA[i]+0.000001)," IDPOP 10 ",sep="")
	s[2]<-paste(" -ej ",as.character(stdEA[i])," 10 9",sep="") 
	s[3]<-paste(" -en ",as.character(stdEA[i])," 9 ",as.character(snE[i]),sep="")
	s[4]<-paste(" -em ",as.character(tm1[i])," 8 9 ",as.character(m1_89[i]),sep="")
	s[5]<-paste(" -em ",as.character(tm1[i])," 9 8 ",as.character(m1_98[i]),sep="")
	s[6]<-paste(" -em ",as.character(tm1[i])," 9 11 ",as.character(m1_911[i]),sep="")
	s[7]<-paste(" -em ",as.character(tm1[i])," 11 9 ",as.character(m1_119[i]),sep="")
	s[8]<-paste(" -es ",as.character(staDP[i])," 3 ",as.character(paDP[i])," -ej ",as.character(staDP[i]+0.000001)," IDPOP 11 ",sep="") 
	s[9]<-paste(" -es ",as.character(staARP[i])," 1 ",as.character(paARP[i])," -ej ",as.character(staARP[i]+0.000001)," IDPOP 11 ",sep="") 
	s[10]<-paste(" -es ",as.character(staNEA[i])," 5 ",as.character(paNEA[i])," -ej ",as.character(staNEA[i]+0.000001)," IDPOP 9",sep="")
	s[11]<-paste(" -en ",as.character(stOAbot2[i])," 9 ",as.character(rEA[i]),sep="") 
	s[12]<-paste(" -ej ",as.character(stdOA2[i])," 9 8 " ,sep="")  
	s[13]<-paste(" -es ",as.character(staNG2[i])," 5 ",as.character(paNG2[i])," -ej ",as.character(staNG2[i]+0.000001)," IDPOP 8 ",sep="")   
	s[14]<-paste(" -en ",as.character(stOAbot1[i])," 11 ",as.character(rP[i]),sep="")
	s[15]<-paste(" -ej ",as.character(stdOA1[i])," 11 7",sep="") 
	s[16]<-paste(" -ej ",as.character(stdNNR[i])," 5 4",sep="") 
	s[17]<-paste(" -en ",as.character(stdNNR[i])," 4 ",as.character(snNNR[i]),sep="") 
	s[18]<-paste(" -ej ",as.character(stdYG2[i])," 8 6 ",sep="")
	s[19]<-paste(" -ej ",as.character(stdYG1[i])," 7 6",sep="")
	s[20]<-paste(" -en ",as.character(stdYG1[i])," 6 ",as.character(snYG[i]),sep="") 
	s[21]<-paste(" -ej ",as.character(stdDDR[i])," 3 2",sep="") 
	s[22]<-paste(" -en ",as.character(stdDDR[i])," 2 ",as.character(snDDR[i]),sep="") 
	s[23]<-paste(" -ej ",as.character(stdDN[i])," 4 2 ",sep="")
	s[24]<-paste(" -en ",as.character(stdDN[i])," 2 ",as.character(snDN[i]),sep="") 
	s[25]<-paste(" -ej ",as.character(stdADN[i])," 2 1 ",sep="")
	s[26]<-paste(" -en ",as.character(stdADN[i])," 1 ",as.character(snDN[i]),sep="") 
	s[27]<-paste(" -ej ",as.character(stdAM[i])," 6 1 ",sep="")
	s[28]<-paste(" -en ",as.character(stdAM[i])," 1 ",as.character(snAM[i]),sep="")
	s1<-c()
	s1[1]<-staNA[i]
	s1[2]<-stdEA[i]
	s1[3]<-stdEA[i]
	s1[4]<-tm1[i]
	s1[5]<-tm1[i]
	s1[6]<-tm1[i]
	s1[7]<-tm1[i]
	s1[8]<-staDP[i]
	s1[9]<-staARP[i]
	s1[10]<-staNEA[i]
	s1[11]<-stOAbot2[i]
	s1[12]<-stdOA2[i]
	s1[13]<-staNG2[i]
	s1[14]<-stOAbot1[i]
	s1[15]<-stdOA1[i]
	s1[16]<-stdNNR[i]
	s1[17]<-stdNNR[i]
	s1[18]<-stdYG2[i]
	s1[19]<-stdYG1[i]
	s1[20]<-stdYG1[i]
	s1[21]<-stdDDR[i]
	s1[22]<-stdDDR[i]
	s1[23]<-stdDN[i]
	s1[24]<-stdDN[i]
	s1[25]<-stdADN[i]
	s1[26]<-stdADN[i]
	s1[27]<-stdAM[i]
	s1[28]<-stdAM[i]

	sid<-sort(s1,index.return=T)
	s_sort<-s[sid$ix]
	ig<-grep("IDPOP",s_sort)
	y<-11
	for (k in ig){
		y<-y+1
		s_sort[k]<-sub("IDPOP",as.character(y),s_sort[k])
	}

	part1<-paste(s_sort,collapse="")

	li1<-paste(ms," ",as.character(6*as.numeric(nchr))," ",as.character(nloci)," -t ",as.character(tnAR[i])," -r ",as.character(srec[i])," ",as.character(ll)," -IT 3 11 ","0 0 0 0 0 0 ",nchr," 0 0 ",nchr," ",nchr," ",nchr," ",as.character(sstD[i])," 0 ",nchr," 0 0 0 0 0 0 0 0 0 ",as.character(sstN[i])," 0 0 0 ",nchr," 0 0 0 0 0 0 0 -n 6 ",as.character(snY[i])," -n 7 ",as.character(snG1[i])," -n 8 ",as.character(snG2[i])," -n 9 ",as.character(snE[i])," -n 10 ",as.character(snA[i])," -n 11 ",as.character(snP[i])," -m 6 7 ",as.character(m67[i])," -m 7 6 ",as.character(m76[i])," -m 7 8 ",as.character(m78[i])," -m 8 7 ",as.character(m87[i])," -m 6 8 ",as.character(m68[i])," -m 8 6 ",as.character(m86[i])," -m 8 9 ",as.character(m89[i])," -m 9 8 ",as.character(m98[i])," -m 9 10 ",as.character(m910[i])," -m 10 9 ",as.character(m109[i])," -m 10 11 ",as.character(m1011[i])," -m 11 10 ",as.character(m1110[i]),part1,sep="")
	#print (c(staNEA[i],staDP[i],staARP[i]))
	if (i==1){
		system(paste(li1," | ",cpd,"compute_ss.py -np 6 -nc ",nchr," -w 30 -b 50 -s > ",out,".tab",sep=""))
	}else{
		system(paste(li1," | ",cpd,"compute_ss.py -np 6 -nc ",nchr," -w 30 -b 50 -s >> ",out,".tab",sep=""))
	}
}
