#### SCRIPTS - SIMULATIONS MULTIPLE DISPERSAL MODEL

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


#Rscript --vanilla ./ooa_m2.r <nchr> <locuslength> <nlociXsim> <recrate=1.12e-8> <nsimulations>

args<-commandArgs(trailingOnly=TRUE)
ms<-"/software/msdir/ms"
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
tbot<-2900 #Duration Bottleneck
stN<-85735 #Sample time Neandertal
stD<-67570 #Sample time Denisova

#Effective popultion sizes
nAR<-sample(500:50000,nsims,replace=T)
nD<-sample(500:50000,nsims,replace=T)
nD1<-sample(500:50000,nsims,replace=T)
nD2<-sample(500:50000,nsims,replace=T)
nN<-sample(500:50000,nsims,replace=T)
nNR<-sample(500:50000,nsims,replace=T)
nY<-sample(500:50000,nsims,replace=T)
nG1<-sample(500:50000,nsims,replace=T)
nBE<-sample(500:50000,nsims,replace=T)
nG2<-sample(500:50000,nsims,replace=T)
nE<-sample(500:50000,nsims,replace=T)
nA<-sample(500:50000,nsims,replace=T)
nP<-sample(500:50000,nsims,replace=T)
nYG<-sample(500:50000,nsims,replace=T)
nNNR<-sample(500:50000,nsims,replace=T)
nDDR<-sample(500:50000,nsims,replace=T)
nDN<-sample(500:50000,nsims,replace=T)
nADN<-sample(500:50000,nsims,replace=T)
nAM<-sample(500:50000,nsims,replace=T)

#Migration modern populations
m78<-runif(nsims, 10^-6, 10^-3)*4*nAR
m87<-runif(nsims, 10^-6, 10^-3)*4*nAR
m810<-runif(nsims, 10^-6, 10^-3)*4*nAR
m108<-runif(nsims, 10^-6, 10^-3)*4*nAR
m1011<-runif(nsims, 10^-6, 10^-3)*4*nAR
m1110<-runif(nsims, 10^-6, 10^-3)*4*nAR
m1112<-runif(nsims, 10^-6, 10^-3)*4*nAR
m1211<-runif(nsims, 10^-6, 10^-3)*4*nAR
m1213<-runif(nsims, 10^-6, 10^-3)*4*nAR
m1312<-runif(nsims, 10^-6, 10^-3)*4*nAR

#Migrations Eurasia
m1_1011<-runif(nsims, 10^-6, 10^-3)*4*nAR
m1_1110<-runif(nsims, 10^-6, 10^-3)*4*nAR
m1_1113<-runif(nsims, 10^-6, 10^-3)*4*nAR
m1_1311<-runif(nsims, 10^-6, 10^-3)*4*nAR

#Bottleneck intensity
rP<-1/runif(nsims,min=2,max=100) #Papua
rEA<-1/runif(nsims,min=2,max=100) #Eurasia

#Divergence time Ancient lineages (Neanderthal, Denisova and Archaic)
tdNNR<-110000
tdDDR<-393000
tdDD1<-tdDDR
tdDD2<-tdDDR
tdDN<-495000
tdADN<-580000
tdAM<-638000

#Divergence time modern lineages
tdYG<-sample(50000:145000,nsims,replace=T)
tdYG1<-tdYG
tdYG2<-tdYG

tdOA1<-samp_int_vec(45000,tdYG1)
tOAbot1<-(tdOA1-tbot)

tdG2BE<-samp_int_vec(50000,tdYG2)
tmin<-cbind()
i<-1
for (i in 1:nsims){
  tmin[i]<-min(tdOA1[i],tdG2BE[i],tdNNR-tgen)}

tdOA2<-samp_int_vec(40000,tmin)
tOAbot2<-(tdOA2-tbot)

tdEA<-samp_int_vec(30000,tOAbot2)

#Time and Rate Admixture events
#Denisova-Asia
taD2A<-samp_int_vec(20000,tdEA)
paD2A<-runif(nsims, 10^-3, 10^-1)
paD2A<-1-paD2A

#Basal Europe-Europe
taBEE<-samp_int_vec(10000,tdEA)
paBEE<-runif(nsims, 0.05, 0.5)
paBEE<-1-paBEE

#Denisova-Papua
taD1P<-samp_int_vec(30000,tOAbot1)
paD1P<-runif(nsims, 10^-3, 10^-1)
paD1P<-1-paD1P

#Archaic Papua
taARP<-samp_vec_vec(taD1P,tOAbot1)
paARP<-runif(nsims, 10^-3, 10^-1)
paARP<-1-paARP

#Neanderthal-Asia
taNEA<-samp_vec_vec(tdEA,tOAbot2)
paNEA<-runif(nsims, 10^-3, 10^-1)
paNEA<-1-paNEA

#Neanderthal-African Ghost
taNG2<-samp_vec_vec(tdOA2,tmin)
paNG2<-runif(nsims, 10^-3, 10^-1)
paNG2<-1-paNG2

#Parameters Scaled for ms
tnAR<-4*nAR*mu*ll #theta

snD<-nD*4*mu*ll/tnAR
snD1<-nD1*4*mu*ll/tnAR
snD2<-nD2*4*mu*ll/tnAR
snN<-nN*4*mu*ll/tnAR
snNR<-nNR*4*mu*ll/tnAR
snY<-nY*4*mu*ll/tnAR
snG1<-nG1*4*mu*ll/tnAR
snG2<-nG2*4*mu*ll/tnAR
snBE<-nBE*4*mu*ll/tnAR
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
stdYG1<-(tdYG1/tgen)/(4*nAR)
stdYG2<-(tdYG2/tgen)/(4*nAR)
stdG2BE<-(tdG2BE/tgen)/(4*nAR)
stdOA1<-(tdOA1/tgen)/(4*nAR)
stdOA2<-(tdOA2/tgen)/(4*nAR)
stOAbot1<-(tOAbot1/tgen)/(4*nAR)
stOAbot2<-(tOAbot2/tgen)/(4*nAR)
stdEA<-(tdEA/tgen)/(4*nAR)

staD2A<-(taD2A/tgen)/(4*nAR)
staBEE<-(taBEE/tgen)/(4*nAR)
staD1P<-(taD1P/tgen)/(4*nAR) 
staARP<-(taARP/tgen)/(4*nAR)
staNEA<-(taNEA/tgen)/(4*nAR)
staNG2<-(taNG2/tgen)/(4*nAR)

sstN<-(stN/tgen)/(4*nAR)
sstD<-(stD/tgen)/(4*nAR)

stdNNR<-(tdNNR/tgen)/(4*nAR)
stdDDR<-(tdDDR/tgen)/(4*nAR)
stdDD1<-(tdDD1/tgen)/(4*nAR)
stdDD2<-(tdDD2/tgen)/(4*nAR)
stdDN<-(tdDN/tgen)/(4*nAR)
stdADN<-(tdADN/tgen)/(4*nAR)
stdAM<-(tdAM/tgen)/(4*nAR)
tm1<-(tdEA/tgen)/(4*nAR)

srec<-4*nAR*(recomb*(ll-1))

#### OUTPUT PARAMETER's FILES
partable<-cbind(nAR,nD,nD1,nD2,nN,nNR,nY,nG1,nG2,nBE,nE,nA,nP,
                nYG,nNNR,nDDR,nDN,nADN,nAM,
                rYG,rNNR,rDDR,rNDN,rADN,rAM,rP,rEA,
                tdYG1,tdYG2,tdOA1,tOAbot1,tdOA2,tOAbot2,tdG2BE,tdEA,
                taD2A,paD2A,taBEE,paBEE,taD1P,paD1P,taARP,paARP,taNEA,paNEA,taNG2,paNG2,
                stN,stD,tdNNR,tdDDR,tdDN,tdADN,tdAM)
partablescaled<-cbind(tnAR,snD,snD1,snD2,snN,snNR,snY,snG1,snG2,snBE,snE,snA,snP,
                      snYG,snNNR,snDDR,snDN,snADN,snAM,srec,
                      stdYG1,stdYG2,stdOA1,stOAbot1,stdOA2,stOAbot2,stdG2BE,stdEA,
                      staD2A,staBEE,staD1P,staARP,staNEA,staNG2,
                      sstN,sstD,stdNNR,stdDDR,stdDN,stdADN,stdAM,
                      m78,m87,m810,m108,m1011,m1110,m1112,m1211,m1213,m1312,tm1,m1_1011,m1_1110,m1_1113,m1_1311)
write.table(partable,paste(out,".param",sep=""),row.names=F,col.names=T,quote=F,sep="\t")
write.table(partablescaled,paste(out,".paramscaled",sep=""),row.names=F,col.names=T,quote=F,sep="\t")

#### Creation ms command line:
i<-1
for (i in 1:nsims){
	s0<-c()
	#Modern populations
	s0[1]<-paste(" -I 13 0 0 0 0 0 0 ",nchr," 0 0 0 ",nchr," ",nchr," ",nchr,sep="")
	#Modern effective population size
	s0[2]<-paste(" -n 7 ",as.character(snY[i]),sep="")
	s0[3]<-paste(" -n 8 ",as.character(snG1[i]),sep="")
	s0[4]<-paste(" -n 9 ",as.character(snBE[i]),sep="")
	s0[5]<-paste(" -n 10 ",as.character(snG2[i]),sep="")
	s0[6]<-paste(" -n 11 ",as.character(snE[i]),sep="")
	s0[7]<-paste(" -n 12 ",as.character(snA[i]),sep="")
	s0[8]<-paste(" -n 13 ",as.character(snP[i]),sep="")
	#Migration between modern populations
	s0[9]<-paste(" -m 7 8 ",as.character(m78[i]),sep="")
	s0[10]<-paste(" -m 8 7 ",as.character(m87[i]),sep="")
	s0[11]<-paste(" -m 8 10 ",as.character(m810[i]),sep="")
	s0[12]<-paste(" -m 10 8 ",as.character(m108[i]),sep="")
	s0[13]<-paste(" -m 10 11 ",as.character(m1011[i]),sep="")
	s0[14]<-paste(" -m 11 10 ",as.character(m1110[i]),sep="")
	s0[15]<-paste(" -m 11 12 ",as.character(m1112[i]),sep="")
	s0[16]<-paste(" -m 12 11 ",as.character(m1211[i]),sep="")
	s0[17]<-paste(" -m 12 13 ",as.character(m1213[i]),sep="")
	s0[18]<-paste(" -m 13 12 ",as.character(m1312[i]),sep="")
	#Sampling time Denisova
	s0[19]<-paste(" -eA ",as.character(sstD[i])," 2 ",nchr,sep="")
	#Sampling time Neanderthal
	s0[20]<-paste(" -eA ",as.character(sstN[i])," 5 ",nchr,sep="")
  
	s<-c()
	s1<-c()
	#Admixtures BasalEurope - Europe
	s[1]<-paste(" -es ",as.character(staBEE[i])," 9 ",as.character(paBEE[i])," -ej ",as.character(staBEE[i]+0.000001)," IDPOP 11",sep="")
	s1[1]<-staBEE[i]
	#Admixtures Denisova Related 2 - Asia
	s[2]<-paste(" -es ",as.character(staD2A[i])," 4 ",as.character(paD2A[i])," -ej ",as.character(staD2A[i]+0.000001)," IDPOP 12",sep="")
	s1[2]<-staD2A[i]
	#Divergence Europe-Asia
	s[3]<-paste(" -ej ",as.character(stdEA[i])," 12 11",sep="")
	s1[3]<-stdEA[i]
	#Effective pop. size Eurasia
	s[4]<-paste(" -en ",as.character(stdEA[i])," 11 ",as.character(snE[i]),sep="")
	s1[4]<-stdEA[i]
	#Migrations between Eurasia-Africa and Eurasia-Papua
	s[5]<-paste(" -em ",as.character(tm1[i])," 10 11 ",as.character(m1_1011[i]),sep="")
	s1[5]<-tm1[i]
	s[6]<-paste(" -em ",as.character(tm1[i])," 11 10 ",as.character(m1_1110[i]),sep="")
	s1[6]<-tm1[i]
	s[7]<-paste(" -em ",as.character(tm1[i])," 11 13 ",as.character(m1_1113[i]),sep="")
	s1[7]<-tm1[i]
	s[8]<-paste(" -em ",as.character(tm1[i])," 13 11 ",as.character(m1_1311[i]),sep="")
	s1[8]<-tm1[i]
	#Admixtures Neanderthal Related - Eurasia
	s[9]<-paste(" -es ",as.character(staNEA[i])," 6 ",as.character(paNEA[i])," -ej ",as.character(staNEA[i]+0.000001)," IDPOP 11",sep="")
	s1[9]<-staNEA[i]
	#Admixtures Denisova Related 1 - Papua
	s[10]<-paste(" -es ",as.character(staD1P[i])," 3 ",as.character(paD1P[i])," -ej ",as.character(staD1P[i]+0.000001)," IDPOP 13",sep="")
	s1[10]<-staD1P[i]
	#Admixtures Archaic - Papua
	s[11]<-paste(" -es ",as.character(staARP[i])," 1 ",as.character(paARP[i])," -ej ",as.character(staARP[i]+0.000001)," IDPOP 13",sep="")
	s1[11]<-staARP[i]
	#Effective pop. size Eurasia during bottleneck
	s[12]<-paste(" -en ",as.character(stOAbot2[i])," 11 ",as.character(rEA[i]),sep="")
	s1[12]<-stOAbot2[i]
	#Effective pop. size Papua during bottleneck
	s[13]<-paste(" -en ",as.character(stOAbot1[i])," 13 ",as.character(rP[i]),sep="")
	s1[13]<-stOAbot1[i]
	#Divergence Ghost 2 - Eurasia (OOA 2)
	s[14]<-paste(" -ej ",as.character(stdOA2[i])," 11 10",sep="")
	s1[14]<-stdOA2[i]
	#Effective pop. size Ghost 2 + Eurasia
	s[15]<-paste(" -en ",as.character(stdOA2[i])," 10 ",as.character(snG2[i]),sep="")
	s1[15]<-stdOA2[i]
	#Admixtures Neanderthal Related - Ghost 2
	s[16]<-paste(" -es ",as.character(staNG2[i])," 6 ",as.character(paNG2[i])," -ej ",as.character(staNG2[i]+0.000001)," IDPOP 10",sep="")
	s1[16]<-staNG2[i]
	#Divergence Ghost 2 - Basal Europe
	s[17]<-paste(" -ej ",as.character(stdG2BE[i])," 9 10",sep="")
	s1[17]<-stdG2BE[i]
	#Effective pop. size Ghost 2 + Basal Europe
	s[18]<-paste(" -en ",as.character(stdG2BE[i])," 10 ",as.character(snG2[i]),sep="")
	s1[18]<-stdG2BE[i]
	#Divergence Ghost 1 - Papua (OOA 1)
	s[19]<-paste(" -ej ",as.character(stdOA1[i])," 13 8",sep="")
	s1[19]<-stdOA1[i]
	#Effective pop. size Ghost 1 + Papua
	s[20]<-paste(" -en ",as.character(stdOA1[i])," 8 ",as.character(snG1[i]),sep="")
	s1[20]<-stdOA1[i]
	#Divergence Neanderthal - Neanderthal Related
	s[21]<-paste(" -ej ",as.character(stdNNR[i])," 6 5",sep="")
	s1[21]<-stdNNR[i]
	#Effective pop. size ancestral Neanderthal
	s[22]<-paste(" -en ",as.character(stdNNR[i])," 5 ",as.character(snNNR[i]),sep="")
	s1[22]<-stdNNR[i]
	#Divergence Africa - Ghost 1 
	s[23]<-paste(" -ej ",as.character(stdYG1[i])," 8 7",sep="")
	s1[23]<-stdYG1[i]
	#Divergence Africa - Ghost 2 
	s[24]<-paste(" -ej ",as.character(stdYG2[i])," 10 7",sep="")
	s1[24]<-stdYG2[i]
	#Effective pop. size modern 
	s[25]<-paste(" -en ",as.character(stdYG[i])," 7 ",as.character(snYG[i]),sep="")
	s1[25]<-stdYG[i]
	#Divergence Denisova - Denisova Related 1
	s[26]<-paste(" -ej ",as.character(stdDD1[i])," 3 2",sep="")
	s1[26]<-stdDD1[i]
	#Divergence Denisova - Denisova Related 2
	s[27]<-paste(" -ej ",as.character(stdDD2[i])," 4 2",sep="")
	s1[27]<-stdDD2[i]
	#Effective pop. size ancestral Denisova
	s[28]<-paste(" -en ",as.character(stdDDR[i])," 2 ",as.character(snDDR[i]),sep="")
	s1[28]<-stdDDR[i]
	#Divergence Denisova - Neanderthal
	s[29]<-paste(" -ej ",as.character(stdDN[i])," 5 2",sep="")
	s1[29]<-stdDN[i]
	#Effective pop. size Denisova + Neanderthal
	s[30]<-paste(" -en ",as.character(stdDN[i])," 2 ",as.character(snDN[i]),sep="")
	s1[30]<-stdDN[i]
	#Divergence Archaic - (Denisova + Neanderthal)
	s[31]<-paste(" -ej ",as.character(stdADN[i])," 2 1",sep="")
	s1[31]<-stdADN[i]
	#Effective pop. size ancestral Archaic 
	s[32]<-paste(" -en ",as.character(stdADN[i])," 1 ",as.character(snADN[i]),sep="")
	s1[32]<-stdADN[i]
	#Divergence Archaic - modern 
	s[33]<-paste(" -ej ",as.character(stdAM[i])," 7 1",sep="")
	s1[33]<-stdAM[i]
	#Effective pop. size ancestral
	s[34]<-paste(" -en ",as.character(stdAM[i])," 1 ",as.character(snAM[i]),sep="")
	s1[34]<-stdAM[i]

	sid<-sort(s1,index.return=T)
	s_sort<-s[sid$ix]
	ig<-grep("IDPOP",s_sort)
	y<-13#
	for (k in ig){
		y<-y+1
		s_sort[k]<-sub("IDPOP",as.character(y),s_sort[k])
	}

	part0<-paste(s0,collapse="")
	part1<-paste(s_sort,collapse="")
  
	li1<-paste(ms," ",as.character(4*as.numeric(nchr))," ",
	           as.character(nloci)," -t ",as.character(tnAR[i])," -r ",as.character(srec[i])," ",
	           as.character(ll),part0,part1,sep="")
	##### OUTPUT SUMMARY STATISTICS: FDSS
	if (i==1){
		system(paste(li1," | ",cpd,"../compute_pd.py -np 6 -nc ",nchr," -w 30 -b 50 -s > ",out,".tab",sep=""))
	}else{
		system(paste(li1," | ",cpd,"../compute_pd.py -np 6 -nc ",nchr," -w 30 -b 50 -s >> ",out,".tab",sep=""))
	}
}
