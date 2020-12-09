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

#Rscript --vanilla ./ooa_m1.r <nchr> <locuslength> <nlociXsim> <rec.rate=1.12e-8> <nsimulations>

args<-commandArgs(trailingOnly=TRUE)
ms<-"/software/msdir/ms"
cpd<-"./"
mod<-"ooa_m1"
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

#Effective population sizes
nAR<-sample(500:50000,nsims,replace=T)
nD<-sample(500:50000,nsims,replace=T)
nD1<-sample(500:50000,nsims,replace=T)
nD2<-sample(500:50000,nsims,replace=T)
nN<-sample(500:50000,nsims,replace=T)
nNR<-sample(500:50000,nsims,replace=T)
nY<-sample(500:50000,nsims,replace=T)
nBE<-sample(500:50000,nsims,replace=T)
nG<-sample(500:50000,nsims,replace=T)
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
m79<-runif(nsims, 10^-6, 10^-3)*4*nAR
m97<-runif(nsims, 10^-6, 10^-3)*4*nAR
m910<-runif(nsims, 10^-6, 10^-3)*4*nAR
m109<-runif(nsims, 10^-6, 10^-3)*4*nAR
m1011<-runif(nsims, 10^-6, 10^-3)*4*nAR
m1110<-runif(nsims, 10^-6, 10^-3)*4*nAR
m1112<-runif(nsims, 10^-6, 10^-3)*4*nAR
m1211<-runif(nsims, 10^-6, 10^-3)*4*nAR

#Migration Eurasians
m1_910<-runif(nsims, 10^-6, 10^-3)*4*nAR
m1_109<-runif(nsims, 10^-6, 10^-3)*4*nAR
m1_1012<-runif(nsims, 10^-6, 10^-3)*4*nAR
m1_1210<-runif(nsims, 10^-6, 10^-3)*4*nAR

#Bottleneck intensity
rP<-1/runif(nsims,min=2,max=100) #Papua
rEA<-1/runif(nsims,min=2,max=100) #Eurasia
rG<-1/runif(nsims,min=2,max=100) #African Ghost

#Divergence time Ancient lineages (Neanderthal, Denisova and Archaic)
tdNNR<-110000
tdDDR<-393000
tdDN<-495000
tdADN<-580000
tdAM<-638000

#Divergence time modern lineages
tdYG<-sample(50000:145000,nsims,replace=T)
tGbot<-(tdYG-tbot)
tdGBE<-tGbot
tmin<-cbind()
i<-1
for (i in 1:nsims){
  tmin[i]<-min(tGbot[i],tdNNR-tgen)}

tdOA<-samp_int_vec(45000,tmin)
tOAbot<-(tdOA-tbot)
tdEA<-sample(30000:40000,nsims,replace=T)


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
taD1P<-samp_int_vec(30000,tOAbot)
paD1P<-runif(nsims, 10^-3, 10^-1)
paD1P<-1-paD1P

#Archaic Papua
taARP<-samp_vec_vec(taD1P,tOAbot)
paARP<-runif(nsims, 10^-3, 10^-1)
paARP<-1-paARP

#Neanderthal-Asia
taNEA<-samp_vec_vec(tdEA,tOAbot)
paNEA<-runif(nsims, 10^-3, 10^-1)
paNEA<-1-paNEA

#Neanderthal-African Ghost
taNG<-samp_vec_vec(tdOA,tmin)
paNG<-runif(nsims, 10^-3, 10^-1)
paNG<-1-paNG

#Parameters Scaled for ms

tnAR<-4*nAR*mu*ll #theta

snD<-nD*4*mu*ll/tnAR
snD1<-nD1*4*mu*ll/tnAR
snD2<-nD2*4*mu*ll/tnAR
snN<-nN*4*mu*ll/tnAR
snNR<-nNR*4*mu*ll/tnAR
snY<-nY*4*mu*ll/tnAR
snG<-nG*4*mu*ll/tnAR
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
stGbot<-(tGbot/tgen)/(4*nAR)
stdGBE<-(tdGBE/tgen)/(4*nAR)
stdOA<-(tdOA/tgen)/(4*nAR)
stOAbot<-(tOAbot/tgen)/(4*nAR)
stdEA<-(tdEA/tgen)/(4*nAR)

staD2A<-(taD2A/tgen)/(4*nAR)
staBEE<-(taBEE/tgen)/(4*nAR)
staD1P<-(taD1P/tgen)/(4*nAR) 
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

#### OUTPUT PARAMETER's FILES
partable<-cbind(nAR,nD,nD1,nD2,nN,nNR,nY,nG,nBE,nE,nA,nP,
                nYG,nNNR,nDDR,nDN,nADN,nAM,
                rYG,rGBE,rNNR,rDDR,rNDN,rADN,rAM,rP,rEA,rG,
                tdYG,tGbot,tdGBE,tdOA,tOAbot,tdEA,
                taD2A,paD2A,taBEE,paBEE,taD1P,paD1P,taARP,paARP,taNEA,paNEA,taNG,paNG,
                stN,stD,tdNNR,tdDDR,tdDN,tdADN,tdAM)
partablescaled<-cbind(tnAR,snD,snD1,snD2,snN,snNR,snY,snG,snBE,snE,snA,snP,
                      snYG,snNNR,snDDR,snDN,snADN,snAM,srec,
                      stdYG,stGbot,stdGBE,stdOA,stOAbot,stdEA,
                      staD2A,staBEE,staD1P,staARP,staNEA,staNG,
                      sstN,sstD,stdNNR,stdDDR,stdDN,stdADN,stdAM,
                      m79,m97,m910,m109,m1011,m1110,m1112,m1211,tm1,m1_910,m1_109,m1_1012,m1_1210)
					  
write.table(partable,paste(out,".param",sep=""),row.names=F,col.names=T,quote=F,sep="\t")
write.table(partablescaled,paste(out,".paramscaled",sep=""),row.names=F,col.names=T,quote=F,sep="\t")

#### Creation ms command line:
i<-1
for (i in 1:nsims){
	s0<-c()
	#Modern populations
	s0[1]<-paste(" -I 12 0 0 0 0 0 0 ",nchr," 0 0 ",nchr," ",nchr," ",nchr,sep="")
	#Modern effective population size
	s0[2]<-paste(" -n 7 ",as.character(snY[i]),sep="")
	s0[3]<-paste(" -n 8 ",as.character(snG[i]),sep="")
	s0[4]<-paste(" -n 9 ",as.character(snBE[i]),sep="")
	s0[5]<-paste(" -n 10 ",as.character(snE[i]),sep="")
	s0[6]<-paste(" -n 11 ",as.character(snA[i]),sep="")
	s0[7]<-paste(" -n 12 ",as.character(snP[i]),sep="")
	#Migration between modern populations
	s0[8]<-paste(" -m 7 9 ",as.character(m79[i]),sep="")
	s0[9]<-paste(" -m 9 7 ",as.character(m97[i]),sep="")
	s0[10]<-paste(" -m 9 10 ",as.character(m910[i]),sep="")
	s0[11]<-paste(" -m 10 9 ",as.character(m109[i]),sep="")
	s0[12]<-paste(" -m 10 11 ",as.character(m1011[i]),sep="")
	s0[13]<-paste(" -m 11 10 ",as.character(m1110[i]),sep="")
	s0[14]<-paste(" -m 11 12 ",as.character(m1112[i]),sep="")
	s0[15]<-paste(" -m 12 11 ",as.character(m1211[i]),sep="") 
	#Sampling time Denisova
	s0[16]<-paste(" -eA ",as.character(sstD[i])," 2 ",nchr,sep="")
	#Sampling time Neandertal
	s0[17]<-paste(" -eA ",as.character(sstN[i])," 5 ",nchr,sep="")
	
	s<-c()
	s1<-c()
	#Admixture BasalEurope - Europe
	s[1]<-paste(" -es ",as.character(staBEE[i])," 8 ",as.character(paBEE[i])," -ej ",as.character(staBEE[i]+0.000001)," IDPOP 10",sep="")
	s1[1]<-staBEE[i]
	#Admixture DenisovaRelated 2 - Asia
	s[2]<-paste(" -es ",as.character(staD2A[i])," 4 ",as.character(paD2A[i])," -ej ",as.character(staD2A[i]+0.000001)," IDPOP 11",sep="")
	s1[2]<-staD2A[i]
	#Divergence time Europe-Asia
	s[3]<-paste(" -ej ",as.character(stdEA[i])," 11 10",sep="")
	s1[3]<-stdEA[i]
	#Effective population size Eurasia
	s[4]<-paste(" -en ",as.character(stdEA[i])," 10 ",as.character(snE[i]),sep="")
	s1[4]<-stdEA[i]
	#Migrations between Eurasia-Africa and Eurasia-Papua
	s[5]<-paste(" -em ",as.character(tm1[i])," 9 10 ",as.character(m1_910[i]),sep="")
	s1[5]<-tm1[i]
	s[6]<-paste(" -em ",as.character(tm1[i])," 10 9 ",as.character(m1_109[i]),sep="")
	s1[6]<-tm1[i]
	s[7]<-paste(" -em ",as.character(tm1[i])," 10 12 ",as.character(m1_1012[i]),sep="")
	s1[7]<-tm1[i]
	s[8]<-paste(" -em ",as.character(tm1[i])," 12 10 ",as.character(m1_1210[i]),sep="")
	s1[8]<-tm1[i]
	#Admixtures Denisova Related 1 - Papua
	s[9]<-paste(" -es ",as.character(staD1P[i])," 3 ",as.character(paD1P[i])," -ej ",as.character(staD1P[i]+0.000001)," IDPOP 12",sep="")
	s1[9]<-staD1P[i]
	#Admixtures Arcaici - Papua
	s[10]<-paste(" -es ",as.character(staARP[i])," 1 ",as.character(paARP[i])," -ej ",as.character(staARP[i]+0.000001)," IDPOP 12",sep="")
	s1[10]<-staARP[i]
	#Admixtures Neanderthal Related - Eurasia
	s[11]<-paste(" -es ",as.character(staNEA[i])," 6 ",as.character(paNEA[i])," -ej ",as.character(staNEA[i]+0.000001)," IDPOP 10",sep="")
	s1[11]<-staNEA[i]
	#Effective pop. size Eurasia during bottleneck
	s[12]<-paste(" -en ",as.character(stOAbot[i])," 10 ",as.character(rEA[i]),sep="")
	s1[12]<-stOAbot[i]
	#Effective pop. size Papua during bottleneck
	s[13]<-paste(" -en ",as.character(stOAbot[i])," 12 ",as.character(rP[i]),sep="")
	s1[13]<-stOAbot[i]
	#Divergence Ghost - Eurasiatici (OOA)
	s[14]<-paste(" -ej ",as.character(stdOA[i])," 10 9",sep="")
	s1[14]<-stdOA[i]
	#Divergence Ghost - Papuan (OOA)
	s[15]<-paste(" -ej ",as.character(stdOA[i])," 12 9",sep="")
	s1[15]<-stdOA[i]
	#Admixtures Neanderthal Related - Ghost
	s[16]<-paste(" -es ",as.character(staNG[i])," 6 ",as.character(paNG[i])," -ej ",as.character(staNG[i]+0.000001)," IDPOP 9",sep="")
	s1[16]<-staNG[i]
	#Divergence Ghost - Basal Europe
	s[17]<-paste(" -ej ",as.character(stdGBE[i])," 8 9",sep="")
	s1[17]<-stdGBE[i]
	#Effective pop. size Ghost during bottleneck
	s[18]<-paste(" -en ",as.character(stGbot[i])," 9 ",as.character(rG[i]),sep="")
	s1[18]<-stGbot[i]
	#Divergence Neanderthal - Neanderthal related
	s[19]<-paste(" -ej ",as.character(stdNNR[i])," 6 5",sep="")
	s1[19]<-stdNNR[i]
	#Effective pop. size Neanderthal
	s[20]<-paste(" -en ",as.character(stdNNR[i])," 5 ",as.character(snNNR[i]),sep="")
	s1[20]<-stdNNR[i]
	#Divergence Africa - Ghost
	s[21]<-paste(" -ej ",as.character(stdYG[i])," 9 7",sep="")
	s1[21]<-stdYG[i]
	#Effective pop. size Ancestral modern population
	s[22]<-paste(" -en ",as.character(stdYG[i])," 7 ",as.character(snYG[i]),sep="")
	s1[22]<-stdYG[i]
	#Divergence Denisova - Denisova Related 2
	s[23]<-paste(" -ej ",as.character(stdDDR[i])," 4 2",sep="")
	s1[23]<-stdDDR[i]
	#Divergence Denisova - Denisova Related 1
	s[24]<-paste(" -ej ",as.character(stdDDR[i])," 3 2",sep="")
	s1[24]<-stdDDR[i]
	#Effective pop. size ancestral Denisova
	s[25]<-paste(" -en ",as.character(stdDDR[i])," 2 ",as.character(snDDR[i]),sep="")
	s1[25]<-stdDDR[i]
	#Divergence Denisova - Neanderthal
	s[26]<-paste(" -ej ",as.character(stdDN[i])," 5 2",sep="")
	s1[26]<-stdDN[i]
	#Effective pop. size Ancestral Denisova and Neanderthal
	s[27]<-paste(" -en ",as.character(stdDN[i])," 2 ",as.character(snDN[i]),sep="")
	s1[27]<-stdDN[i]
	#Divergence Unknown Archaic - Ancestral Denisova and Neanderthal pop.
	s[28]<-paste(" -ej ",as.character(stdADN[i])," 2 1",sep="")
	s1[28]<-stdADN[i]
	#Effective pop. size Ancestral ancient population
	s[29]<-paste(" -en ",as.character(stdADN[i])," 1 ",as.character(snADN[i]),sep="")
	s1[29]<-stdADN[i]
	#Divergence Ancient - Modern pop.
	s[30]<-paste(" -ej ",as.character(stdAM[i])," 7 1",sep="")
	s1[30]<-stdAM[i]
	#Effective pop. size Ancestral populations (Ancient+Modern)
	s[31]<-paste(" -en ",as.character(stdAM[i])," 1 ",as.character(snAM[i]),sep="")
	s1[31]<-stdAM[i]
	
	sid<-sort(s1,index.return=T)
	s_sort<-s[sid$ix]
	ig<-grep("IDPOP",s_sort)
	y<-12#
	for (k in ig){
		y<-y+1
		s_sort[k]<-sub("IDPOP",as.character(y),s_sort[k])
	}
	
	part0<-paste(s0,collapse="")
	part1<-paste(s_sort,collapse="")

	li1<-paste(ms," ",as.character(4*as.numeric(nchr))," ",
	           as.character(nloci)," -t ",as.character(tnAR[i])," -r ",as.character(srec[i])," ",
	           as.character(ll),part0,part1,sep="")
	
	if(i==1){
##### OUTPUT SUMMARY STATISTICS: FDSS
	 system(paste(li1," | ",cpd,"../compute_pd.py -np 6 -nc ",nchr," -w 30 -b 100 -s > ",out,".tab",sep=""))
 }else{
		 system(paste(li1," | ",cpd,"../compute_pd.py -np 6 -nc ",nchr," -w 30 -b 100 -s >> ",out,".tab",sep=""))
	 }
}
