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

samp_vec_int<-function(x=1:10,y=1){
#x is a vector, y is an integer
out<-c()
for (i in 1:length(x)){
	if (x[i]!=y){
		out[i]<-sample(x[i]:y,1,replace=T)
	}else{
		out[i]<-x[i]
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
ms<-"/opt/software/genetics/ms/ms"
cpd<-"./"
mod<-"mod_2b"
nchr<-as.character(args[1])
tgen<-25
mu<-1.5e-8
recomb<-as.numeric(args[4])
ll<-as.numeric(args[2])#locus length
nsims<-as.numeric(args[5])#number of ABC simulations
nloci<-as.numeric(args[3])#loci to simulate in each sim
out<-paste(mod,"_ll",as.character(ll),"_nl",as.character(nloci),"_r",as.character(recomb),"_nc",nchr,sep="")

##PARAMETERS

#Ne Present Time
Ne1BO<-sample(300:32000,nsims,replace=T) 
Ne2BO<-sample(300:32000,nsims,replace=T) 
Ne3BO<-sample(300:32000,nsims,replace=T) 
Ne4BO<-sample(300:32000,nsims,replace=T) 
NeST<-sample(300:32000,nsims,replace=T) 
Ne1NT<-sample(300:32000,nsims,replace=T) 
Ne2NT<-sample(300:32000,nsims,replace=T) 

#Migrations 

#Mig SUBPOP BO; SUBPOP NT
MigBO<-runif(nsims,min=exp(log(10^-4)), max=exp(log(0.1))) #loguniform
MigNT<-runif(nsims,min=exp(log(10^-4)), max=exp(log(0.1)))

#Mig ST-subPop NT 
Mig56<-runif(nsims,min=exp(log(10^-5)), max=exp(log(0.1)))
Mig65<-runif(nsims,min=exp(log(10^-5)), max=exp(log(0.1)))
Mig57<-runif(nsims,min=exp(log(10^-5)), max=exp(log(0.1)))
Mig75<-runif(nsims,min=exp(log(10^-5)), max=exp(log(0.1)))

#Mig ST popAnc NT
MigSTNT<-runif(nsims,min=exp(log(10^-5)), max=exp(log(0.1)))
MigNTST<-runif(nsims,min=exp(log(10^-5)), max=exp(log(0.1)))

#Mig BO-ST 
MigBOST<-runif(nsims,min=exp(log(10^-6)), max=exp(log(10^-2)))
MigSTBO<-runif(nsims,min=exp(log(10^-6)), max=exp(log(10^-2)))

#Bottleneck Intensity Borneo
NeancBO<-samp_vec_int(Ne1BO,320000)
rBO<-Ne1BO/NeancBO


#Ne Ancient
NeancST<-samp_vec_int(NeST,100000)
Neanc1NT<-samp_vec_int(Ne1NT,320000)
Neanc2NT<-samp_vec_int(Ne2NT,320000)
NeancNT<-sample(1000:100000,nsims,replace=T)

#Events Times
tsep4BO<-sample(8750:400000,nsims,replace=T) #T end Bott. BO 

BottDur<-sample(250:100000,nsims,replace=T) #Duration of Bottleneck
tbottend<-tsep4BO+BottDur

tsepBOST<-sample(1500000:4000000,nsims,replace=T)



tStopMig<-samp_vec_vec(tbottend,tsepBOST)

tBotNT<-sample(250:100000,nsims,replace=T) 
tBotST<-tBotNT

tstrNT<-sample(100000:1500000,nsims,replace=T)
tsepNTST<-samp_vec_vec(tstrNT,tsepBOST)


##SCALED PARAMETERS

theta<-4*Ne1BO*mu*ll
srec<-4*Ne1BO*(recomb*(ll-1))

sNe1BO<-Ne1BO*4*mu*ll/theta
sNe2BO<-Ne2BO*4*mu*ll/theta
sNe3BO<-Ne3BO*4*mu*ll/theta
sNe4BO<-Ne4BO*4*mu*ll/theta
sNeST<-NeST*4*mu*ll/theta
sNe1NT<-Ne1NT*4*mu*ll/theta
sNe2NT<-Ne2NT*4*mu*ll/theta


sMigBO<-MigBO*4*Ne1BO
sMigNT<-MigNT*4*Ne1BO

sMig56<-Mig56*4*Ne1BO
sMig65<-Mig65*4*Ne1BO  
sMig57<-Mig57*4*Ne1BO       
sMig75<-Mig75*4*Ne1BO

sMigSTNT<-MigSTNT*4*Ne1BO
sMigNTST<-MigNTST*4*Ne1BO

sMigBOST<-MigBOST*4*Ne1BO
sMigSTBO<-MigSTBO*4*Ne1BO

sNeancBO<-NeancBO*4*mu*ll/theta
sNeancST<-NeancST*4*mu*ll/theta
sNeancNT<-NeancNT*4*mu*ll/theta
sNeanc1NT<-Neanc1NT*4*mu*ll/theta
sNeanc2NT<-Neanc2NT*4*mu*ll/theta

stsep4BO<-(tsep4BO/tgen)/(4*Ne1BO)


stbottend<-(tbottend/tgen)/(4*Ne1BO)

stsepBOST<-(tsepBOST/tgen)/(4*Ne1BO)

stStopMig<-(tStopMig/tgen)/(4*Ne1BO)

stBotNT<-(tBotNT/tgen)/(4*Ne1BO)
stBotST<-stBotNT
ststrNT<-(tstrNT/tgen)/(4*Ne1BO)
stsepNTST<-(tsepNTST/tgen)/(4*Ne1BO)



partable<-cbind(Ne1BO,Ne2BO,Ne3BO,Ne4BO,NeST,Ne1NT,Ne2NT,MigBO,MigNT,Mig56,Mig65,Mig57,Mig75,MigSTNT,MigNTST,MigBOST,MigSTBO,NeancBO,rBO,NeancST,Neanc1NT,Neanc2NT,
tsep4BO,BottDur,tbottend,tsepBOST,tStopMig,tBotNT,tBotST,tstrNT,tsepNTST)
colnames(partable)<-c("Ne1BO","Ne2BO","Ne3BO","Ne4BO","NeST","Ne1NT","Ne2NT","MigBO","MigNT","Mig56","Mig65","Mig57","Mig75","MigSTNT","MigNTST","MigBOST","MigSTBO","NeancBO","rBO","Neanc1NT","Neanc2NT","NeancST",
"tsep4BO","BottDur","tbottend","tsepBOST","tStopMig","tBotNT","tBotST","tstrNT","tsepNTST")
partablescaled<-cbind(sNe1BO,sNe2BO,sNe3BO,sNe4BO,sNeST,sNe1NT,sNe2NT,sMigBO,sMigNT,sMig56,sMig65,sMig57,sMig75,sMigSTNT,sMigNTST,sMigBOST,sMigSTBO,sNeancBO,rBO,sNeanc1NT,sNeanc2NT,sNeancST,stsep4BO,stbottend,stsepBOST,stStopMig,stBotNT,stBotST,ststrNT,stsepNTST)

write.table(partable,paste(out,".1.param",sep=""),row.names=F,quote=F,sep="\t")
write.table(partablescaled,paste(out,".1.paramscaled",sep=""),row.names=F,col.names=T,quote=F,sep="\t")

i<-1
for (i in 1:nsims){	
	s<-c()
	s[1]<-paste(" -ej ",as.character(stsep4BO[i])," 4 1 ",sep="") 
	s[2]<-paste(" -ej ",as.character(stsep4BO[i])," 3 1 ",sep="")
	s[3]<-paste(" -ej ",as.character(stsep4BO[i])," 2 1 ",sep="")
	s[4]<-paste(" -en ",as.character(stsep4BO[i])," 1 ", as.character(rBO[i]),sep="")

	s[5]<-paste(" -em ",as.character(stStopMig[i])," 1 5 ", as.character(sMigBOST[i]),sep="")
	s[6]<-paste(" -em ",as.character(stStopMig[i])," 5 1 ", as.character(sMigSTBO[i]),sep="")
	
	s[7]<-paste(" -en ",as.character(stBotST[i])," 5 ", as.character(sNeancST[i]),sep="")
	s[8]<-paste(" -en ",as.character(stBotNT[i])," 6 ", as.character(sNeanc1NT[i]),sep="")
	s[9]<-paste(" -en ",as.character(stBotNT[i])," 7 ", as.character(sNeanc2NT[i]),sep="")

	s[10]<-paste(" -en ",as.character(stbottend[i])," 1 ", as.character(sNeancBO[i]),sep="")
	s[11]<-paste(" -ej ",as.character(stsepBOST[i])," 5 1 ",sep="")
	s[12]<-paste(" -ej ",as.character(ststrNT[i])," 7 6 ",sep="")
	s[13]<-paste(" -en ",as.character(ststrNT[i])," 6 ", as.character(sNeancNT[i]),sep="")
	s[14]<-paste(" -em ",as.character(ststrNT[i])," 5 6 ", as.character(sMigSTNT[i]),sep="")
	s[15]<-paste(" -em ",as.character(ststrNT[i])," 6 5 ", as.character(sMigNTST[i]),sep="")
	s[16]<-paste(" -ej ",as.character(stsepNTST[i])," 6 5 ",sep="")

	s1<-c()
	
	s1[1]<-stsep4BO[i]
	s1[2]<-stsep4BO[i]
	s1[3]<-stsep4BO[i]
	s1[4]<-stsep4BO[i]

	s1[5]<-stStopMig[i]
	s1[6]<-stStopMig[i]
	
	s1[7]<-stBotST[i]
	s1[8]<-stBotNT[i]
	s1[9]<-stBotNT[i]

	s1[10]<-stbottend[i]
	s1[11]<-stsepBOST[i]
	s1[12]<-ststrNT[i]
	s1[13]<-ststrNT[i]
	s1[14]<-ststrNT[i]
	s1[15]<-ststrNT[i]
	s1[16]<-stsepNTST[i]
	
	sid<-sort(s1,index.return=T)
	s_sort<-s[sid$ix]
	
		

	part1<-paste(s_sort,collapse="")
	
	li1<-paste(ms," ",as.character(7*as.numeric(nchr))," ",as.character(nloci)," -t ",as.character(theta[i])," -r ",as.character(srec[i])," ",as.character(ll)," -I 7 ",nchr," ",nchr," ",nchr," ",nchr," ",nchr," ",nchr," ",nchr," -n 1 ",as.character(sNe1BO[i])," -n 2 ",as.character(sNe2BO[i])," -n 3 ",as.character(sNe3BO[i])," -n 4 ",as.character(sNe4BO[i])," -n 5 ",as.character(sNeST[i])," -n 6 ",as.character(sNe1NT[i])," -n 7 ",as.character(sNe2NT[i])," -m 1 2 ",as.character(sMigBO[i])," -m 2 1 ",as.character(sMigBO[i])," -m 1 3 ",as.character(sMigBO[i])," -m 3 1 ",as.character(sMigBO[i])," -m 1 4 ",as.character(sMigBO[i])," -m 4 1 ",as.character(sMigBO[i])," -m 2 3 ",as.character(sMigBO[i])," -m 3 2 ",as.character(sMigBO[i])," -m 2 4 ",as.character(sMigBO[i])," -m 4 2 ",as.character(sMigBO[i])," -m 3 4 ",as.character(sMigBO[i])," -m 4 3 ",as.character(sMigBO[i])," -m 5 6 ",as.character(sMig56[i])," -m 6 5 ",as.character(sMig65[i])," -m 5 7 ",as.character(sMig57[i])," -m 7 5 ",as.character(sMig75[i])," -m 6 7 ",as.character(sMigNT[i])," -m 7 6 ",as.character(sMigNT[i]), part1, sep="") 

print(i)
#print(li1)

	if (i==1){
		system(paste(li1," | ",cpd,"compute_ss.py -np 7 -nc ",nchr," -w 30 -b 100 -s > ",out,".1.tab",sep=""))
	}else{
		system(paste(li1," | ",cpd,"compute_ss.py -np 7 -nc ",nchr," -w 30 -b 100 -s >> ",out,".1.tab",sep=""))
	}
}

