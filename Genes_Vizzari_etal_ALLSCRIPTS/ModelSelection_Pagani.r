require(abcrf)
args<-commandArgs(trailingOnly=TRUE)
nsim<-100000
pods<-0
ntrees<-500
fdir1<-"/ooa_m1/"
fdir2<-"/ooa_m2/"
fdir3<-"/observed/pagani/"

ll<-as.character(args[1]) #locus length
nl<-args[2] #Numeber of loci
r<-as.character(args[3]) #Recombination rate
chr<-as.character(args[4]) #Number of chromosomes per population

parcomb<-length(ll)*length(nl)*length(r)*length(chr)
pc<-c()
q<-0
for (locusl in ll){
  for (nloci in nl){
    for (rec in r){
      for (cr in chr){
        a<-c()
        q<-q+1
        print(paste("Parameter combination: ",as.character(q),"/",as.character(parcomb),sep=""))
        pc[q]<-paste("ll",locusl,"_nl",nloci,"_r",rec,"_nc",cr,sep="")
        
	      #####Summary statistics SD model
        name<-paste(fdir1,"ooa_m1","_ll",locusl,"_nl",nloci,"_r",rec,"_nc",cr,".tab",sep="")
        if(file.exists(name)){
          ooa_m1<-read.table(name)
        }else{
          next
        }
	     #####Summary statistics MD model
        name<-paste(fdir2,"ooa_m2","_ll",locusl,"_nl",nloci,"_r",rec,"_nc",cr,".tab",sep="")
        if(file.exists(name)){
          ooa_m2<-read.table(name)
        }else{
          next
        }

        s<-rbind(ooa_m1[(pods+1):nsim,],ooa_m2[(pods+1):nsim,])    
        i<-factor(c(rep(1,nsim-pods),rep(2,nsim-pods)))
        f<-apply(s,2,var)!=0
	
	     #####Reference Table: simulatated data
        da<-data.frame(i,s[,f])
	      
	     #####Observed summary statistics
	  obs<-list()
          for (i in 1:6){
          name<-paste(fdir3,"obs",as.character(i),"_Pagani_2chr_500bp_10k.txt",sep="")
          if(file.exists(name)){
            obs[[i]]<-read.table(name)
          }else{
            next
          }
        t<-rbind(obs[[1]],obs[[2]],obs[[3]],obs[[4]],obs[[5]],obs[[6]])
		
	     #####Training forest for model selection
        a<-abcrf(i~.,data=da,lda=T,ntree=ntrees,paral=T,ncores=10)
		  
             #####Model selection based on observed data
 	c<-predict(object=a,obs=t[,f],training=da,ntree=ntrees,paral=T,paral.predict=T,ncores=10,ncores.predict=10)
		
        write.table(a$model.rf$confusion.matrix,paste("modsel_pagani_ll",locusl,"_nl",nloci,"_r",rec,"_nc",cr,".confMatrix",sep=""))
        print(c)
      }
    }
  }
}
