#Usage: Rscript --vanilla ./analysis_ecc.r <locuslength> <nrloci> <recrate> <nrchr>
require(abcrf)
args<-commandArgs(trailingOnly=TRUE)
nsim<-100000
pods<-0
ntrees<-500
fdir1<-"/ooa_m1/"
fdir2<-"/ooa_m2/"
/observed/pagani/<-"/observed/pagani/"

ll<-as.character(args[1])
nl<-args[2]
r<-as.character(args[3])
chr<-as.character(args[4])

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
        
         name<-paste(fdir1,"ooa_m1","_ll",locusl,"_nl",nloci,"_r",rec,"_nc",cr,".tab",sep="")
        if(file.exists(name)){
          ooa_m1<-read.table(name)
        }else{
          next
        }
        name<-paste(fdir2,"ooa_m2","_ll",locusl,"_nl",nloci,"_r",rec,"_nc",cr,".tab",sep="")
        if(file.exists(name)){
          ooa_m2<-read.table(name)
        }else{
          next
        }

        s<-rbind(ooa_m1[(pods+1):nsim,],ooa_m2[(pods+1):nsim,])    
        i<-factor(c(rep(1,nsim-pods),rep(2,nsim-pods)))
        f<-apply(s,2,var)!=0
		
        da<-data.frame(i,s[,f])
		
        t<-rbind(obs[,f],obs1[,f],obs2[,f],obs3[,f],obs4[,f],obs5[,f])
		
        a<-abcrf(i~.,data=da,lda=T,ntree=ntrees,paral=T,ncores=10)
        c<-predict(object=a,obs=t,training=da,ntree=ntrees,paral=T,paral.predict=T,ncores=10,ncores.predict=10)
		
        write.table(a$model.rf$confusion.matrix,paste("modsel_pagani_ll",locusl,"_nl",nloci,"_r",rec,"_nc",cr,".confMatrix",sep=""))
        print(c)
      }
    }
  }
}
