require(abcrf)
args<-commandArgs(trailingOnly=TRUE)
nsim<-100000
pods<-0
ntrees<-500
fdir1<-"/ooa_m2/"
fdir2<-"/observed/"

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
        q<-1
        pc[q]<-paste("ll",locusl,"_nl",nloci,"_r",rec,"_nc",cr,sep="")
        name<-paste(fdir1,"ooa_m2","_ll",locusl,"_nl",nloci,"_r",rec,"_nc",cr,".param",sep="")
        print(name)
        if(file.exists(name)){
          par<-read.table(name,h=T)
        }else{
          next
        }
        name<-paste(fdir1,"ooa_m2","_ll",locusl,"_nl",nloci,"_r",rec,"_nc",cr,".paramscaled",sep="")
        print(name)
        if(file.exists(name)){
          parsc<-read.table(name,h=T)
        }else{
          next
        }
        name<-paste(fdir1,"ooa_m2","_ll",locusl,"_nl",nloci,"_r",rec,"_nc",cr,".tab",sep="")
        print(name)
        if(file.exists(name)){
          sust<-read.table(name)
        }else{
          next
        }
        name<-paste(fdir1,"prior_ooa_m2.txt",sep="")
        print(name)
        if(file.exists(name)){
          prior<-read.table(name)
        }else{
          next
        }
        o<-list()
        name<-paste(fdir2,"malaspinas/obs5_Malaspinas_2chr_500bp_10k.txt",sep="")
        print(name)
        if(file.exists(name)){
          o[[1]]<-read.table(name)
        }else{
          next
        }
        name<-paste(fdir2,"pagani/obs6_500_n2_10k.txt",sep="")
        print(name)
        if(file.exists(name)){
          o[[2]]<-read.table(name)
        }else{
          next
        }
        obs<-rbind(o[[1]],o[[2]])
        par.t<-cbind(par[1],par[,7:19],par[,26:47],parsc[,42:51],parsc[,53:56])
        for (i in 37:50) par.t[i]=par.t[i]/(par.t[1]*4)
        for (i in 26:36){
          if (i==26 | i==28 | i==30 | i==32 | i==34 | i==36) par.t[i]=1-par.t[i]
        }
        col<-c(colnames(par.t))
		
        ris2<-matrix(ncol=9,nrow=ncol(par.t))
        colnames(ris2)<-c("Parameter","expectation","median","variance","variance.cdf",
                          "quantile90m","quantile90M","quantile50m","quantile50M")
        ris3<-matrix(ncol=9,nrow=ncol(par.t))
        colnames(ris3)<-c("Parameter","expectation","median","variance","variance.cdf",
                          "quantile90m","quantile90M","quantile50m","quantile50M")
       
	    par.s<-par.t[((pods+1):(nsim)),]
        f<-apply(sust,2,var)!=0
        Sust.s<-sust[((pods+1):(nsim)),f]
        data.obs<-sust[(1:pods),f]
        data.s<-list()
       
	    p<-c()
        for (i in 1:length(par.s)){
          data.s[[i]]<-data.frame(par.s[i], Sust.s)
        
		  p[i]<-c(paste(col[i],"~.",sep=""))
         
   		  a<-regAbcrf(as.formula(p[i]),data=data.s[[i]],ntree=ntrees,paral=T,ncores=10)
          
		  c1<-predict(object=a,obs=obs,training=data.s[[i]],ntree=ntrees,quantiles=c(0.05, 0.95),
                      paral=T,paral.predict=T,ncores=5,ncores.predict=5)
          c2<-predict(object=a,obs=obs,training=data.s[[i]],ntree=ntrees,quantiles=c(0.25, 0.75),
                      paral=T,paral.predict=T,ncores=5,ncores.predict=5)
					  
          ris2[i,]<-c(col[i],c1$expectation[1],c1$med[1],c1$variance[1],c1$variance.cdf[1],
                      c1$quantiles[1,1],c1$quantiles[1,2],c2$quantiles[1,1],c2$quantiles[1,2])
          ris3[i,]<-c(col[i],c1$expectation[2],c1$med[2],c1$variance[2],c1$variance.cdf[2],
                      c1$quantiles[2,1],c1$quantiles[2,2],c2$quantiles[2,1],c2$quantiles[2,2])
          
		  save(a,file=paste("./parametri/",col[i],".allen",sep=""))
          
		  pdf(file=(paste("./parametri/",col[i],"_mal_ooa_m2_ll",locusl,"_nl",nloci,"_r",rec,"_nc",cr,".pdf",sep="")))
          densityPlot(object=a,obs=rbind(obs[1,],obs[1,]),training=data.s[[i]],
                      xlab="",ylab="Density",main=paste("Posterior density of ",col[i]," Malaspinas",sep=""),
                      paral=T,ncores=10)
          abline(v=c1$expectation[1],col="red")
          abline(v=c1$med[1],col="blue")
          abline(h=1/(prior[i,2]-prior[i,1]),col="darkgreen",lty=5)
          dev.off()
         
		 pdf(file=(paste("./parametri/",col[i],"_pag_ooa_m2_ll",locusl,"_nl",nloci,"_r",rec,"_nc",cr,".pdf",sep="")))
          densityPlot(object=a,obs=rbind(obs[2,],obs[2,]),training=data.s[[i]],
                      xlab="",ylab="Density",main=paste("Posterior density of ",col[i]," Pagani",sep=""),
                      paral=T,ncores=10)
          abline(v=c2$expectation[2],col="red")
          abline(v=c2$med[2],col="blue")
          abline(h=1/(prior[i,2]-prior[i,1]),col="darkgreen",lty=5)
          dev.off()
        
		write.table(ris2,paste("ooa_m2_mal_ll",locusl,"_nl",nloci,"_r",rec,"_nc",cr,".summary",sep=""),
                      row.names=F,quote=F,sep="\t")
          write.table(ris3,paste("ooa_m2_pag_ll",locusl,"_nl",nloci,"_r",rec,"_nc",cr,".summary",sep=""),
                      row.names=F,quote=F,sep="\t")
        }
      }
    }
  }
}
