#Usage: Rscript --vanilla ./Power_humans.r <locuslength> <nrloci> <recrate> <nrchr>
require(abcrf)
args<-commandArgs(trailingOnly=TRUE)
nsim<-50000
pods<-1000
ntrees<-500
fdir1<-"ms_ooa1/"
fdir2<-"ms_ooa2/"

ll<-as.character(args[1])
nl<-args[2]
r<-as.character(args[3])
chr<-as.character(args[4])

parcomb<-length(ll)*length(nl)*length(r)*length(chr)
ris<-matrix(ncol=7,nrow=parcomb)#rows=parcomb, columns=classification error, prior error rate, power and mean posterior probability
colnames(ris)<-c("sdm_ClErr","mdm_ClErr","PriorErrRate","sdm_pow","mdm_pow","sdm_mpp","mdm_mpp")
ris1<-matrix(ncol=4,nrow=pods)#model preferred and posterior probability for each model
colnames(ris1)<-c("sdm_modind","sdm_postpr","mdm_modind","mdm_postpr")
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
                                print(name)
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
                                
                                #Trining forest
                                a<-abcrf(i~.,data=da,lda=T,ntree=ntrees,paral=T,ncores=10)
                                
                                #Power SDM model
                                b<-predict(object=a,obs=ooa_m1[1:pods,f],training=da,ntree=ntrees,paral=T,paral.predict=T,ncores=5,ncores.predict=5)
                                ris1[,1]<-b$allocation
                                ris1[,2]<-b$post.prob
                                v1<-sum(b$allocation==1)/length(b$allocation)
                                p1<-mean(b$post.prob[b$allocation==1])
                                
                                #Power MDM model
                                b<-predict(object=a,obs=ooa_m2[1:pods,f],training=da,ntree=ntrees,paral=T,paral.predict=T,ncores=5,ncores.predict=5)
                                ris1[,3]<-b$allocation
                                ris1[,4]<-b$post.prob
                                v2<-sum(b$allocation==2)/length(b$allocation)
                                p2<-mean(b$post.prob[b$allocation==2])
                                
                                #Output file power Out of Africa model
                                k<-a$model.rf$confusion[,"class.error"]
                                ris[q,]<-c(k[1],k[2],a$prior.err,v1,v2,p1,p2)
                                write.table(ris1,paste("ll",locusl,"_nl",nloci,"_r",rec,"_nc",cr,"_FDSS.raw",sep=""),row.names=F,quote=F,sep="\t")
                        }
                }
        }
}
rownames(ris)<-pc
write.table(ris,paste("ll",locusl,"_nl",nloci,"_r",rec,"_nc",cr,"_FDSS.summary",sep=""),quote=F,sep="\t")
