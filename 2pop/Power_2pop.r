#Usage: Rscript --vanilla ./POwer_2pop.r <locuslength> <nrloci> <rec.rate> <nrchr>
require(abcrf)
args<-commandArgs(trailingOnly=TRUE)
nsim<-50000
pods<-1000
ntrees<-500
fdir<-"ms2pop/"

ll<-as.character(args[1])
nl<-args[2]
r<-as.character(args[3])
chr<-as.character(args[4])

parcomb<-length(ll)*length(nl)*length(r)*length(chr)
ris<-matrix(ncol=10,nrow=parcomb)#rows=parameter combination, columns=classification error, prior error rate, true positive rate(pods) and mean posterior probability(pods)
colnames(ris)<-c("m1_ClErr","m2_ClErr","m3_ClErr","PriorErrRate","m1_pow","m2_pow","m3_pow","m1_mpp","m2_mpp","m3_mpp")
ris1<-matrix(ncol=6,nrow=pods)#model preferred and posterior probability for each model (pods)
colnames(ris1)<-c("m1_modind","m1_postpr","m2_modind","m2_postpr","m3_modind","m3_postpr")
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
                                name<-paste(fdir,"m1","_ll",locusl,"_nl",nloci,"_r",rec,"_nc",cr,".tab",sep="")
                                print(name)
                                if(file.exists(name)){
                                        m1<-read.table(name)
                                }else{
                                        next
                                }
                                name<-paste(fdir,"m2","_ll",locusl,"_nl",nloci,"_r",rec,"_nc",cr,".tab",sep="")
                                if(file.exists(name)){
                                        m2<-read.table(name)
                                }else{
                                        next
                                }

                                name<-paste(fdir,"m3","_ll",locusl,"_nl",nloci,"_r",rec,"_nc",cr,".tab",sep="")
                                if(file.exists(name)){
                                        m3<-read.table(name)
                                }else{
                                        next
                                }
                               s<-rbind(m1[(pods+1):nsim,],m2[(pods+1):nsim,],m3[(pods+1):nsim,])
                               i<-factor(c(rep(1,nsim-pods),rep(2,nsim-pods),rep(3,nsim-pods)))
                               f<-apply(s,2,var)!=0
                               da<-data.frame(i,s[,f])
                               
                               #Training Forest
                               a<-abcrf(i~.,data=da,lda=T,ntree=ntrees,paral=T,ncores=6)

                               #Power Diverence model
                               b<-predict(object=a,obs=m1[1:pods,f],training=da,ntree=ntrees,paral=T,paral.predict=T,ncores=3,ncores.predict=3)
                               ris1[,1]<-b$allocation
                               ris1[,2]<-b$post.prob
                               v1<-sum(b$allocation==1)/length(b$allocation)
                               p1<-mean(b$post.prob[b$allocation==1])
                               
                               #Power Diverence model with migration
                               b<-predict(object=a,obs=m2[1:pods,f],training=da,ntree=ntrees,paral=T,paral.predict=T,ncores=3,ncores.predict=3)
                               ris1[,3]<-b$allocation
                               ris1[,4]<-b$post.prob
                               v2<-sum(b$allocation==2)/length(b$allocation)
                               p2<-mean(b$post.prob[b$allocation==2])
                               
                               #Power Diverence model with pulse of admixture
                               b<-predict(object=a,obs=m3[1:pods,f],training=da,ntree=ntrees,paral=T,paral.predict=T,ncores=3,ncores.predict=3)
                               ris1[,7]<-b$allocation
                               ris1[,8]<-b$post.prob
                               v3<-sum(b$allocation==3)/length(b$allocation)
                               p3<-mean(b$post.prob[b$allocation==3])
                               
                               #Output files power analysis two population models
                               k<-a$model.rf$confusion[,"class.error"]
                               ris[q,]<-c(k[1],k[2],k[3],a$prior.err,v1,v2,v3,p1,p2,p3)
                               write.table(ris1,paste("ll",locusl,"_nl",nloci,"_r",rec,"_nc",cr,"_FDSS.raw",sep=""),row.names=F,quote=F,sep="\t")
                       }
               }
       }
}
rownames(ris)<-pc
print(ris)
write.table(ris,paste("ll",locusl,"_nl",nloci,"_r",rec,"_nc",cr,"_FDSS.summary",sep=""),quote=F,sep="\t")

