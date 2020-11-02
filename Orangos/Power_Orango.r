#Usage: Rscript --vanilla ./Power_Orango.r <locuslength> <nrloci> <recrate> <nrchr>
require(abcrf)
args<-commandArgs(trailingOnly=TRUE)
nsim<-50000
pods<-1000
ntrees<-500
fdir<-"Orango/"

ll<-as.character(args[1])
nl<-as.character(args[2])
r<-as.character(args[3])
chr<-as.character(args[4])

parcomb<-length(ll)*length(nl)*length(r)*length(chr)
ris<-matrix(ncol=13,nrow=parcomb)#rows=parcomb, columns=classification error, prior error rate, power and mean posterior probability (pods)
colnames(ris)<-c("mod_1a_ClErr","mod_2a_ClErr","mod_1b_ClErr","mod_2b_ClErr","PriorErrRate","mod_1a_pow","mod_2a_pow","mod_1b_pow","mod_2b_pow","mod_1a_mpp","mod_2a_mpp","mod_1b_mpp","mod_2b_mpp")
ris1<-matrix(ncol=8,nrow=pods)#model preferred and posterior probability for each model (pods)
colnames(ris1)<-c("mod_1a_modind","mod_1a_postpr","mod_2a_modind","mod_2a_postpr","mod_1b_modind","mod_1b_postpr","mod_2b_modind","mod_2b_postpr")
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
                                name<-paste(fdir,"mod_1a","_ll",locusl,"_nl",nloci,"_r",rec,"_nc",cr,".tab",sep="")
                                if(file.exists(name)){
                                        m1a<-read.table(name)
                                }else{
                                        next
                                }
                                name<-paste(fdir,"mod_2a","_ll",locusl,"_nl",nloci,"_r",rec,"_nc",cr,".tab",sep="")
                                if(file.exists(name)){
                                        m2a<-read.table(name)
                                }else{
                                        next
                                }
                                name<-paste(fdir,"mod_1b","_ll",locusl,"_nl",nloci,"_r",rec,"_nc",cr,".tab",sep="")
                                if(file.exists(name)){
                                        m1b<-read.table(name)
                                }else{
                                        next
                                }
                                name<-paste(fdir,"mod_2b","_ll",locusl,"_nl",nloci,"_r",rec,"_nc",cr,".tab",sep="")
                                if(file.exists(name)){
                                        m2b<-read.table(name)
                                }else{
                                        next
                                }
                                s<-rbind(m1a[(pods+1):nsim,],m2a[(pods+1):nsim,],m1b[(pods+1):nsim,],m2b[(pods+1):nsim,])
                                i<-factor(c(rep(1,nsim-pods),rep(2,nsim-pods),rep(3,nsim-pods),rep(4,nsim-pods)))
                                f<-apply(s,2,var)!=0
                                da<-data.frame(i,s[,f])
                                
                                #Training forest
                                a<-abcrf(i~.,data=da,lda=F,ntree=ntrees,paral=T,ncores=16)
                                
                                #Power model 1a
                                b<-predict(object=a,obs=m1a[1:pods,f],training=da,ntree=ntrees,paral=T,paral.predict=T,ncores=8,ncores.predict=8)
                                ris1[,1]<-b$allocation
                                ris1[,2]<-b$post.prob
                                v1<-sum(b$allocation==1)/length(b$allocation)
                                p1<-mean(b$post.prob[b$allocation==1])
                                
                                 #Power model 2a
                                b<-predict(object=a,obs=m2a[1:pods,f],training=da,ntree=ntrees,paral=T,paral.predict=T,ncores=8,ncores.predict=8)
                                ris1[,3]<-b$allocation
                                ris1[,4]<-b$post.prob
                                v2<-sum(b$allocation==2)/length(b$allocation)
                                p2<-mean(b$post.prob[b$allocation==2])
                                
                                #Power model 1b
                                b<-predict(object=a,obs=m1b[1:pods,f],training=da,ntree=ntrees,paral=T,paral.predict=T,ncores=8,ncores.predict=8)
                                ris1[,5]<-b$allocation
                                ris1[,6]<-b$post.prob
                                v3<-sum(b$allocation==3)/length(b$allocation)
                                p3<-mean(b$post.prob[b$allocation==3])
                                
                                #Power model 2b
                                b<-predict(object=a,obs=m2b[1:pods,f],training=da,ntree=ntrees,paral=T,paral.predict=T,ncores=8,ncores.predict=8)
                                ris1[,7]<-b$allocation
                                ris1[,8]<-b$post.prob
                                v4<-sum(b$allocation==4)/length(b$allocation)
                                p4<-mean(b$post.prob[b$allocation==4])
                                
                                #Output files power Orango models
                                k<-a$model.rf$confusion[,"class.error"]
                                ris[q,]<-c(k[1],k[2],k[3],k[4],a$prior.err,v1,v2,v3,v4,p1,p2,p3,p4)
                                write.table(ris1,paste("ll",locusl,"_nl",nloci,"_r",rec,"_nc",cr,"_FDSS.raw",sep=""),row.names=F,quote=F,sep="\t")
                        }
                }
        }
}
rownames(ris)<-pc
#print(ris)
write.table(ris,paste("ll",locusl,"_nl",nloci,"_r",rec,"_nc",cr,"_FDSS.summary",sep=""),quote=F,sep="\t")
