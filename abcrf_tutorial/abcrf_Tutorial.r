##Toy example: model1= Two divergent populations; model2: two divergent populations exchanging migrants.

require(abcrf)

### MODEL SELECTION ###

m1_sim_SuSt<-read.table("model1_Simulated.SummaryStatistics.tab")
m2_sim_SuSt<-read.table("model2_Simulated.SummaryStatistics.tab")
Obs.Data<-read.table("Observed.SummaryStatistics.tab")


##Reference Table
sim.data.all<-rbind(m1_sim_SuSt,m2_sim_SuSt)
id_models<-factor(c(rep(1,nrow(m1_sim_SuSt)),rep(2,nrow(m2_sim_SuSt))))
f<-apply(sim.data.all,2,var)!=0

Reference.Table<-data.frame(id_models,sim.data.all[,f])

##Training forest for model selection
Forest.ModelSelection<-abcrf(id_models~.,data=Reference.Table,lda=T,ntree=500)
print(Forest.ModelSelection)

##Model selection based on observed data
Predict.ModelSelection<-predict(object=Forest.ModelSelection,obs=rbind(Obs.Data[,f],Obs.Data[,f]),training=Reference.Table,ntree=500)
print(Predict.ModelSelection)

##Plot model selection results
pdf("Plot_ModelSelection_Results.pdf")
plot(Forest.ModelSelection,obs=rbind(Obs.Data[,f],Obs.Data[,f]),training=Reference.Table)
dev.off()


### PARAMETERS ESTIMATION ###
m1_Params<-read.table("model1_Demographic.Parameters.param",header=T)
par.names<-colnames(m1_Params)

data.RefTables<-list()
formula<-c()

## Reference tables for parameters estimation: one for each demographic parameters
for (i in 1:ncol(m1_Params)) {
  data.RefTables[[i]]<-data.frame(as.numeric(m1_Params[,i]),m1_sim_SuSt[,f])
  colnames(data.RefTables[[i]])[1]<-par.names[i]
  formula[i]<-c(paste(par.names[i],"~.",sep=""))
}

###Estimate one of the demographic parameters: n1 (effective population size of the first population)

##Training forest for parameter estimation
Forest.ParamEstimation<-regAbcrf(as.formula(formula[1]), data=data.RefTables[[1]], ntree=500)
print(Forest.ParamEstimation$model.rf)

##Parameter estimation based on observed data
Predict.ParamEstimation<-predict(object=Forest.ParamEstimation,obs=rbind(Obs.Data[,f],Obs.Data[,f]),training=data.RefTables[[1]],ntree=500)
print(Predict.ParamEstimation)

##Plot posterior density
pdf("Plot_PosteriorDensity_Results.pdf")
densityPlot(object=Forest.ParamEstimation,obs=rbind(Obs.Data[,f],Obs.Data[,f]),training=data.RefTables[[1]],xlab=par.names[1],ylab="Density")
dev.off()

