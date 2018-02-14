source("FunctionsBTPS.R")
source("GraphSetUp.R")
source("SimulationSetUp.R")


####
#Simulate data
mySimData<-yieldDataSim(1500,YieldMeanFunc="Linear",YieldVarFunc="NonConst",YieldError="Normal",p=4,meas_error=.000001)

##Fit data
myEst<-fitBTPS(mySimData$Obs_Prem,mySimData$Cov_Rate,mySimData$Yield_Hist,knots =c(2,9),penalty = c(3,3),degree = c(5,3),newXLim = c(.55,.95) ,newZLim =c(100,300))

myVar<-VarianceEstimator(myEst) 

###monte carlo fit of penalized BTPB
##Values used for the set up presented in the paper:reps=1000,seed=9389,samplesize=500 ,pVal=4,MeasurementError=.01,BTPBknots=c(2,10),BTPBpenalty=c(3,3),BTPBdegree=c(5,3),XLim =c(.55,.95) ,ZLim =c(100,300)
##Values used for the creation of the toy data: reps=1000,seed=9389,samplesize=1500,Yield_Mean="Linear",VarFunc="NonConst",Yield_Error="Normal"  ,pVal=4,MeasurementError=.001,BTPBknots=c(2,9),BTPBpenalty=c(3,3),BTPBdegree=c(5,3),XLim =c(.55,.95) ,ZLim =c(100,300)
mcBTPBtest<-simStudyFunc(reps=1000,seed=9389,samplesize=1500,Yield_Mean="Linear",VarFunc="NonConst",Yield_Error="Normal"  ,pVal=4,MeasurementError=.001,BTPBknots=c(2,9),BTPBpenalty=c(3,3),BTPBdegree=c(5,3),tol2nd = .00005,XLim =c(.55,.95) ,ZLim =c(100,300))

####monte carlo fit of kernel Density estimation (takes a while because of the jackknife variance)
##Values used for the set up presented in the paper: Nreps=1000,seed=9389,samplesize=500,zValues=seq(100,300, by=10),wRange=c(0,450),numberOfW=100
##Values used for the creation of the toy data: Nreps=1000,seed=9389,samplesize=500,zValues=seq(100,300, by=10),wRange=c(0,450),numberOfW=100,Yield_Mean="Linear",VarFunc="NonConst",Yield_Error="Normal" 
mcKernelTest<-simStudyKernel(Nreps=1000,seed=9389,samplesize=500,zValues=seq(100,300, by=10),wRange=c(0,450),numberOfW=100,Yield_Mean="Linear",VarFunc="NonConst",Yield_Error="Normal"  )


##load in toy data using the above set ups(1000 replicates to show how the graphs work)
myMCBTPBest<-readRDS("BTPBmc")
myMCkernel<-readRDS("KernelMC")

##data set up to show how graphs 1 and 2 are produced
myGraphs(myMCBTPBest$Estimates,myMCBTPBest$Variance)

##data set up to show how graphs 1 and 2 of the supplimental file are produced
myGraphsKern(myMCkernel$Estimates,myMCkernel$Variance,Yield_Mean="Linear",VarFunc="NonConst",Yield_Error="Normal" )

##show graphs of empirical study
##create simulated data to take the place of the real data
set.seed(9389)
mySimData<-yieldDataSim(1500,YieldMeanFunc="Linear",YieldVarFunc="NonConst",YieldError="Normal",p=4,meas_error=.0001)

#fit it using the penalized BTPS
myEst<-fitBTPS(mySimData$Obs_Prem,mySimData$Cov_Rate,mySimData$Yield_Hist,knots =c(2,9),penalty = c(3,3),degree = c(5,3),newXLim = c(.55,.95) ,newZLim =c(100,300))
myVar<-VarianceEstimator(myEst) 


#Fit the simulated data using the kernel density estimator
fit<-fitKernel(mySimData$Yield_Curr,mySimData$Yield_Hist,knots=10,wRange=c(0,450),fixedZ = seq(100,300, by=10),NumW = 1000)

##Code for Graph 3
par( mfrow = c( 2, 3 ) )
par(mar=c(4.5,2,1,1))
for(i in 1:5){
value<-c(120,160,200,240,280)[i]
BTPBvalue<-c(120,161,202,239,280)[i]
#plot the penalized BTPB fit
plot(value*myEst$newX[round(myEst$newZ,0)==BTPBvalue],myEst$secondDerv[round(myEst$newZ,0)==BTPBvalue]/value^2/4+1.96*sqrt(myVar$variance[round(myEst$newZ,0)==BTPBvalue])/value^2/4,main=paste0("z=",value),type='l',col="green",ylab="",lty=2,lwd=2.5,xlab="",cex.main=1.75,cex.axis=1.5,cex.lab=1.5)
lines(value*myEst$newX[round(myEst$newZ,0)==BTPBvalue],myEst$secondDerv[round(myEst$newZ,0)==BTPBvalue]/value^2/4,col="black",lty=1,lwd=2.5)
lines(value*myEst$newX[round(myEst$newZ,0)==BTPBvalue],myEst$secondDerv[round(myEst$newZ,0)==BTPBvalue]/value^2/4-1.96*sqrt(myVar$variance[round(myEst$newZ,0)==BTPBvalue])/value^2/4,col="blue",lty=2,lwd=2.5)

#plot the kernel density fit
columnValue<-c(2:ncol(fit$estimate))[as.numeric(dimnames(fit$estimate)[[2]][-1])==value]
lines(fit$estimate[,1],fit$estimate[,columnValue],col="green",lty=1,lwd=2.5)
lines(fit$estimate[,1],fit$estimate[,columnValue]+fit$SD[,columnValue],col="green",lty=2,lwd=2.5)
lines(fit$estimate[,1],fit$estimate[,columnValue]-fit$SD[,columnValue],col="green",lty=2,lwd=2.5)
}

plot(1, type = "n", axes=FALSE, xlab="", ylab="")

legend("center",c("Penalized BTPB","95% CI Using BTPB","Kernel Density","95% CI Using K. Den"),col=c("black","blue","green","green"),lty=c(1,2,1,2), lwd=c(2),cex=1.3)

##Code for Graph 4
par( mfrow = c( 2, 3 ) )
par(mar=c(4.5,4.5,4,1))
for(i in 1:5){
  value<-c(160,180,200,220,240)[i]
  location<-c(1:length(fit$fixedZvalue))[fit$fixedZvalue==value]
  
#plot the kernel density estimation fit
  plot(fit$estimate[,1],fit$estimate[,location+1]+1.96*fit$SD[,location],col="Gray40",lty=3,lwd=3,main=paste0("z=",value),type='l',ylab="",xlab="Yield",ylim=c(0,.05),cex.main=1.75,cex.axis=1.5,cex.lab=1.5)
  lines(fit$estimate[,1],fit$estimate[,location+1]-1.96*fit$SD[,location],col="Gray40",lty=3,lwd=3)
  lines(fit$estimate[,1],fit$estimate[,location+1],type="l",col="Gray50",lty=1,lwd=1)
  
  

}

plot(1, type = "n", axes=FALSE, xlab="", ylab="")

legend("center",c("Kernel Density","95% CI"),col=c("Gray50","gray40"),lty=c(1,3,3), lwd=c(3),cex=1.8)


##Code for graph 5, to incorperate subsidy values from table 1, remove the # (this gets graph 6)

  
  
  par( mfrow = c( 2, 3 ) )
  par(mar=c(4.5,4.5,4,1))
  for(i in 1:5){
    
    value<-c(160,180,200,220,240)[i]
    
    
    x_i<-c(.50,.55,.6,.65,.7,.75,.8,.85)
    mySubsidy<-data.frame(xvalue=c(.50,.55,.60,.65,.70,.75,.80,.85),subsidy=1-c(.67,.64,.64,.59,.59,.55,.48,.38))
    yield<-x_i*value
  
  iadj<-value-5
  iadj2<-value+5
  
  myFinalSub<-merge(data.frame(origX=round(myEst$origX,2),origY=myEst$origY,origZ=myEst$origZ),mySubsidy,by.x="origX",by.y="xvalue",all.x=T)

  
   myFinalSub$Paid<-myFinalSub$origY#*myFinalSub$subsidy
  plot(myFinalSub$origX[myFinalSub$origZ>iadj&myFinalSub$origZ<iadj2],myFinalSub$Paid[myFinalSub$origZ>iadj&myFinalSub$origZ<iadj2],col="black",pch=16,cex=1,main=paste("Land Quality",value,"APH"),ylab="Dollars per Acre",xlab="Coverage Rate", cex.main=1.5, cex.sub=1.5,cex.lab=1.5, cex.axis=1.5)

  lines(myEst$newX[round(myEst$newZ,0)%in% c((value):(value+1))],myEst$Fit[round(myEst$newZ,0)%in%c((value):(value+1))],col="gray50",lwd=3)

  }
  
#Table 2 something similar

##