
##Now the sim Study
simStudyFunc<-function(reps=1000,seed=9389,samplesize=1500,Yield_Mean="Linear",VarFunc="Const",Yield_Error="Normal"  ,alphaVal=5,betaVal=3,pVal=4,MeasurementError=.001,BTPBknots=c(2,9),BTPBpenalty=c(3,3),BTPBdegree=c(5,3),tol2nd = .00005,XLim =c(.55,.95) ,ZLim =c(100,300)){
  
  set.seed(seed)
  
  #Create a list of the names
  repNames<-list("True")
 
  for(j in 2:(reps+1)){
    repNames[j]<-paste("Rep",(j-1),sep="")
    
  }
  
 
  Estimates<-array(,dim=c(100,100,(reps+1)))
  SDEst<-array(,dim=c(100,100,reps))
  

  
  for(reps in 2:1001){
    sim.dat<-yieldDataSim(samplesize,YieldMeanFunc=Yield_Mean,YieldVarFunc=VarFunc,YieldError=Yield_Error,alphaE=alphaVal,betaE=betaVal,p=pVal,meas_error=MeasurementError)
    test<-fitBTPS(sim.dat$Obs_Prem,sim.dat$Cov_Rate,sim.dat$Yield_Hist,knots=BTPBknots,penalty=BTPBpenalty,degree=BTPBdegree,tol2nd = tol2nd,newXLim =XLim ,newZLim =ZLim)
    
    mine<-VarianceEstimator(test) 
    overallSig<-rep(sim.dat$OverallSigma[1],length(test$newX))
    Estimates[,,reps]<-matrix(test$secondDerv/(pVal*test$newZ^2),nrow=100)
    SDEst[,,reps-1]<-matrix(sqrt(mine$variance)/(pVal*test$newZ^2),nrow=100) 
    
    
    
    print(reps)
    
    
  }
  dimnames(Estimates)<-list(unique(test$newX),unique(test$newZ),repNames)
  dimnames(SDEst)<-list(unique(test$newX),unique(test$newZ),repNames[-1])
  
  #True Values (Depends on Set Up)
  {
    
    if(Yield_Mean=="Linear"){
      reg_Yield<-(-25+1.3*test$newZ)

    }
    
    if(Yield_Mean=="Quad"){
      reg_Yield<-(1.2*(test$newZ-50)+(test$newZ-150)^2/200)
    }
    
    if(YieldError=="Normal"){
    if(VarFunc=="Const"){
      overallSig<-rep(10,length(test$newX))
      
    }
      if(VarFunc=="NonConst"){
        overallSig<-10*abs(reg_Yield)^.1
        
      }  
      Estimates[,,1]<-matrix(dnorm((test$newX*test$newZ-reg_Yield)/overallSig)/overallSig,nrow=100) 
    }
    
    if(YieldError=="Beta"){
      if(VarFunc=="Const"){
        overallSig<-rep(50,length(test$newX))
        
      }
      if(VarFunc=="NonConst"){
        overallSig<-50*abs(reg_Yield)^.1
        
      }  
      Estimates[,,1]<-matrix(dbeta((test$newX*test$newZ-reg_Yield)/overallSig+(alphaVal/(alphaVal+betaVal)),alphaVal,betaVal)/overallSig,nrow=100)
      
    }
    
    
  #overallSig<-rep(55,length(test$newX))
  Estimates[,,1]<-matrix(dnorm((test$newX*test$newZ-(-25+1.3*test$newZ))/overallSig)/overallSig,nrow=100)
  

  
  
}

  
  return(list(Estimates=Estimates,StdError=SDEst))
}




#Now the sim study for Kernel
simStudyKernel<-function(Nreps=1000,seed=9389,zValues=seq(100,300, by=10),wRange=c(0,450)){
  set.seed(seed)
  
  repNames<-c()

  for(j in 1:Nreps){
    repNames[j]<-paste("Rep",(j),sep="")
    
  }
  fixedZvalue<-zValues
  numberFixed<-length(zValues)
  beginW<-wRange[1]
  endW<-wRange[2]
  Estimates<-array(,dim=c(1000,numberFixed,1000),dimnames=list(seq(beginW,endW,length.out=1000),fixedZvalue,repNames))
  VarEst<-array(,dim=c(1000,numberFixed,1000),dimnames=list(seq(beginW,endW,length.out=1000),fixedZvalue,repNames))

  for(reps in 1:Nreps){
    
    ###now just need to set up the data sim the same
    sim.dat<-yieldDataSim(1500,YieldMeanFunc="Quad",YieldError="Beta",YieldVarFunc = "Prop")
    
    fit<-fitKernel(sim.dat$Yield_Curr,sim.dat$Yield_Hist,knots=10,wRange=wRange,fixedZ = fixedZvalue)
    

    Estimates[,,reps]<-fit[[1]][,-1]
    VarEst[,,reps]<-fit[[2]][,-1]

  }
 return(list(Estimates=Estimates,Variance=VarEst)) 
}





