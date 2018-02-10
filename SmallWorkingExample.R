##Mini test of BTPS

##Create x and Z values
x<-c(0:100)
z<-c(0:100)

##get true values
my.dataframe<-expand.grid(x,z)
y<-with(my.dataframe,Var1^2+Var2^2+Var1^3*Var2)+rnorm(nrow(my.dataframe),0,50)
x<-my.dataframe$Var1
z<-my.dataframe$Var2

#get random sample
set.seed(9389)
randomSample<-sample(10201,size=500)


##fit the BTPS spline
test<-fitBTPS(y[randomSample],x[randomSample],z[randomSample],knots=c(4,4))
true1derv<-2*test$newX+3*test$newX^2*test$newZ 
true2derv<-2+6*test$newX*test$newZ


##Show it fits the derivatives in the easy case
plot3d(c(test$newX,x),c(test$newZ,z),c(test$Fit,y),col=c(rep("Black",length(test$newX)),rep("Red",length(y))))
plot3d(c(test$newX,test$newX),c(test$newZ,test$newZ),c(true1derv,test$firstDerv),col=c(rep("Black",length(test$newX)),rep("Red",length(test$newX))))
plot3d(c(test$newX,test$newX),c(test$newZ,test$newZ),c(true2derv,test$secondDerv),col=c(rep("Black",length(test$newX)),rep("Red",length(test$newX))))


########Mini test of simulated data and the fit (Maybe add in p)
set.seed(9389)
##simulate insurance data
sim.dat<-yieldDataSim(1000,YieldMeanFunc="Linear",YieldVarFunc="NonConst",YieldError="Normal")

##plot the insurance data
plot3d(sim.dat$Cov_Rate,sim.dat$Yield_Hist,sim.dat$Obs_Prem)

##fit the insurance data
test<-fitBTPS(sim.dat$Obs_Prem,sim.dat$Cov_Rate,sim.dat$Yield_Hist,knots=c(2,9),penalty=c(3,3),degree=c(5,3),tol2nd = .0005,newXLim =c(.55,.95) ,newZLim =c(100,300))

##show the values fit
plot3d(c(test$newX,sim.dat$Cov_Rate),c(test$newZ,sim.dat$Yield_Hist),c(test$Fit,sim.dat$Obs_Prem),col=c(rep("Red",length(test$newX)),rep("Black",length(sim.dat$Obs_Prem))))


newOverallSigma<-sim.dat$OverallSigma[1]*abs(-25+1.3*test$newZ)^.1
##fit the second partial derivative
plot3d(c(test$newX,test$newX),c(test$newZ,test$newZ),c(test$secondDerv/test$newZ^2,dnorm((test$newX*test$newZ-(-25+1.3*test$newZ))/newOverallSigma)/newOverallSigma),col=c(rep("Red",length(test$newX)),rep("Black",length(test$newX))))


#sixWayPlot(test$newX,test$newZ,test$secondDerv,(-25+1.3*test$newZ),newOverallSigma)

##Add on the confidence intervals
mine<-VarianceEstimator(test) 

Est<-test$secondDerv/test$newZ^2
SDEst<-sqrt(mine$variance)/test$newZ^2
plot3d(c(test$newX,test$newX,test$newX,test$newX),c(test$newZ,test$newZ,test$newZ,test$newZ),c(Est,dnorm((test$newX*test$newZ-(-25+1.3*test$newZ))/newOverallSigma)/newOverallSigma,Est+1.96*SDEst,Est-1.96*SDEst),col=c(rep("Red",length(test$newX)),rep("Black",length(test$newX)),rep("Green",2*length(test$newX))))

##maybe get the plot working to look
