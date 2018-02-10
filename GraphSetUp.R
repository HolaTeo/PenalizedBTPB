####Graphs for BTPB



myGraphs<-function(simDataEst,simDatVar,quantileV=c(.1,.25,.50,.75,.9),Title=NULL){
  myMeans<-apply(simDataEst[,,-1], c(1,2), mean,na.rm=T)
  myEst<-melt(myMeans)
  myTrue<-melt(simDataEst[,,1])

  

  
  par(mfrow=c(2,3),mar = c(2,1.15,.6,-1) + 2)

  
  roundedZ<-round(myEst[,2],4)
  myZValues<-sort(roundedZ)[round(length(roundedZ)*quantileV)]
  percentile95<-apply(simDataEst[,,-1], c(1,2), quantile,na.rm=T,p=.95)
  percentile5<-apply(simDataEst[,,-1], c(1,2), quantile,na.rm=T,p=.05)
  myEst95<-melt(percentile95)
  myEst5<-melt(percentile5)
  
  for(i in 1:length(quantileV)){
    par(mar = c(2,2,5,1))
    maxValue<-max(myEst95[roundedZ==myZValues[i],3])
    minValue<-min(myEst5[roundedZ==myZValues[i],3])
    
    plot(myEst[,1][roundedZ==myZValues[i]]*myZValues[i],myEst[,3][roundedZ==myZValues[i]],main=paste0(quantileV[i]*100,"th Percentile (z=",round(myZValues[i],-1),")"),type='l',xlab="",col="gray50",ylim=c(minValue,maxValue),ylab="Yield",lty=1,lwd=3,cex.axis=1.3,cex.title=1.5, cex.main=1.5, cex.sub=1.5)

    
    lines(myEst5[,1][roundedZ==myZValues[i]]*myZValues[i],myEst5[,3][roundedZ==myZValues[i]],col="gray60",lty=2,lwd=3)
    
    
    lines(myEst95[,1][roundedZ==myZValues[i]]*myZValues[i],myEst95[,3][roundedZ==myZValues[i]],col="gray60",lty=2,lwd=3)
    
    lines(myTrue[,1][roundedZ==myZValues[i]]*myZValues[i],myTrue[,3][roundedZ==myZValues[i]],col="black",lwd=2)
    
  }
  
  
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  
  legend("center",c("True","MC Mean","2.5 & 97.5\nPercentiles"),col=c("black","gray50","gray60"),lty=c(1,1,2), lwd=c(2,2,2),cex=1.5)
  
  
  
  
  myVar<-apply(simDataEst[,,-1], c(1,2), var,na.rm=T)
  myVarEst<-melt(myVar)
  myInd<-array(,dim = c(100,100,1000))
  for(i in 2:1001){
    
    myInd[,,i-1]<-simDataEst[,,i]+1.96*sqrt(myVar)<simDataEst[,,1]|simDataEst[,,i]-1.96*sqrt(myVar)>simDataEst[,,1]
    
    
  }
  myCoverage<-apply(myInd, c(1,2), mean,na.rm=T)
  
  dimnames(myCoverage)<-dimnames(myMeans)
  
  myCovMelt<-melt(myCoverage)
  
  var1Quant<-sort(myCovMelt$Var1)[c(1,round(length(myCovMelt$Var1)*seq(0,1,length.out=8),0))]
  var2Quant<-sort(myCovMelt$Var2)[c(1,round(length(myCovMelt$Var2)*seq(0,1,length.out=8),0))]
  
  myCovMeltSimp<-subset(myCovMelt,Var1%in%var1Quant&Var2%in%var2Quant)
  names(myCovMeltSimp)[3]<-"Coverage"
  myCovMeltSimp$Coverage<-(1-myCovMeltSimp$Coverage)*100
  
  b <- c(80,85,90,95,100)
  ggplot(myCovMeltSimp, aes(Var1, Var2)) +
    geom_tile(aes(fill = Coverage)) +
    geom_text(size=7,aes(label = (round(Coverage, 1)))) +
    scale_fill_gradientn(limits=c(80,100),
                         colours=c("navyblue", "blue", "light blue","white", "red"),
                         breaks=b, labels=format(b),name="Coverage\nRate")+theme_bw()+ scale_x_continuous(name="Coverage Level") +
    scale_y_continuous(name="Land Quality (APH)")+ theme(text = element_text(size=16) , legend.title=element_text(size=18) , legend.text=element_text(size=14))#+ ggtitle( paste(Title ,"Coverage Rate"))
  
  
  
  
  
}



pdf("HeatMapLinearNorm.pdf",width=8,height=6,useDingbats=F)
myGraphs(myRecall,myRecallVar)

dev.off()  

for(i in 1:21){
  myAPH<- as.numeric(dimnames(QuadNormNon)[[2]])[i]
  MYw<- as.numeric(dimnames(QuadNormNon)[[1]])
  trueQuadNormNon[,i]<- dnorm(MYw,  1.2*(myAPH-50)+(myAPH-150)^2/200,20/2*myAPH^.2)
  
}

###Kernel Heat Maps
myGraphsKern<-function(simDataEst,simDatVar,true,quantileV=c(.1,.25,.50,.75,.9),Title=NULL){
true<-array(,dim=c(1000,21))
for(i in 1:21){
  myAPH<- as.numeric(dimnames(simDataEst)[[2]])[i]
  MYw<- as.numeric(dimnames(simDataEst)[[1]])
}
LowerMatrix<-simDataEst-1.96*simDatVar
UpperMatrix<-simDataEst+1.96*simDatVar

temp<-array(,dim=c(1000,21,1000))
for(j in 1:1000){
  temp[,,j]<-(LowerMatrix[,,j]<true)&(UpperMatrix[,,j]>true)
}

Coverage<-apply(temp,c(1,2),mean)
dimnames(Coverage)<-list(dimnames(simDataEst)[[1]],dimnames(simDataEst)[[2]])

test<-melt(Coverage,id.var="CovRate")
test<-test[test$Var1>(1.2*(test$Var2-50)+(test$Var2-150)^2/200-40*(test$Var2^.2)/2)&test$Var1<(1.2*(test$Var2-50)+(test$Var2-150)^2/200+40*(test$Var2^.2)/2),]
##a StdYield FUnction

test$StdYield<-round(test$Var1-1.2*(test$Var2-50)-(test$Var2-150)^2/200,0)
FinalTest<-ddply(test,.(Var2,StdYield),summarize,covRate=mean(value))
FinalTest<-subset(FinalTest,Var2%in%round(seq(100,300,by = 20),-1)&StdYield%in%round(seq(-40,40,by=8)))
FinalTest<-subset(FinalTest,abs(StdYield)<40)
b <- c(60,70,80,90,95,100)

pdf("KernHeatMapQuadNorm.pdf",width=8,height=6,useDingbats=F)
ggplot(FinalTest, aes((StdYield), (Var2))) +
  geom_tile(aes(fill = covRate*100))+geom_text(size=6,aes(label = (round(covRate*100, 1))))+  scale_fill_gradientn(limits=c(60,100),
                                                                                                                   colours=c("Gray 50","blueviolet", "blue", "light blue","white", "red"),
                                                                                                                   breaks=b, labels=format(b),name="Coverage\nRate")+theme_bw()+ scale_x_continuous(name="Standardized Yield") +
  scale_y_continuous(name="Land Quality (APH)")+ theme(text = element_text(size=16) , legend.title=element_text(size=18) , legend.text=element_text(size=14))
dev.off()  


percentile95<-apply(simDataEst, c(1,2), quantile,na.rm=T,p=.975)
percentile5<-apply(simDataEst, c(1,2), quantile,na.rm=T,p=.025)
myMean<-apply(simDataEst, c(1,2), mean,na.rm=T)


pdf("KernMCQuadBeta.pdf",width=8,height=6,useDingbats=F)
par(mfrow=c(2,3),mar = c(2,0,.6,-1) + 2)
par(mar = c(4,4,4,1))
for(i in 1:5){
  value<-c(120,160,200,240,280)[i]
  percentage<-c(.1,.25,.5,.75,.9)[i]
  where<-c(1:21)[as.numeric(dimnames(simDataEst)[[2]])==value]
  plot(MYw,percentile95[,where],col="Gray50",type="l",lwd=3,lty=2,ylim=c(0,.025),main=paste0(percentage*100,"th Percentile (z=",round(value,-1),")"),xlab="",ylab="",cex.main=1.8,cex.axis=1.5)
  lines(MYw,myMean[,where],col="Gray50",lwd=3,lty=1)
  lines(MYw,percentile5[,where],col="Gray50",lwd=3,lty=2)
  
  lines(MYw,true[,where],lwd=3,lty=1)
}
plot(1, type = "n", axes=FALSE, xlab="", ylab="")

legend("center",c("True","MC Mean","2.5th & 97.5th\nPercentile"),col=c("black","gray50","gray60"),lty=c(1,1,3), lwd=c(2,4,4),cex=1.8)
dev.off()    
}

#####Empirical 

pdf("BTPS_KernelDensity2.pdf",width=8,height=6,useDingbats=F)
par( mfrow = c( 2, 3 ) )
price_p<-4
myZValues<-c(136.4,147,165.4,182.3,192.9)
for(i in 1:5){
  value<-c(135,145,165,180,190,195)[i]
  #plot(form$coverage[value==form$APH_Yield]*value,form$secondD[value==form$APH_Yield]/value^2/price_p,main=value,type='l',ylim=c(0,max(form$secondD,na.rm=T)/value^2))
  
  plot(hi$newX[roundedZ==myZValues[i]]*myZValues[i],hi$secondDerv[roundedZ==myZValues[i]]/myZValues[i]^2/price_p,main=paste0(quantileV[i]*100,"th Percentile (z=",value,")"),type='l',col="black",ylab="",ylim=c(0,max(form$secondD,na.rm=T)/value^2),lwd=2.5,xlab="",cex.main=1.75,cex.axis=1.5,cex.lab=1.5)
  sdValues<-myZValues[i]
  
  
  #lines(hi$newX[roundedZ==myZValues[i]]*myZValues[i],hi$secondDerv[roundedZ==myZValues[i]]/myZValues[i]^2+standardError,main=paste(quantileV[i]*100," Percentile\n Historical Yield=",myZValues[i]),type='l',lty=2)
  roundedZsd<-round(varianceBTPS$zValues,1)
  standardError<-1.96*sqrt(varianceBTPS$variance[roundedZsd==sdValues])*price_p
  lines(varianceBTPS$xValues[roundedZsd==sdValues]*sdValues,hi$secondDerv[roundedZ==sdValues]/sdValues^2/price_p-standardError,type='l',lty=2,col="black",lwd=2.5)
  lines(varianceBTPS$xValues[roundedZsd==sdValues]*sdValues,hi$secondDerv[roundedZ==sdValues]/sdValues^2/price_p+standardError,type='l',lty=2,col="black",lwd=2.5)
  lines(hi2$estimate[,1],hi2$estimate[,i+1],col="Gray50",lty=1,lwd=2.5)
  lines(hi2$estimate[,1],hi2$estimate[,i+1]+1.96*hi2$SD[,i+1],col="Gray50",lty=2,lwd=2.5)
  lines(hi2$estimate[,1],hi2$estimate[,i+1]-1.96*hi2$SD[,i+1],col="Gray50",lty=2,lwd=2.5)
}

plot(1, type = "n", axes=FALSE, xlab="", ylab="")

legend("center",c("Penalized BTPB","95% CI Using BTPB","Kernel Density","95% CI Using K. Den"),col=c("black","Black","gray50","gray60"),lty=c(1,2,1,2), lwd=c(2),cex=1.3)
dev.off()

pdf("True_KernelDensity_GrayDottedThis.pdf",width=8,height=6,useDingbats=F)
par( mfrow = c( 2, 3 ) )
price_p<-4
myZValues<-c(136.4,147,165.4,182.3,192.9)
for(i in 1:5){
  value<-c(135,145,165,180,190,195)[i]
  #plot(form$coverage[value==form$APH_Yield]*value,form$secondD[value==form$APH_Yield]/value^2/price_p,main=value,type='l',ylim=c(0,max(form$secondD,na.rm=T)/value^2))
  
  #plot(hi$newX[roundedZ==myZValues[i]]*myZValues[i],hi$secondDerv[roundedZ==myZValues[i]]/myZValues[i]^2/price_p,main=paste0(quantileV[i]*100,"th Percentile (z=",value,")"),type='l',col="black",ylab="",ylim=c(0,max(form$secondD,na.rm=T)/value^2),lwd=3)
  
  
  plot(hi2$estimate[,1],hi2$estimate[,i+1]+1.96*hi2$SD[,i+1],col="Gray40",lty=3,lwd=3,main=paste0(quantileV[i]*100,"th Percentile (z=",value,")"),type='l',ylab="",xlab="Yield",ylim=c(0,.035),cex.main=1.75,cex.axis=1.5,cex.lab=1.5)
  
  lines(hi2$estimate[,1],hi2$estimate[,i+1]-1.96*hi2$SD[,i+1],col="Gray40",lty=3,lwd=3)
  lines(hi2$estimate[,1],hi2$estimate[,i+1],col="Gray50",lty=1,lwd=1)
}

plot(1, type = "n", axes=FALSE, xlab="", ylab="")

legend("center",c("Kernel Density","95% CI"),col=c("Gray50","gray40"),lty=c(1,3,3), lwd=c(3),cex=1.8)
dev.off()

##Reverse

#################Actual Premium Values
pdf("PremiumValue.pdf",width=8,height=6,useDingbats=F)
estimated<-c()
par( mfrow = c( 2, 3 ) )
par(mar=c(4.5,4.5,4,1))
for(i in c(135,145,165,180,190,195)){
  
  #for(i in seq(130,205,by=5)){
  
  x_i<-c(.50,.55,.6,.65,.7,.75,.8,.85)
  mySubsidy<-data.frame(xvalue=c(.50,.55,.6,.65,.7,.75,.8,.85),subsidy=1-c(.67,.64,.64,.59,.59,.55,.48,.38))
  yield<-x_i*i
  
  yieldRows<-round(as.numeric(dimnames(integralValues)[[1]]))
  
  premium<-c()
  for(k in 1:length(x_i)){
    
    value<-NA
    if(sum(yieldRows==yield[k])>0){
      value<-median(integralValues[yieldRows==yield[k],dimnames(integralValues)[[2]]==i])
    }else if(sum(yieldRows==(yield[k]-.1))>0){
      value<-integralValues[yieldRows==(yield[k]-.1),dimnames(integralValues)[[2]]==i]
    }else if(sum(yieldRows==(yield[k]-.2))>0){value<-integralValues[yieldRows==(yield[k]-.2),dimnames(integralValues)[[2]]==i]
    }else{
      value<-median(integralValues[round(yieldRows)==round(yield[k]),dimnames(integralValues)[[2]]==i])
    }
    premium<-c(premium,as.numeric(value) )
    
  }
  
  
  
  
  iadj<-i-5
  iadj2<-i+5
  
  myFinalSub<-merge(data.frame(origX=hi$origX,origY=hi$origY,origZ=hi$origZ),mySubsidy,by.x="origX",by.y="xvalue",all.x=T)
  myFinalSub$Paid<-myFinalSub$origY#*myFinalSub$subsidy
  plot(myFinalSub$origX[myFinalSub$origZ>iadj&myFinalSub$origZ<iadj2],myFinalSub$Paid[myFinalSub$origZ>iadj&myFinalSub$origZ<iadj2],col="black",ylim=c(0,max(premium,myFinalSub$origY[myFinalSub$origZ>iadj&myFinalSub$origZ<iadj2])),pch=16,cex=1,main=paste("Land Quality",i,"APH"),ylab="Dollars per Acre",xlab="Coverage Rate", cex.main=1.5, cex.sub=1.5,cex.lab=1.5, cex.axis=1.5)
  # lines(hi$origX[hi$origZ==i],hi$fitOrig[hi$origZ==i],col="blue",type="l")
  lines(x_i,premium,col="gray50",lwd=3)
  estimated<-rbind(estimated,data.frame(CovRate=x_i,LandQ=i,Premium=premium))
}
dev.off()  
