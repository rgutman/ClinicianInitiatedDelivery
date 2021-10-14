# This requires mitss_func.R
source("./mitss_func.R")


MITSS = function(treatmentInd, propScore, balanceOn = "all",
                 outcomeVec, dataMat=NA, Resp.Curve="two",
                 estimandFunc=estimandMean, numEstimands=1,
                 isOutcomeBinary=1, nNumImpute=40, nNumSub= 6,
                 obsUse=NA, confidence=0.90){
  
# Find which observation recieved treatment v. control
trtObs = which(treatmentInd == 1)
conObs = which(treatmentInd == 0)

# Define which observations to use.
# Since we left obsUse NA, the else will run.
if(length(obsUse) == 1){
  whichUse = rep(1,length(propScore))	
} else{
  whichUse = rep(0,length(propScore))
  whichUse[obsUse] = 1
}

# Define which observations are used in
# each of the treatment groups.
# Since we have subset the data to the matched set,
# we use the entire data set, just split into treat
# and control.
whichUseTrt = whichUse[trtObs] 
whichUseCon = whichUse[conObs]

#if there is at least one covariate
if(length(dataMat) > 1){
  # removing the orthogonlizing all of the covariates
  # we keep the names to ensure the downstream functions
  # can use the same notation.
  ortMatAll <- dataMat
  ortMatAll <- as.data.frame(ortMatAll)
  #defining the covariate matrices of the control/treatment group
  ortTrtMat = ortMatAll[trtObs,]
  ortConMat = ortMatAll[conObs,]
    }else{
  ortTrtMat = NA
  ortConMat = NA}

#propensity score values of all the subjects in the treated group
curXtrt2 = propScore[trtObs]

#propensity score values of all the subjects in the control group
curXcon2 = propScore[conObs]

#outcome values of all the subjects in the treatment group
curYtrt2 = subset(outcomeVec, treatmentInd == 1)

#outcome values of all the subjects in the control group
curYcon2 = subset(outcomeVec, treatmentInd == 0)

#finding the boundaries of the propensity score values
# in the treatment group and the control group
minC = min(curXcon2)
minT = min(curXtrt2) # This was mislabeled.
maxC = max(curXcon2)
maxT = max(curXtrt2) # This was mislabeled.

#If obsUse is empty using all of the observations.
# This won't run since obsUse is NA.
if(length(obsUse) == 1){
  obsUse = seq(1:length(propScore))
}

# Defining the breakpoints/knots on the propensity score
# the following lines are done in order to make sure
# that we have at least one treated and one control within
# each propensity score subclass starting with the number
# that was defined nNumSub.
# * Note that nNumSub is defined as 6 in the beginning.
if(balanceOn == "all"){
  subClassVecR = cut(propScore[obsUse],
                     breaks = quantile(propScore[obsUse],
                                       prob = (0:nNumSub)/nNumSub)-c(0.01,rep(0,nNumSub)), labels=FALSE)
}else if(balanceOn == "treatment"){
  # treat.obs <- obsUse[exam[obsUse,]$noins==1]
  subClassVecR = cut(propScore[obsUse],breaks = quantile(propScore[trtObs], prob = (0:nNumSub)/nNumSub)-c(Inf,rep(0,nNumSub-1),-Inf), labels=FALSE)
}else if(balanceOn == "control"){
  control.obs <- obsUse[exam[obsUse,]$noins==0]
  subClassVecR = cut(propScore[obsUse],breaks = quantile(propScore[control.obs], prob = (0:nNumSub)/nNumSub)-c(Inf,rep(0,nNumSub-1),-Inf), labels=FALSE)
}

curQuantR = nNumSub

# Finding how many subjects recieve the treatment in each subclass
# and how many recieve the control in each subclass.
countVec <- tapply(treatmentInd[obsUse], subClassVecR, sum) # Treat
countVecL <- tapply((1 - treatmentInd[obsUse]),
                    subClassVecR, sum) # Control

# As long as there is one Subclass where there is no subjects
# in either the control or treated group continue.
# Running this line as is reduced the knots/breaks to 1.
# I expect it was meant to work until each subclass had at least
# one subject which it did. So I will change this statement
# to NOT (!).
# (1a) The first truth statement seeks for more than 3 in each
# subclass for treatment ) OR
# (1b) The second truth statement seeks for more than 3 in each
# subclass for control.
# (2) The third truth statement ensure we still have at least one
# break point at the end of this iteration.
while((length(which(countVec < 3)) > 0  |
       length(which(countVecL < 3)) > 0) & curQuantR > 1){
  # Reduce the number of subclasses
  curQuantR<-curQuantR - 1 
  
  #recalculate for each unit the subclass that it belong to
  if(curQuantR > 1){
    if(balanceOn == "all"){
      subClassVecR = cut(propScore[obsUse],
                         breaks = quantile(propScore[obsUse],
                                           prob = (0:curQuantR)/curQuantR)-c(0.01,rep(0,curQuantR)), labels=FALSE)
    }else if(balanceOn == "treatment"){
     # treat.obs <- obsUse[exam[obsUse,]$noins==1]
      subClassVecR = cut(propScore[obsUse],breaks = quantile(propScore[trtObs], prob = (0:curQuantR)/curQuantR)-c(Inf,rep(0,curQuantR-1),-Inf), labels=FALSE)
    }else if(balanceOn == "control"){
      control.obs <- obsUse[exam[obsUse,]$noins==0]
      subClassVecR = cut(propScore[obsUse],breaks = quantile(propScore[control.obs], prob = (0:curQuantR)/curQuantR)-c(Inf,rep(0,curQuantR-1),-Inf), labels=FALSE)
    }
  } else{
    subClassVecR<-rep(1, length(propScore[obsUse]))}
  
  #finding how many subjects recieve the treatment in each subclass and how many recieve the control in each subclass
  countVec<-tapply(treatmentInd[obsUse], subClassVecR, sum);
  countVecL<-tapply((1 - treatmentInd[obsUse]), subClassVecR, sum);
}

# Defining the knots of the spline based on the subclass	
if(balanceOn == "all"){
  knots = quantile(propScore[obsUse], prob = (0:curQuantR)/curQuantR)-c(0.01,rep(0,curQuantR))
}else if(balanceOn == "treatment"){
  # treat.obs <- obsUse[exam[obsUse,]$noins==1]
  knots = quantile(propScore[trtObs], prob = (0:curQuantR)/curQuantR)-c(0.01,rep(0,curQuantR))
}else if(balanceOn == "control"){
  control.obs <- obsUse[exam[obsUse,]$noins==0]
  knots = quantile(propScore[control.obs], prob = (0:curQuantR)/curQuantR)-c(0.01,rep(0,curQuantR))
}	

nNumKnots = length(knots)

# The choice of using one response curve or 2 response curves
# (1 treatment & 1 Control)
if(Resp.Curve!="one"){
  knotsT = knots
  knotsC = knots
  
  # setting the boundary knots of the spline for the
  # treatment group at the maximal and minimal values of
  # the propensity score value for the treatment group
  knotsT[1] = min(curXtrt2)
  knotsT[nNumKnots] = max(curXtrt2)
  
  # setting the boundary knots of the spline for the control
  # group at the maximal and minimal values of the propensity
  # score value for the control gorup
  knotsC[1] = min(curXcon2)
  knotsC[nNumKnots] = max(curXcon2)
  
  
  # Creating the cubic spline basis over the propensity score
  # for the treated units using the knots defined earlier
  # and adding the ortogonlized covariates if they exist
  lenKnotsT = length(knotsT)
  if(length(knots)>2){
    basisT = ns(curXtrt2, knots = knotsT[2:(lenKnotsT -1)])
  } else{
    basisT = ns(curXtrt2)
    
  }
  
  # * Line above did not originally run when we allowed the
  # subclasses to reduce to one.
  # basisT = ns(curXtrt2)

  if(dim(as.matrix(ortTrtMat))[2] > 1){
    xxT = model.matrix(~ as.matrix(predict(basisT, curXtrt2) + ortTrtMat));
  } else{
    xxT = model.matrix(~ predict(basisT, curXtrt2));
  }
  
  # Creating the cubic spline basis over the propensity score
  # for the control units using the knots defined earlier
  # and adding the ortogonlized covariates if the exists
  lenKnotsC = length(knotsC)
  
  if(length(knots)>2){
    basisC = ns(curXcon2, knots = knotsC[2:(lenKnotsC -1)])
  } else{
    basisC = ns(curXcon2)
    
  }
  
  # * Line above also did not run when we allowed the algorithm
  # to iteratively reduce the subclasses to 1.
  # basisC = ns(curXcon2)
  
  if(dim(as.matrix(ortConMat))[2] > 1){
    xxC = model.matrix(~ as.matrix(predict(basisC, curXcon2) + ortConMat));
  } else{
    xxC = model.matrix(~ predict(basisC, curXcon2));
  }
  
  # Calculating outcome model for the treatment units.
  # * This is the point where Roee thougt we might see
  # the issues.
  if(isOutcomeBinary == 1){
    #sampling coefficients for the treatment model
    bValTrtAS = sampleCoeffB(curYtrt2, xxT, nNumImpute)
    
    #sampling coefficients for the control model
    bValConAS = sampleCoeffB(curYcon2, xxC, nNumImpute)
    }else{
    #sampling coefficients for the treatment model
    bValTrtAS = sampleCoeffC(curYtrt2, xxT, nNumImpute)
    
    #sampling coefficients for the control model
    bValConAS = sampleCoeffC(curYcon2, xxC, nNumImpute)
    }
  }else{
  knotsAll = knots
  curXall = c(curXcon2,curXtrt2) 
  curYall = c(curYcon2,curYtrt2)
  ortMat = as.matrix(rbind(ortTrtMat,ortConMat))
  
  #setting the boundary knots of the spline 
  knotsAll[1] = min(curXall)
  knotsAll[nNumKnots] = max(curXall)
  
  #creating the cubic spline basis over the propensity score
  #and adding the ortogonlized covariates if they exist
  lenKnotsAll = length(knotsAll)
  basisAll = ns(curXall, knots = knotsAll[2:(lenKnotsAll -1)])
  if(dim(as.matrix(ortMat))[2] > 1){
    xxAll = model.matrix(~ predict(basisAll, curXall) + ortMat);
  } else{
    xxAll = model.matrix(~ predict(basisAll, curXall));
  }
  
  #calculating outcome model for all units
  if(isOutcomeBinary == 1){
    #sampling coefficients for model
    bValAllAS = sampleCoeffB(curYall, xxAll, nNumImpute)
  } else{
    #sampling coefficients for the treatment model
    bValAllAS = sampleCoeffC(curYall, xxAll, nNumImpute)
  }
}

#Matrix to save point estimate from each imputation
imputeMatASM = matrix(0,nNumImpute,numEstimands)
imputeMatASM.ATT = matrix(0,nNumImpute,numEstimands)
imputeMatASM.ATC = matrix(0,nNumImpute,numEstimands)

#Matrix to save covariance matrix from each imputation
imputeMatASSig = array(0,dim=c(nNumImpute,numEstimands,numEstimands))
imputeMatASSig.ATT = array(0,dim=c(nNumImpute,numEstimands,numEstimands))
imputeMatASSig.ATC = array(0,dim=c(nNumImpute,numEstimands,numEstimands))

#estimating the causal effect only on the population for each overlap exists
obsUseN = which(c(whichUseCon,whichUseTrt) == 1)
obsUseN.ATT = which(c(rep(0,length(whichUseCon)),whichUseTrt) == 1)
#obsUseN.ATT = trtObs
obsUseN.ATC = which(c(whichUseCon,rep(0,length(whichUseTrt))) == 1)
#obsUseN.ATC = conObs

#the imputation phase
if(Resp.Curve!="one"){
  #Imputation phase using two response curve models (one control, one treatment)
  for(numImpute in 1:nNumImpute){
    if(isOutcomeBinary == 1){
      #obtaining the linear predictor of the outcome for the predicted treated potential outcome for those that were in the control
      prob = estimationPredict(bValTrtAS$coef[,numImpute], curXcon2, ortConMat, basisT)
      #sampling the outcome from binary distribution with linear predictors 
      imputeTrt= rbinom(length(prob), 1, 1/(1+exp(-prob)))
      
      #obtaining the linear predictor of the outcome for the predicted control potential outcome for those that were in the treated
      prob = estimationPredict(bValConAS$coef[,numImpute], curXtrt2, ortTrtMat, basisC)
      #sampling the outcome from binary distribution with linear predictors 
      imputeCon= rbinom(length(prob), 1, 1/(1+exp(-prob)))
      } else{
      #if each potential outcome is scalar
      if(bValTrtAS$numOut == 1){
        #the standard error of the parameters is just a scalar
        sigMatT = matrix(bValTrtAS$sig[nNumImpute],bValTrtAS$numOut,bValTrtAS$numOut)
        sigMatC = matrix(bValConAS$sig[nNumImpute],bValConAS$numOut,bValConAS$numOut)
      } else{
        #the standard error of the parameters is just a matrix
        sigMatT = matrix(bValTrtAS$sig[,nNumImpute],bValTrtAS$numOut,bValTrtAS$numOut)
        sigMatC = matrix(bValConAS$sig[,nNumImpute],bValConAS$numOut,bValConAS$numOut)
        }
      
      #obtaining the linear predictor of the outcome for the predicted treated potential outcome for those that were in the control 
      #using multivariate normal distribution
      randMatT = matrix(rnorm(length(curXcon2)*bValTrtAS$numOut, 0, 1),length(curXcon2),bValTrtAS$numOut)
      imputeTrt= estimationPredict(matrix(bValTrtAS$coef[,numImpute],bValTrtAS$numCov,bValTrtAS$numOut), curXcon2, ortConMat, basisT) + randMatT%*%chol(sigMatT)
      
      
      #obtaining the linear predictor of the outcome for the predicted control potential outcome for those that were in the treated
      #using multivariate normal distribution
      randMatC = matrix(rnorm(length(curXtrt2)*bValConAS$numOut, 0, 1),length(curXtrt2),bValConAS$numOut)
      imputeCon= estimationPredict(matrix(bValConAS$coef[,numImpute],bValConAS$numCov,bValConAS$numOut), curXtrt2, ortTrtMat, basisC) + randMatC%*%chol(sigMatC)
    }
    
    #using the estimand function to calculate the estimands on the observations that should be included
    if(bValTrtAS$numOut == 1)
    {
      estimVal = estimandFunc(as.matrix(c(imputeTrt,curYtrt2)), as.matrix(c(curYcon2,imputeCon)),obsUseN)
      estimVal.ATT = estimandFunc(as.matrix(c(imputeTrt,curYtrt2)), as.matrix(c(curYcon2,imputeCon)),obsUseN.ATT)
      estimVal.ATC = estimandFunc(as.matrix(c(imputeTrt,curYtrt2)), as.matrix(c(curYcon2,imputeCon)),obsUseN.ATC)
    } else{				
      estimVal = estimandFunc(rbind(imputeTrt,curYtrt2), rbind(curYcon2,imputeCon),obsUseN)
      estimVal.ATT = estimandFunc(rbind(imputeTrt,curYtrt2), rbind(curYcon2,imputeCon),obsUseN.ATT)
      estimVal.ATC = estimandFunc(rbind(imputeTrt,curYtrt2), rbind(curYcon2,imputeCon),obsUseN.ATC)
    }
    
    #Saving the point estimate for this imputation
    imputeMatASM[numImpute, ]<- estimVal$Q
    imputeMatASM.ATT[numImpute, ]<- estimVal.ATT$Q
    imputeMatASM.ATC[numImpute, ]<- estimVal.ATC$Q
    
    #saving the variance covariance matrix for this imputation
    imputeMatASSig[numImpute, ,]<-estimVal$U
    imputeMatASSig.ATT[numImpute, ,]<-estimVal.ATT$U
    imputeMatASSig.ATC[numImpute, ,]<-estimVal.ATC$U
  }
                }else{
  #Inputation phase using one response curve
  for(numImpute in 1:nNumImpute)
  {
    if(isOutcomeBinary == 1)
    {
      #obtaining the linear predictor of the outcome for the potential outcomes
      prob = estimationPredict(bValAllAS$coef[,numImpute], curXall, ortMat, basisAll)
      #sampling the outcome from binary distribution with linear predictors 
      imputeAll= rbinom(length(prob), 1, 1/(1+exp(-prob)))
      imputeTrt = imputeAll[1:length(curYcon2)]
      imputeCon = imputeAll[(length(curYcon2)+1):length(c(curYcon2,curYtrt2))]
    } else {
      
      #if each potential outcome is scalar
      if(bValAllAS$numOut == 1)
      {
        #the standard error of the parameters is just a scalar
        sigMatAll = matrix(bValAllAS$sig[nNumImpute],bValAllAS$numOut,bValAllAS$numOut)
      } else
      {
        #the standard error of the parameters is just a matrix
        sigMatAll = matrix(bValAllAS$sig[,nNumImpute],bValAllAS$numOut,bValAllAS$numOut)
      }
      
      #obtaining the linear predictor of the outcome for potential outcomes 
      #using multivariate normal distribution
      randMatAll = matrix(rnorm(length(curXall)*bValAllAS$numOut, 0, 1),length(curXall),bValAllAS$numOut)
      imputeAll = estimationPredict(matrix(bValAllAS$coef[,numImpute],bValAllAS$numCov,bValAllAS$numOut), curXall, ortMat, basisAll) + randMatAll%*%chol(sigMatAll)
      imputeTrt = imputeAll[1:length(curYcon2)]
      imputeCon = imputeAll[(length(curYcon2)+1):length(c(curYcon2,curYtrt2))]
    }
    
 
    #using the estimand function to calculate the estimands on the observations that should be included
    if(bValAllAS$numOut == 1){
      estimVal = estimandFunc(as.matrix(c(imputeTrt,curYtrt2)), as.matrix(c(curYcon2,imputeCon)),obsUseN)
      estimVal.ATT = estimandFunc(as.matrix(c(imputeTrt,curYtrt2)), as.matrix(c(curYcon2,imputeCon)),obsUseN.ATT)
      estimVal.ATC = estimandFunc(as.matrix(c(imputeTrt,curYtrt2)), as.matrix(c(curYcon2,imputeCon)),obsUseN.ATC)
    } else{				
      estimVal = estimandFunc(rbind(imputeTrt,curYtrt2), rbind(curYcon2,imputeCon),obsUseN)
      estimVal.ATT = estimandFunc(rbind(imputeTrt,curYtrt2), rbind(curYcon2,imputeCon),obsUseN.ATT)
      estimVal.ATC = estimandFunc(rbind(imputeTrt,curYtrt2), rbind(curYcon2,imputeCon),obsUseN.ATC)
          }
    
    #Saving the point estimate for this imputation
    imputeMatASM[numImpute, ]<- estimVal$Q
    imputeMatASM.ATT[numImpute, ]<- estimVal.ATT$Q
    imputeMatASM.ATC[numImpute, ]<- estimVal.ATC$Q
    
    #saving the variance covariance matrix for this imputation
    imputeMatASSig[numImpute, ,]<-estimVal$U
    imputeMatASSig.ATT[numImpute, ,]<-estimVal.ATT$U
    imputeMatASSig.ATC[numImpute, ,]<-estimVal.ATC$U
  }
}

#calculating the point estimate across imputations
meanVal = apply(imputeMatASM,2,mean)
meanVal.ATT = apply(imputeMatASM.ATT,2,mean)
meanVal.ATC = apply(imputeMatASM.ATC,2,mean)

#calculating between imputation variance
Bm = var(imputeMatASM)
Bm.ATT = var(imputeMatASM.ATT)
Bm.ATC = var(imputeMatASM.ATC)

#calculating within imputation variance
Tm = apply(imputeMatASSig,2:3,mean)
Tm.ATT = apply(imputeMatASSig.ATT,2:3,mean)
Tm.ATC = apply(imputeMatASSig.ATC,2:3,mean)

#calculating super population variance
varAll = Tm + (1 + 1/nNumImpute)*Bm
varAll.ATT = Tm.ATT + (1 + 1/nNumImpute)*Bm.ATT
varAll.ATC = Tm.ATC + (1 + 1/nNumImpute)*Bm.ATC

#calculating finite population variance estimate
varFin = (1 + 1/nNumImpute)*Bm
varFin.ATT = (1 + 1/nNumImpute)*Bm.ATT
varFin.ATC = (1 + 1/nNumImpute)*Bm.ATC

#degrees of freedom for T distribution of multiple imputation (super population)
tValDf = (nNumImpute - 1) * ((1 + 1/nNumImpute) * crossprod(as.vector(Bm), as.vector(t(solve(varAll))))/numEstimands)^-2
tValDf.ATT = (nNumImpute - 1) * ((1 + 1/nNumImpute) * crossprod(as.vector(Bm.ATT), as.vector(t(solve(varAll.ATT))))/numEstimands)^-2
tValDf.ATC = (nNumImpute - 1) * ((1 + 1/nNumImpute) * crossprod(as.vector(Bm.ATC), as.vector(t(solve(varAll.ATC))))/numEstimands)^-2

#degrees of freedom for T distribution of multiple imputation (finite population)
fin.tValDf = (nNumImpute - 1) * ((1 + 1/nNumImpute) * crossprod(as.vector(Bm), as.vector(t(solve(varFin))))/numEstimands)^-2
fin.tValDf.ATT = (nNumImpute - 1) * ((1 + 1/nNumImpute) * crossprod(as.vector(Bm.ATT), as.vector(t(solve(varFin.ATT))))/numEstimands)^-2
fin.tValDf.ATC = (nNumImpute - 1) * ((1 + 1/nNumImpute) * crossprod(as.vector(Bm.ATC), as.vector(t(solve(varFin.ATC))))/numEstimands)^-2

#Margin of error for confidence intervals (super population)
M = qt(.5*(1+confidence),tValDf)*sqrt(varAll)
M.ATT = qt(.5*(1+confidence),tValDf.ATT)*sqrt(varAll.ATT)
M.ATC = qt(.5*(1+confidence),tValDf.ATT)*sqrt(varAll.ATC)

#Confidence intervals (super poulation)
conit = c(meanVal-M,meanVal+M)
conit.ATT = c(meanVal.ATT-M.ATT,meanVal.ATT+M.ATT)
conit.ATC = c(meanVal.ATC-M.ATC,meanVal.ATC+M.ATC)

#Margin of error for confidence intervals (finite population)
fin.M = qt(.5*(1+confidence),fin.tValDf)*sqrt(varFin)
fin.M.ATT = qt(.5*(1+confidence),fin.tValDf.ATT)*sqrt(varFin.ATT)
fin.M.ATC = qt(.5*(1+confidence),fin.tValDf.ATT)*sqrt(varFin.ATC)

#Confidence intervals (finte poulation)
fin.conit = c(meanVal-fin.M,meanVal+fin.M)
fin.conit.ATT = c(meanVal.ATT-fin.M.ATT,meanVal.ATT+fin.M.ATT)
fin.conit.ATC = c(meanVal.ATC-fin.M.ATC,meanVal.ATC+fin.M.ATC)

output<-matrix(NA,nrow=3,ncol=8)
output[,1]<-c(meanVal,meanVal.ATT,meanVal.ATC)
output[,2]<-c(varAll,varAll.ATT,varAll.ATC)
output[,3]<-sqrt(c(varAll,varAll.ATT,varAll.ATC))
output[,4]<-round(output[,1]/output[,3],digits=4)
output[,5]<-c(length(obsUseN),length(obsUseN.ATT),length(obsUseN.ATC))
output[,6]<-round(c(tValDf,tValDf.ATT,tValDf.ATC),digits=2)
output[,7]<-c(conit[1],conit.ATT[1],conit.ATC[1])
output[,8]<-c(conit[2],conit.ATT[2],conit.ATC[2])

row.names(output)<-c("ATE","ATT","ATC")
colnames(output)<-c("estimate","variance","std.error","t.stat","size","deg.F",
                    paste(confidence*100,".LB",sep=""),
                    paste(confidence*100,".UB",sep=""))

output2<-matrix(NA,nrow=3,ncol=8)
output2[,1]<-c(meanVal,meanVal.ATT,meanVal.ATC)
output2[,2]<-c(varFin,varFin.ATT,varFin.ATC)
output2[,3]<-sqrt(c(varFin,varFin.ATT,varFin.ATC))
output2[,4]<-round(output2[,1]/output2[,3],digits=4)
output2[,5]<-c(length(obsUseN),length(obsUseN.ATT),length(obsUseN.ATC))
output2[,6]<-round(c(fin.tValDf,fin.tValDf.ATT,fin.tValDf.ATC),digits=2)
output2[,7]<-c(fin.conit[1],fin.conit.ATT[1],fin.conit.ATC[1])
output2[,8]<-c(fin.conit[2],fin.conit.ATT[2],fin.conit.ATC[2])

row.names(output2)<-c("ATE","ATT","ATC")
colnames(output2)<-c("estimate","variance","std.error","t.stat","size","deg.F",
                     paste(confidence*100,".LB",sep=""),
                     paste(confidence*100,".UB",sep=""))
output3<-matrix(NA,nrow=3,ncol=length(countVec)+1)
output3[1,]<-c(countVec,sum(countVec))
output3[2,]<-c(countVecL,sum(countVecL))
output3[3,]<-c(countVec+countVecL,sum(countVec+countVecL))

colnames(output3)<-c(paste("subclass",1:length(countVec)),"Total")
rownames(output3)<-c("Treatment","Control","Total")


if(length(knots)>2){
  return(list(super.pop=output,finite.pop=output2,SubClass.count=output3, subclass_vec = subClassVecR))
} else{
  return(list(super.pop=output,finite.pop=output2,SubClass.count=output3, knots = "No knots", subclass_vec = subClassVecR))
  
  }

}

