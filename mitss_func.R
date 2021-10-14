# mitss_func.R is a file that contains all functions that are called upon
# when running mitss.R


# Function that checks if the column is a factor column or binary
IsBinaryOrFactor <- function(column){
  if(length(unique(column)) == 2){
    return(1)}
  if(is.factor(column) == TRUE){
    return(1)}
  return(0)
}

#StepwisePropCalc is a function to estimate the propensity score based on Imbens & Rubin chapter 13 stepwise procedure
#CovMat is a N by K matrix with N people and K covariates. MustInclude is a vecor of size K of {0,1} values that defines
#covariates that must be included in the propensity score. Cl is scalar value defining the threshhold for including main effect
#covariate. Cq is a scalar value defining the threshold for including an interaction of covariates
StepwisePropCalc = function(TreatIndicator, CovMat, MustInclude, Cl=1, Cq=2.71, main.only = FALSE)
{
  #saving the name of the covariates
  colNames = names(CovMat)
  
  #CovMat = as.matrix(CovMat)
  
  #addTem is a function that adds a covaraiate/interaction to the curren model. termtoAdd is the column of the covariate/interaction considered.
  #curModel is the model in the current state. matName is the name of the matrix of covariates
  addTerm = function(termToAdd, curModel, matName)
  {
    #adding a new covariate to the string    
    strAddFormula = sprintf("~.+%s[,%d]",matName, termToAdd)
    
    #generating and R formula object from the new string
    fadd = formula(strAddFormula)
    
    #generating a new model with the added covariate/interaction
    fadd = update(formula(curModel), fadd)
    
    #comparing the previous model to the new model using the Chi^2 test       
    aovVal = add1(curModel, fadd, scale = 0, test ="Chisq", k = 0)
    
    #returing the Chi^2 value of the comparison between new model to old model
    return(aovVal[2,4])
  }
  
  #creating a formula with all of the covariates that must be in the model
  if(sum(MustInclude) >0)
  {
    strFormula = sprintf("TreatIndicator~%s]",paste("CovMat[,", which(MustInclude == 1), sep = "",collapse="]+"))
  }
  else
  {
    strFormula = sprintf("TreatIndicator~1")
  }
  
  #checking the formula
  
  #Calculating propensity score model for all of the covariates that must be in the model
  curModel = glm(formula(strFormula),family=binomial(logit))
  
  #checking which of the covariates was not included in the model yet
  notYetIncluded = which(MustInclude == 0)
  
  #for of the covariates not included in the model find their Chi^2 value for adding the covariate
  addRes = sapply(notYetIncluded, addTerm, curModel=curModel, matName = "CovMat")
  
  #find the covariate with the highest Chi^2 value
  maxAdd = which.max(addRes)
  
  #as long as there are covariates that have Chi^2 value greater than threshold Cl 
  #add covariate one at a time  	
  while(as.vector(addRes)[maxAdd] >= Cl)
  {
    #update the old formula with additional covariate over the threshold
    strAddFormula = sprintf("~.+CovMat[,%d]",notYetIncluded[maxAdd])
    
    #generate a new R formula object
    fadd = formula(strAddFormula)
    
    #update the old model with the new model
    fadd = update(formula(curModel), fadd)  
    
    #calculate new propensity score model   	
    curModel = update(curModel, fadd)
    
    #change the status of the added covariate to be included in the model
    MustInclude[notYetIncluded[maxAdd]] = 1
    notYetIncluded = which(MustInclude == 0)
    
    #for all of the covarites not yet incuded in the model calculate the Chi^2 value to add
    #to the current model
    if(length(notYetIncluded) > 0)
    {
      addRes = sapply(notYetIncluded, addTerm, curModel=curModel, matName = "CovMat")
      maxAdd = which.max(addRes)
    } else
    {
      break;
    }
  }
  
  #after adding all of the main effect covariates. Now working on interactions and sqaured main effects
  
  #generating a formula to calculate all of the second order interactions
  strFormula = sprintf("~(%s])^2-1",paste("CovMat[,", which(MustInclude == 1), sep = "",collapse="]+"))
  
  
  #checking which columns are continuous so that their sqaured will be examined
  isBinaryOrFactor = lapply(CovMat,IsBinaryOrFactor)
  #generating a formula to calculate all of the sqaured main effect covariates
  strFormula = sprintf("%s+%s]^2)",strFormula,paste("I(CovMat[,", which(MustInclude == 1 & isBinaryOrFactor == 0), sep = "",collapse="]^2)+"))
  
  #removing the main effects from the formula
  strFormula = sprintf("%s-%s]",strFormula,paste("CovMat[,", which(MustInclude == 1), sep = "",collapse="]-"))
  
  #generating a amtrix that only includes second order interactions and squared main effects
  CovMatSec = model.matrix(formula(strFormula))
  
  #creating temp dataset that mimicks covariates names
  strFormula2 = sprintf("~(%s)^2-1",paste(colNames[MustInclude==1], sep = "",collapse="+"))
  strFormula2 = sprintf("%s+%s^2)",strFormula2,paste("I(",colNames[MustInclude==1 & isBinaryOrFactor == 0], sep = "",collapse="^2)+"))
  strFormula2 = sprintf("%s-%s",strFormula2,paste(colNames[MustInclude==1], sep = "",collapse="-"))
  temDat = data.frame(CovMat[1:2,])
  colnames(temDat) = colNames
  temDat = model.matrix(formula(strFormula2),data=temDat)
  #end temp dataset
  
  #all of the second order and sqaured main effects are not included at first
  shouldInclude = rep(0, length(CovMatSec[1,]))
  notYetIncluded = which(shouldInclude == 0)
  
  #calculating the Chi^2 value for adding each of the interaction or squared effect
  addRes = sapply(notYetIncluded, addTerm, curModel=curModel, matName = "CovMatSec")
  
  #finding the max Chi^2 value among the interactions and sqaured effect
  maxAdd = which.max(addRes)
  
  if(main.only==FALSE)
  {
    #Continue to add interactions or sqaured effects in a stepwise fashion as long as they are smaller than the 
    #threshold Cq
    while(addRes[maxAdd] >= Cq)
    {
      #update the propensity score model with the new interactio/squared effect
      strAddFormula = sprintf("~.+CovMatSec[,%d]",notYetIncluded[maxAdd])
      fadd = formula(strAddFormula)
      fadd = update(formula(curModel), fadd)     	
      curModel = update(curModel, fadd)
      
      #change the status of the interaction/sqaured effect to included
      shouldInclude[notYetIncluded[maxAdd]] = 1
      
      #now calculate the Chi^2 for the interactions/sqaured effects that were not included yet in the propesnity score
      notYetIncluded = which(shouldInclude == 0)
      addRes = sapply(notYetIncluded, addTerm, curModel=curModel, matName = "CovMatSec")
      maxAdd = which.max(addRes)
    }
  }
  varnames = c(colNames[MustInclude==1], colnames(temDat)[shouldInclude==1])
  fmla<-as.formula(paste("TreatIndicator~",paste(varnames,sep="",collapse="+")))
  curModel = glm(formula=fmla,family=binomial(logit),data=CovMat)
  
  #return the final propensity score model
  return(curModel)		
  
}



estimateMean = function(y, propScoreUse, x, knots, boundKnots)
{
  valReturn = NULL
  basis = NULL
  lenKnots = length(knots)
  if(lenKnots > 2)
  {
    basis = ns(propScoreUse, knots = knots[2:(lenKnots -1)], Boundary.knots = boundKnots)
    valReturn = lm(y ~ predict(basis, propScoreUse) + x);
  }
  else
  {
    basis = ns(propScoreUse)
    valReturn = lm(y ~ predict(basis, propScoreUse) + x)
  }
  
  return(list(reg=valReturn, basis = basis));
}


estimationPredict = function(coeff, propScore, x, basis)
{
  if(dim(as.matrix(x))[2] > 1)
  {
    Xall = model.matrix(~predict(basis, propScore) + x)
  } else
  {
    Xall = model.matrix(~predict(basis, propScore))
  }
  
  return (Xall %*% as.matrix(coeff));	
}

riwish2 = function(a,v,M) {riwish(v,M)}

betaSamp = function(matNum,sigMat,Bhat,numOut,XtX)
{
  return(rmnorm(1,matrix(Bhat,length(Bhat),1),kronecker(matrix(sigMat,numOut,numOut),XtX)))
}

#a function that samples from the posterior distribution of the outcome model parameters
#Y is the outcome
#X is vector of covariates
#nNumSamp is number of samples (or imputed datasets)
sampleCoeffB = function(Y, X, nNumSamp)
{
  #transforming X to be a matrix
  X = as.matrix(X)
  #identifying how many units are in the dataset
  numObs = length(X[,1])
  
  #number of coefficients per outcome
  numCoef = length(X[1,])
  
  
  
  
  #applying the logistic regression formula to calculate the science or the just Y given X
  outForm = sprintf("Y~.-1")
  lmObj = bayesglm(formula(outForm),data=data.frame(X),family=binomial(),maxit=1000)
  
  #finding which of the parameters is not null
  noNa = which(!is.na(coef(lmObj)))
  
  #obtaining the asymptotic variance-covariance of the model  (coefficients)
  #this is the normal approximation matrix to the model paramters
  sigOpt = vcov(lmObj)
  
  #a vector with the size of the number of parameters filled with 1s
  sigVec = rep(1, nNumSamp)
  
  #checking that the asymptotic variance-covariance is positive definite
  #(it is not when there are some numerical issues)
  if(!is.positive.definite(sigOpt))
  {
    #making sure that the variance-covariance matrix is positive definite
    sigOpt = make.positive.definite(sigOpt)
  }
  
  #number of coeffiecients that can be estimated
  nCoef = length(noNa)
  
  #total number of coefficients
  totalCoef = length(coef(lmObj))
  
  #Generating samples from the posterior distribution of the coefficients (here from the normal approximation)
  coefMat = matrix(rnorm(nNumSamp*nCoef) %*% kronecker(diag(sqrt(sigVec)), chol(sigOpt))  + coef(lmObj)[noNa], nCoef, nNumSamp)
  
  #for parameter that can't be estimated filling the values with 0s
  if(totalCoef - nCoef > 0)
  {
    coefMat = rbind(coefMat, matrix(0, totalCoef - nCoef, nNumSamp))
  }
  
  #returning the samples from the asymptotic distribution of the coefficients
  return(list(sig=NULL, coef=coefMat, numCov = numCoef, numOut = 1))
}

#function to sample from inverse wishart distribution with v degrees of freedom and scale matrix M
riwish2 = function(a,v,M) {riwish(v,M)}

#function for sampling from the posterior distribution of the outcomes using linear regression models when Y is scalar a vector of contionuous variables
#Y is the outcome
#X is vector of covariates
#nNumSamp is number of samples
sampleCoeffC = function(Y, X, nNumSamp)
{
  #transforming X to be a matrix
  X = as.matrix(X)
  #identifying how many units are in the dataset
  numObs = length(X[,1])
  
  #number of coefficients per outcome
  numCoef = length(X[1,])
  
  #number of outcome values
  numOut = length(as.matrix(Y)[1,])
  
  #QR decomposition of X
  qrMat = qr(X)
  
  #identifying the rank of X and making sure the matrix is invertabile
  rankX = qrMat$rank
  if(numObs - numOut - rankX + 1 < numOut)
  {
    rankX = numObs - 2*numOut + 1
  }
  
  #using only the covariates that can be inverted
  pivotX = qrMat$pivot[1:rankX]
  XX = X[,pivotX]
  
  #calculating the linear regression solution of Y given X (X^{t}X)^{-1}X^{t}Y
  XtX = t(XX)%*%XX
  if(!is.positive.definite(XtX))
  {
    XtX = make.positive.definite(XtX)
  }
  XtXinv = solve(XtX)
  
  #obtaining the coefficient of the linear regression
  Beta = XtXinv%*%t(XX)%*%Y
  
  #obtaining the variance covariance matrix
  SE = t(Y-XX%*%Beta)%*%(Y-XX%*%Beta) 
  
  #sampling sigma from inverse wishart posterior
  sigMat = sapply(1:nNumSamp, riwish2, v=numObs - numOut - rankX + 1, M=SE) 
  
  #sampling the coefficient from multivariate normal distribution
  coefMat = sapply(1:nNumSamp, betaSamp, Bhat = Beta, sigMat = sigMat, numOut = numOut, XtX = XtXinv)
  
  #filling in zero for the coefficients that can't be estimated
  if(rankX < numCoef)
  {
    tempMat = matrix(0,numCoef*numOut,nNumSamp)
    updateNum = rep(pivotX,numOut)+rep(seq(0,numCoef*numOut - 1,numCoef),each=length(pivotX))
    tempMat[updateNum,] = coefMat
    coefMat = tempMat
  }
  
  return(list(sig=sigMat, coef=coefMat, numCov = numCoef, numOut = numOut))
}

estimandProp = function(Ytrt, Ycon, obsUse)
{
  
  temp.obsUse  = obsUse
  tempVecIm <-mean(as.matrix(Ytrt)[temp.obsUse,1]) -  mean(as.matrix(Ycon)[temp.obsUse,1])
  Q = tempVecIm
  numA = length(which(as.matrix(Ytrt)[temp.obsUse,1] == 1 & as.matrix(Ycon)[temp.obsUse,1] == 1))
  numB = length(which(as.matrix(Ytrt)[temp.obsUse,1] == 1 & as.matrix(Ycon)[temp.obsUse,1] == 0))
  numC = length(which(as.matrix(Ytrt)[temp.obsUse,1] == 0 & as.matrix(Ycon)[temp.obsUse,1] == 1))
  numD = length(which(as.matrix(Ytrt)[temp.obsUse,1] == 0 & as.matrix(Ycon)[temp.obsUse,1] == 0))
  U = sqrt((numA+numD)*(numB+numC) + 4*numB*numC) /(length(temp.obsUse)*sqrt(length(temp.obsUse))) 
  
  return(list(Q=QQ,U=matrix(UU)))  
}

estimandMean = function(Ytrt, Ycon, obsUse)
{
  temp.obsUse <-obsUse
  tempVecIm <-as.matrix(Ytrt)[temp.obsUse,] -  as.matrix(Ycon)[temp.obsUse,]
  Q <- apply(as.matrix(tempVecIm),2,mean)
  U <- (var(tempVecIm) / length(temp.obsUse))
  
  QQ = Q
  UU = U
  
  return(list(Q=QQ,U=matrix(UU)))  
}


#function that recieve coefficients (coeff), matrix of covariates (xx), propensity score (propScore) and cubic spline basis (basis). 
#The function generates the predicted linear value of for every row of the covariates
estimationPredict = function(coeff, propScore, x, basis)
{
  if(dim(as.matrix(x))[2] > 1)
  {
    Xall = model.matrix(~as.matrix(predict(basis, propScore) + x))
  } else
  {
    Xall = model.matrix(~predict(basis, propScore))
  }
  
  return (Xall %*% as.matrix(coeff));	
}

# T statistic calculation.
t.prob <- function(t, n){2*pt(-abs(t),df=n)}