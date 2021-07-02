#The mk function is used to test significant main and interaction effects in high-dimensional (generalized) linear model.
#R version: R/3.6.3.

#Description
#This function can be used to perform inference to test significant main and interaction effects in high-dimensional (generalized) linear model. HDI package is needed to be installed before running this function. 

#library(hdi)

#Usage
mk<- function (X, Y, cov, w, nse, family = c("gaussian","binomial"), method = c("lasso", "ridge"),standardize=c("TRUE","FALSE"),adjust =c("TRUE","FALSE"), ...){	
  
  #Arguments
  #X   main effect matrix that can be gene expression data.
  #Y   response vector that can be either continuous or binary.
  #cov   data.frame or matrix of covariates. 
  #w   vector of weights used to analyze the data with a mixed kernel. p-values are combined with the Cauchy comination method under different w values. 
  #nse   number of seed that been set before HDI.
  #family   either "gaussian" or "binomial", relying on the type of response. 
  #method   either "lasso" or "ridge" to estimate the main and interaction effects.
  #standardize   Should design matrix be standardized to unit column standard deviation.
  #adjust   Should the p-values obtained from HDI be adjusted via multiple test adjustment.
  # ... : other arguments passed to hdi.
  
  #Value 
  #Caup   Vetcor of p-values that are cauchy combination of p-values obtained under different weights.
  
  ########Step-1  define linear and nonlinear interaction#########
  #Choosing the ith and jth X variables, we use a linear kernel to define the linear interaction b/w them. 
  temp_totall<-matrix(0,nrow(X),)
  for (i in 1:(ncol(X)-1)){
    for (j in (i+1):(ncol(X))){
      a1<- X[,i]
      a2<- X[,j]
      templ<- cbind(a1,a2)
      #define linear interaction     
      templ1<- apply(templ,1,function(x) x[1]*x[2])
      templ1<- as.matrix(templ1)
      #name the interaction terms as Xi_Xj
      colnames(templ1)<- paste0(colnames(X)[i],"_",colnames(X)[j])
      temp_totall<- cbind(temp_totall,templ1)
    } 
  }
  #scale the variables
  kernell<- scale(temp_totall[,-1])
  
  #Choosing the ith and jth X variables, we use a Gaussian kernel to define the non-linear interaction b/w them.
  temp_totalnl<- matrix(0,nrow(X),)
  for (i in 1:(ncol(X)-1)){
    for (j in (i+1):(ncol(X))){
      a1<- X[,i]
      a2<- X[,j]
      tempnl<- cbind(a1,a2)
      #use a Gaussian kernel to define the non-linear interaction
      tempnl1<- apply(tempnl,1,function(x) exp(-0.5*((x[1]-x[2])^2)))
      tempnl1<- as.matrix(tempnl1)
      #name the interaction terms as Xi_Xj
      colnames(tempnl1)<- paste0(colnames(X)[i],"_",colnames(X)[j])
      temp_totalnl<- cbind(temp_totalnl,tempnl1)
    } 
  }
  kernelnl<- scale(temp_totalnl[,-1])
  
  if (adjust){
    ########Step-2   High-dimensional Inference (HDI)#########
    #Based on the linear and nonlinear interaction matrices, we constitute the mixed interaction matrix under different w.
    if (is.null(cov)) {
      #construct a p-value matrix containing p-values of X and its interaction terms      
      pval=matrix(0, length(w), ncol(X)+ncol(kernell)) 
      pvaladj=matrix(0, length(w), ncol(X)+ncol(kernell)) 
      for (k in 1: length(w)){
        #establishing the mixed interaction matrix according to each w.
        kernelm = w[k]* kernell+ (1-w[k])* kernelnl
        Eff <- cbind (scale(X), kernelm)
        if (ncol(X)==2){
          #we use glm to test significant main and interaction effects when d=2
          #glm
          fit<-glm(Y~Eff,family = family)
          pval[k,]<-summary(fit)$coefficients[2:(ncol(Eff)+1),4]
          pvaladj[k,]<-pval[k,]*mta(Eff)
        } else #we use hdi to test significant main and interaction effects
        {
          #HDI
          set.seed(nse)
          if (method == "lasso") fit <- lasso.proj (Eff, Y, family = family,standardize=standardize) else fit <- ridge.proj(Eff, Y, family = family,standardize=standardize)
          pval[k,]<- fit$pval
          pvaladj[k,]<-pval[k,]*mta(Eff)
          cat("w =",(k-1)/(length(w)-1),"\n")
        }
      }
      #Name the column of p-value matrix obtained from HDI
      colnames(pvaladj)<-c(colnames(X),colnames(kernell))
    } else {
      #If there are covariates in the model, then main effects, interaction effects and covariates are all fitted.
      pval=matrix(0, length(w), ncol(cov)+ncol(X)+ncol(kernell))
      pvaladj=matrix(0, length(w), ncol(cov)+ncol(X)+ncol(kernell))
      for (k in 1 :length(w)){ 
        kernelm = w[k]* kernell+ (1-w[k])* kernelnl 
        Eff<- cbind (scale(cov), scale(X), kernelm)
        Eff1<-cbind(scale(X), kernelm)
        if (ncol(X)==2){
          #glm
          fit<-glm(Y~Eff,family = family)
          pval[k,]<-summary(fit)$coefficients[2:(ncol(Eff)+1),4]
          pvaladj[k,]<-pval[k,]*mta(Eff)
        } else {
          #HDI
          set.seed(nse)
          if (method == "lasso") fit <- lasso.proj (Eff, Y, family = family,standardize=standardize) else fit <- ridge.proj(Eff, Y, family = family,standardize=standardize)
          pval[k,]<- fit$pval
          pvaladj[k,]<-pval[k,]*mta(Eff)
          cat("w =",(k-1)/(length(w)-1),"\n")
        }
      }
      #Name the column of p-value matrix obtained from HDI
      colnames(pvaladj)<-c(colnames(cov),colnames(X),colnames(kernell))
    }
    #There are part of p-values larger than 1 after multiple testing adjustment, we let them be 1. 
    pvaladj[which(pvaladj>=1)]=1
    
    
    ########STEP-3   P-value combination#########
    #P-value obtained from HDI according to different w need to be combined for further inference.
    Caup<- as.matrix(apply(pvaladj,2,delp))
    return(Caup)
  } else {
    if (is.null(cov)) {
      pval=matrix(0, length(w), ncol(X)+ncol(kernell)) 
      for (k in 1: length(w)){
        kernelm = w[k]* kernell+ (1-w[k])* kernelnl
        Eff <- cbind (scale(X), kernelm)
        if (ncol(X)==2){
          #glm
          fit<-glm(Y~Eff,family = family)
          pval[k,]<-summary(fit)$coefficients[2:(ncol(Eff)+1),4]
        } else {
          #HDI
          set.seed(nse)
          if (method == "lasso") fit <- lasso.proj (Eff, Y, family = family,standardize=standardize) else fit <- ridge.proj(Eff, Y, family = family,standardize=standardize)
          pval[k,]<- fit$pval
          cat("w =",(k-1)/(length(w)-1),"\n")
        }
      }
      colnames(pval)<-c(colnames(X),colnames(kernell))
    } else {
      pval=matrix(0, length(w), ncol(cov)+ncol(X)+ncol(kernell))
      for (k in 1 :length(w)){ 
        kernelm = w[k]* kernell+ (1-w[k])* kernelnl 
        Eff<- cbind (scale(cov), scale(X), kernelm)
        if (ncol(X)==2){
          #glm
          fit<-glm(Y~Eff,family = family)
          pval[k,]<-summary(fit)$coefficients[2:(ncol(Eff)+1),4]
        } else {
          set.seed(nse)
          if (method == "lasso") fit <- lasso.proj(Eff, Y, family = family,standardize=standardize) else fit <- ridge.proj(Eff, Y, family = family,standardize=standardize)
          pval[k,]<- fit$pval
          cat("w =",(k-1)/(length(w)-1),"\n")
        }
      }
      colnames(pval)<-c(colnames(cov),colnames(X),colnames(kernell))
    }
    
    ########STEP-3   P-value combination#########
    #P-value obtained from HDI according to different w need to be combined for further inference.
    Caup<- as.matrix(apply(pval,2,CauchyP))
    return(Caup)
  } 
}


#Multiple testing adjustment

mta<-function(x)
{
  corx<-cor(x)
  ecl<-eigen(corx)
  ecl1<-ecl$values
  ecl11<-abs(ecl1)
  ecl111<-ifelse(ecl11>=1,1,0)
  ecl112<-ecl11-floor(ecl11)
  fx<-ecl111+ecl112
  sum.fx<-sum(fx)
  return(sum.fx)
}

#Cauchy combination
CauchyP = function(p)
{
  cct = sum(tan(pi*(0.5-p))/length(p))
  return(1-pcauchy(cct))
}
#If w p-values are all equal to 1 in a variable, then they are not combined by Cauchy. Their combination is 1. 
#If w p-values are all less than 1, then they are combined by Cauchy. If part of p-values are less then 1 and the others are equal to 1, then we just use Cauchy to combine the ones less than 1.
delp<-function(x){
  if (length(x)==1){
    out<-x
  } else {
    index1<-which(x<1)
    index2<-which(x==1)
    if(length(index1)==length(x)){
      out<-CauchyP(x)
    } else if (length(index2)==length(x)){
      out<-CauchyP(x)
    } else {
      x1<-c(x[-index2])
      out<-CauchyP(x1)
    }
  }
  return(out)
} 

