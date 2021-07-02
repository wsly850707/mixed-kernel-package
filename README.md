# mixed-kernel-package
# Installation
     library(devtools)  
     install_github("wsly850707/mixed-kernel-package")

# Requirements
mk requires the following packages- hdi, MASS, mvtnorm

# Introduction
mk is a function that combines both the linear and Gaussian kernel with different weights to capture linear or nonlinear interaction of two genes. In a high-dimensional setup, a Cauchy transformation of the p-values obtained using the HDI package under different weights is implemented to obtain an aggregate p-values for testing each interaction effect. 

# Functions
Three functions are loaded when the mk package is installed:
1. CauchyP--Combine the p-values obtained from the HDI package under different weights.
2. mta--Multiple testing adjustment step when fitting multiple variables in the interaction model. An effective number of tests can be estimated. The adjusted p-values are obtained by multiplying this number by the p-values obtained from HDI. 
3. delp--Organize the p-values obtained from mta. If w p-values are all equal to 1 for a variable, they are not combined by the CauchyP step. The final p-value is set as 1. If w p-values are all less than 1, they are combined by CauchyP. Those p-values larger than 1 do not contriubte to the final aggregrated p-values.

# Datasets
Five datasets are loaded when the mk package is installed. They are the datasets of the lung cancer example analyzed in case study 2 in the manuscript with ending=lung cancer status, expression =gene expression, pro=probe name, genename=gene name and pheno=phenotype.

# Case study Example
We show the implementation of the case study 2 shown in our paper. This code is also used as an Example in our package. We just show how to obtain the p-values after multiple test adjustments and the Cauphy combination in pathway hsa05223. The heatmap code can be obtained from the Example of mk package. See code below.

    library(mk)
    library(hdi)
    library(mvtnorm)
    library(MASS)
    
    #deleting the suspicious cancer samples
    ending1<-as.matrix(ending[-c(188:192),])
    #According to the status, we replace the status with 0 and 1. 
    ending1[,2][which(ending1[,2]=="no cancer")]<-0
    ending1[,2][which(ending1[,2]=="cancer")]<-1
    #When loading ending, the rowname of it is regarded as first column, so we adjust it.
    ending2<-as.matrix(ending1[,-1])
    rownames(ending2)<-ending1[,1]

    #The gene expression data are obtained based on probeset, so we need to map probesets to genes.
    proexp<-cbind(pro[,2],expression[,-1])
    #For genes with different probesets , we calculate the average of different probesets as the gene level expression.
    proexp1<-aggregate(x=proexp[,-1],by=list(proexp[,1]),FUN=mean)

    pheno1<-as.matrix(pheno[,-1])
    rownames(pheno1)<-as.matrix(pheno[,1])

    ##Applying a logistic regression to select signficant clinical covariates to include in the model. Samples with missing values were deleted.
    phenona<-is.na(pheno1[,1])
    indexna<-which(phenona==TRUE)
    pheno2<-pheno1[-indexna,]
    state1<-as.numeric(ending2[,1])[-indexna]
    age1<-as.numeric(pheno2[,1])
    packyear1<-as.numeric(pheno2[,5])
    lg<-glm(state1~age1+pheno2[,2]+pheno2[,3]+pheno2[,4]+packyear1+pheno2[,6]+pheno2[,7],"binomial")
    summary(lg)
    #Only age and lymphadenopathy are remained as meaningful covariates. 
    lymphadenopathy<-pheno2[,7]
    xx<-data.frame(age1,lymphadenopathy)
    xxx<-model.matrix(~lymphadenopathy,xx)
    cov<-cbind(age1,xxx[,2])

    ##mapping the KEGG pathway genes and GEO genes
    genename<-as.matrix(genename)
    int<-intersect(proexp1[,1],genename[,2])
    rownames(proexp1)<-proexp1[,1]
    proexp2<-proexp1[,-1]
    #There are 66 genes mapped to the pathway.
    expressiongene<-proexp2[int,]
    #The rows represent genes which need to be transposed. 
    expressiongene1<-t(expressiongene[,-c(188:192)])
    #samples with missing values were deleted
    expressiongene2<-expressiongene1[-indexna,]

    w=seq(0,1,by=0.2)
    nse=2000 # set the see number

    ###The speed of the mk package is limited by the HDI inference step which fits a high-dimensional model when testing each variable. This is the step that takes the most time ### 
    hsa05223<-mk(expressiongene2,state1,cov,w,nse,family="binomial",method="lasso",standardize=TRUE, adjust=TRUE)
    #hsa05223 is a p-value vector which is obtained after the Cauchy p-value combination and multiple test adjustment.

# Simulation Example

First we show the simulation example with the d=2 model and the continuous response. We set coefficient terms=0.2 to assess the testing power of the overall test and interaction terms.  

As in our paper,we set the simulation replicates to 1000. Users need to install the mk package from Github before use it.

    library(devtools)
    install_github("wsly850707/mixed-kernel-package")

    #generating the dataset which contains two main effects, interaction effect and continuous response
    library(mk)
    library(mvtnorm)
    library(MASS)

    n=100 # sample size
    d=2 # number of gene variables
    B=1000 # number of replications
    C=seq(0,1,by=0.2) # grid of the weight to generate the data. C=0, 1 correspond to nonlinear and linear interaction effect, respectively. Any values b/w 0 and 1 give a mixed effect.

    #test_P is the matrix of p-value obtained from B times of GLM. 
    test_Pl=test_Pg=test_Pm=matrix(0,B,3)
    #overalltest_P is the matrix of p-value obtained from B times of overall test.
    overalltest_Pl=overalltest_Pg=overalltest_Pm=matrix(0,B,1)
    #powerl, powerg and powerm are the power and size matrix.
    power1=power2=power3=powerall=matrix(0,length(C),3)

    for (t in 1:length(C))
    {
     c=C[t]
      for (i in 1:B)
       {
        # assuming independence between gene variables
        sigma=diag(1,d)
       
        #The main effect variables are generated from a multivariate normal distribution with mean 0 and covariance sigma.    
        x=mvrnorm(n,rep(0,d),sigma)
        x1=x[,1]
        x2=x[,2]    
        error=rnorm(n,0,1)# Error terms are generated from a standand normal distribution.
        
        #Using linear kernel, gaussian kernel and a kernel that is similar to guassian kernel to define the interaction effect, we combine three types of interaction term with main effects matrix to constuct three types of independent variable matrix.
        z1=scale(x1*x2) #linear kernel
        kernell<-cbind(x1,x2,z1) 
        z2=scale(exp(-0.5*((x1-x2)^2))) #Gaussian kernel
        kernelg<-cbind(x1,x2,z2)
        z.2=scale(exp(-0.5*(abs(x1-x2)))) ##The reason we construc this kernel is that we will use this kernel as a nonlinear kernel to generate response variable to avoid using the same kernel to generate and analyze the data.
        Y=0.5+0.2*x1+0.2*x2+c*0.2*z1+(1-c)*0.2*z.2+error      
   
       #GLM in linear interaction model   
        fitlinear <- glm (Y~kernell)
        test_Pl[i,]<-summary(fitlinear)$coefficients[2:4,4]
        anoval <- anova(fitlinear,test="Chisq")
        overalltest_Pl[i,]=anoval[2,5]    

      #GLM in nonlinear interaction model
        fitgaussian <- glm (Y~kernelg)
        test_Pg[i,]<-summary(fitgaussian)$coefficients[2:4,4]
        anovag <- anova(fitgaussian,test="Chisq")
        overalltest_Pg[i,]=anovag[2,5]

      #mk in mixed interaction model
        w=seq(0,1,by=0.2)
        pval=matrix(0,length(w),ncol(kernell))
        pvaloverall=matrix(0,length(w),1)
          for (o in 1:length(w)){
           #Establishing the mixed interaction vector according to each w, we combine it  with x to constitute the predictor variables to fit glm.
           z3=w[o]*z1+(1-w[o])*z2
           kernelmix<-cbind(x1,x2,z3)
           fitmix <- glm (Y~kernelmix)
           pval[o,]=summary(fitmix)$coefficients[2:4,4]
           anovam <- anova(fitmix,test="Chisq")
           pvaloverall[o,]=anovam[2,5]
         }
        test_Pm[i,]<-mk(x,Y,cov=NULL,w,nse,family="gaussian",method="lasso",standardize=TRUE,adjust=FALSE) 
        overalltest_Pm[i,]<-as.matrix(apply(pvaloverall,2,CauchyP))
           
       }
  
    #calculating the power or power of x1
    power1[t,]<-c(length(which(test_Pl[,1]<0.05)),length(which(test_Pg[,1]<0.05)),length(which(test_Pm[,1]<0.05)))/B
    #calculating the power or power of x2
    power2[t,]<-c(length(which(test_Pl[,2]<0.05)),length(which(test_Pg[,2]<0.05)),length(which(test_Pm[,2]<0.05)))/B
    #calculating the power or power of kernel
    power3[t,]<-c(length(which(test_Pl[,3]<0.05)),length(which(test_Pg[,3]<0.05)),length(which(test_Pm[,3]<0.05)))/B
    #calculating the power or power of model
    powerall[t,]=c(length(which(overalltest_Pl[,1]<0.05)),length(which(overalltest_Pg[,1]<0.05)),length(which(overalltest_Pm[,1]<0.05)))/B
    colnames(power1)<-colnames(power2)<-colnames(power3)<-colnames(powerall)<-c("L","NL","M")  
    rownames(power1)<-rownames(power2)<-rownames(power3)<-rownames(powerall)<-c("c=0","c=0.2","c=0.4","c=0.6","c=0.8","c=1")  
    }

    library(round)
    round(powerall,3) #power for the overall test under different data generating mechanisms and under model L, NL and M 
    round(power1,3) #power for testing x1 main effect under different data generating mechanisms and under model L, NL and M 
    round(power2,3) #power for testing x2 main effect under different data generating mechanisms and under model L, NL and M
    round(power3,3) #power for testing interaction x1x2 effect under different data generating mechanisms and under model L, NL and M

The simulation example with the d>2 scenario and the continuous response. 

    library(mk)
    library(mvtnorm)
    library(MASS)

    n=100 # sample size
    d=10 # number of gene variables
    B=c(1:10) # number of replications, we set it to 1000 in our paper
    C=seq(0,1,by=0.2) # grid of the weight to generate the data. C=0, 1 correspond to nonlinear and linear interaction effect, respectively. Any values b/w 0 and 1 give a mixed effect.

    #defining the column number of interaction effect matrix according to d
    inter_num<-choose(d,2)

    #Test_P is the matrix of p-values obtained from B replications. 
    #test_Pl for linear, test_Pg for Guassian, test_Pm for the mixed kernel
    test_Pl=test_Pg=test_Pm=matrix(0,length(B),(d+inter_num))
    #powerl, powerg and powerm are the power and size matrix for applying the linear, Gaussian and mixed kernel function.
    powerl=powerg=powerm=matrix(0,length(C),(d+inter_num))

    #for each c value, running length(B) times replicates
    #C=c(0, 0.2, 0.4, 0.6, 0.8, 1), when C=0, the interaction effect is purely nonlinear; when C=1, the interacton effect is purely linear; when C takes values between 0 and 1, the effect is a mixture of linear and nonliear. ###

    c=0.2 # can vary c to get different results under different interaction effect conditions 
    for (i in 1:length(B)) 
     {
     # assuming independence between gene variables
     sigma=diag(1,d)
    
     #The main effect variables are generated from a multivariate normal distribution with mean 0 and covariance sigma.
     set.seed(B[i]+1000)
     x=mvrnorm(n,rep(0,d),sigma)
     colnames(x)<-paste0("x",1:ncol(x))
    
     error<-rnorm(n,0,0.8) # Error terms are generated from a normal distribution with mean 0 and sd 0.8.

     #defining linear interaction matrix
     inter_setl<-matrix(0,n,)
     for (m in 1:(ncol(x)-1)){
      for (e in (m+1):(ncol(x))){
        m1<-x[,m]
        m2<-x[,e]
        templ<-cbind(m1,m2)
        #using linear kernel to generate the linear interaction
        templ1<-apply(templ,1,function(x) x[1]*x[2])
        templ1<-as.matrix(templ1)
        #naming the interaction terms as xd-xe
        colnames(templ1)<-paste0(colnames(x)[m],"-",colnames(x)[e])
        #putting each interaction term into the temporaty matrix 
        inter_setl<-cbind(inter_setl,templ1)
       } 
     }
    inter_setl1<-scale(inter_setl[,-1])
    kernell<-scale(cbind(x,inter_setl1))
    
    #defining the nonlinear interaction matrix
    inter_setg<-matrix(0,n,)
    for (f in 1:(ncol(x)-1)){
      for (g in (f+1):(ncol(x))){
        f1<-x[,f]
        f2<-x[,g]
        tempg<-cbind(f1,f2)
        #using RBF kernel to generate the nonlinear interaction
        tempg1<-apply(tempg,1,function(x) exp(-0.5*((x[1]-x[2])^2)))
        tempg1<-as.matrix(tempg1)
        #naming the interaction terms as xf-xg
        colnames(tempg1)<-paste0(colnames(x)[f],"-",colnames(x)[g])
        #putting each interaction term into the temporaty matrix 
        inter_setg<-cbind(inter_setg,tempg1)
      } 
    }
    inter_setg1<-scale(inter_setg[,-1])
    kernelg<-scale(cbind(x,inter_setg1))

    ###Generating the data. ### 
    ###For the nonlinear interaction effect, we used a different kernel function than the RBF one to avoid using the same kernel to generate and analyze the same dataset ###
    inter_setg.<-matrix(0,n,)
    for (j in 1:(ncol(x)-1)){
      for (l in (j+1):(ncol(x))){
        j1<-x[,j]
        j2<-x[,l]
        tempg.<-cbind(j1,j2)
        #using this kernel to generate the nonlinear interactions
        tempg.1<-apply(tempg.,1,function(x) exp(-0.5*(abs(x[1]-x[2]))))
        tempg.1<-as.matrix(tempg.1)
        #naming the interaction terms as xj-xl
        colnames(tempg.1)<-paste0(colnames(x)[j],"-",colnames(x)[l])
        #putting each interaction term into the temporaty matrix 
        inter_setg.<-cbind(inter_setg.,tempg.1)
      } 
    }
    inter_setg.1<-scale(inter_setg.[,-1])
    kernelg.<-scale(cbind(x,inter_setg.1))
    
    alpha<-rep(0,dim(x)[2])
    alpha[1:5]<-rep(0.2,5)
    beta<-rep(0,dim(inter_setg1)[2])
    beta[c(2,11,16,19,27)]<-rep(0.2,5)
    Y<-0.5+x%*%alpha+c*(inter_setl1%*%beta)+(1-c)*(inter_setg.1%*%beta)+error
    
    # Applying the HDI package assuming a linear interaction model
    library(hdi)
    set.seed(1000)
    fitlinear<-lasso.proj(kernell,Y)
    test_Pl[i,]<-fitlinear$pval
      
    # Applying the HDI package assuming a nonlinear interaction model
    set.seed(1000)
    fitgaussian<-lasso.proj(kernelg,Y)
    test_Pg[i,]<-fitgaussian$pval
     
    # Applying the mk package assuming the mixed interaction model
    w=seq(0,1,by=0.2)
    nse=1000 # set the seed number
    test_Pm[i,]<-mk(x,Y,cov=NULL,w,nse,family="gaussian",method="lasso",standardize=TRUE,adjust=FALSE) 
    cat("Simu Rep = ", i, "\n")
    }
    
    ########Power and size calculation########
    #calculating the power and size of the linear interaction model.
    test_Pl1<-test_Pl
    test_Pl1[test_Pl1>=0.05]=2
    test_Pl1[test_Pl1<0.05]=1
    test_Pl1[test_Pl1==2]=0
    powerl<-t(as.matrix(apply(test_Pl1,2,function(x) sum(x)/length(x))))
    #calculating the power and size of the nonlinear interaction model
    test_Pg1<-test_Pg
    test_Pg1[test_Pg1>=0.05]=2
    test_Pg1[test_Pg1<0.05]=1
    test_Pg1[test_Pg1==2]=0
    powerg<-t(as.matrix(apply(test_Pg1,2,function(x) sum(x)/length(x))))
    #calculating the power and size of the mixed interaction model
    test_Pm1<-test_Pm
    test_Pm1[test_Pm1>=0.05]=2
    test_Pm1[test_Pm1<0.05]=1
    test_Pm1[test_Pm1==2]=0
    powerm<-t(as.matrix(apply(test_Pm1,2,function(x) sum(x)/length(x))))
    colnames(powerl)<-colnames(powerg)<-colnames(powerm)<-colnames(kernell)
    rownames(powerl)<-rownames(powerg)<-rownames(powerm)<-c("c=0.2")

    #outputting the power matrix of the linear interaction model
    mainpowerl<-t(as.matrix(powerl[,which(alpha!=0)])) #power matrix of main effects
    interpowerl<-t(as.matrix(powerl[,length(alpha)+which(beta!=0)])) #power matrix of interaction effects
    #outputting the size of main effects and interaction effects of the linear interaction model
    mainsizel<-t(as.matrix(powerl[,which(alpha==0)])) #size matrix of main effects
    intersizel<-t(as.matrix(powerl[,length(alpha)+which(beta==0)])) #size matrix of interaction effects
    #calculating the average size of main effect and interaction effect
    mainavsizel<-as.matrix(apply(mainsizel,1,mean)) #average size matrix of main effects
    colnames(mainavsizel)<-"mainsize" 
    interavsizel<-as.matrix(apply(intersizel,1,mean)) #average size matrix of interaction effects
    colnames(interavsizel)<-"intersize"
    #combining power matrix and average size matrix 
    poweravsizeL<-cbind(mainpowerl,mainavsizel,interpowerl,interavsizel)

    #outputting the power matrix of the nonlinear interaction model
    mainpowerg<-t(as.matrix(powerg[,which(alpha!=0)])) #power matrix of main effects
    interpowerg<-t(as.matrix(powerg[,length(alpha)+which(beta!=0)])) #power matrix of interaction effects
    #calcuating the size of main effects and interaction effects of the nonlinear interaction model
    mainsizeg<-t(as.matrix(powerg[,which(alpha==0)])) #size matrix of main effects
    intersizeg<-t(as.matrix(powerg[,length(alpha)+which(beta==0)])) #size matrix of interaction effects
    #calculating the average size of main effects and interaction effects
    mainavsizeg<-as.matrix(apply(mainsizeg,1,mean)) #average size matrix of main effects
    colnames(mainavsizeg)<-"mainsize"
    interavsizeg<-as.matrix(apply(intersizeg,1,mean)) #average size matrix of interaction effects
    colnames(interavsizeg)<-"intersize"
    #combining power matrix and average size matrix 
    poweravsizeNL<-cbind(mainpowerg,mainavsizeg,interpowerg,interavsizeg)

    #outputting the power matrix of the mixed interaction model
    mainpowerm<-t(as.matrix(powerm[,which(alpha!=0)])) #power matrix of main effects
    interpowerm<-t(as.matrix(powerm[,length(alpha)+which(beta!=0)])) #power matrix of interaction effects
    #calculating the size of main effects and interaction effects of the mixed interaction model
    mainsizem<-t(as.matrix(powerm[,which(alpha==0)])) #size matrix of main effects
    intersizem<-t(as.matrix(powerm[,length(alpha)+which(beta==0)])) #size matrix of interaction effects
    #calculating the average size of main effects and interaction effects
    mainavsizem<-as.matrix(apply(mainsizem,1,mean)) #average size matrix of main effects
    colnames(mainavsizem)<-"mainsize"
    interavsizem<-as.matrix(apply(intersizem,1,mean)) #average size matrix of interaction effects
    colnames(interavsizem)<-"intersize"
    #combining power matrix and average size matrix 
    poweravsizeM<-cbind(mainpowerm,mainavsizem,interpowerm,interavsizem)
    rownames(poweravsizeL)<-rownames(poweravsizeNL)<-rownames(poweravsizeM)<-c("c=0.2")

    library(round)
    round(poweravsizeL,3) #power and average size for main and interaction effects under different data generating mechanisms and under model L 
    round(poweravsizeNL,3) #power and average size for main and interaction effects under different data generating mechanisms and under model NL 
    round(poweravsizeM,3) #power and average size for main and interaction effects under different data generating mechanisms and under model M 



