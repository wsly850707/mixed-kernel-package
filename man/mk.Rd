\name{mk}
\alias{mk}

\title{
Mixed Kernel Model
}
\description{
This function can be used to perform inference to test significant main and interaction effects in high-dimensional (generalized) linear model. HDI package is needed to be installed before running this function.
}
\usage{
mk(X, Y, cov, w, nse, family = c("gaussian", "binomial"), method = c("lasso", "ridge"), 
standardize = c("TRUE", "FALSE"), adjust = c("TRUE", "FALSE"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
main effect matrix that can be gene expression data.
}
  \item{Y}{
response vector that can be either continuous or binary.
}
  \item{cov}{
 data.frame or matrix of covariates.
}
  \item{w}{
vector of weights used to analyze the data with a mixed kernel. p-values are combined with the Cauchy comination method under different w values.
}
  \item{nse}{
number of seed that need to be set before HDI.
}
  \item{family}{
either 'gaussian' or 'binomial', relying on the type of response.
}
  \item{method}{
either "lasso" or "ridge" to estimate the main and interaction effects.
}
  \item{standardize}{
Should design matrix be standardized to unit column standard deviation.
}
  \item{adjust}{
Should the p-values obtained from HDI be adjusted via multiple test adjustment.
}
  \item{\dots}{
other arguments passed to hdi.
}
}
\details{
Choosing the ith and jth X variables, we use a linear kernel to define the linear interaction  and gaussian kernel to define the nonlinear interaction between them.

Based on the linear and nonlinear interaction matrices, we constitute the mixed interaction matrix under different w.

P-value obtained from HDI according to different w are combined by Cauchy combination.

If need, the multiple testing adjustment is needed to be done to get the adjusted-values.
}
\value{
\item{Caup}{
Vetcor of p-values that are cauchy combination of p-values obtained under different weights.
}
}
\references{
Zhang CH, Zhang SS.(2014) Confidence intervals for low dimensional parameters in high dimensional linear models. J. R. Stat. Soc. Ser. B Stat. Methodol. 2014; 76:217-242.

Dezeure R, B P, Meier L, et al. High-dimensional inference: Confidence intervals, P-values and R-software hdi. Stat. Sci. 2015; 30: 533-558.

Van De Geer S, B P, Ritov Y, et al. On asymptotically optimal confidence regions and tests for high-dimensional models. Ann. Stat. 2014; 42:1166-1202.
}
\author{
Yan Liu
}

\examples{
#The gene expression data of lung cancer was downloaded from the [GEO] repository 
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE4115 
#corresponding to the case study 2 in the paper.
#Cancer status, clinical characters and gene expression are included in the data. 
#We separated the dataset into three files named ending, charac and expression respectively.
#deleting the suspicious cancer samples
ending1<-as.matrix(ending[-c(188:192),])
#According to the status, we replace the status with 0 and 1. 
ending1[,2][which(ending1[,2]=="no cancer")]<-0
ending1[,2][which(ending1[,2]=="cancer")]<-1
#When loading ending, the rowname of it is regarded as first column, so we adjust it.
ending2<-as.matrix(ending1[,-1])
rownames(ending2)<-ending1[,1]
#The gene expression data are obtained based on probeset, so we need to map probesets to genes.
#The "gene and probeID" file can also be get from the GEO website. 
#Only gene and probeID are picked from the file.
proexp<-cbind(pro[,2],expression[,-1])
#For genes with different probesets , we calculated the average of 
#different probesets as the gene level expression.
proexp1<-aggregate(x=proexp[,-1],by=list(proexp[,1]),FUN=mean)
#input clinical characters data. Seven clinical characters containing age, gender, race, 
#smoking status, pack years of smoking, 
#hemopytsis and lymphadenopathy were picked up into analysis.
pheno1<-as.matrix(pheno[,-1])
rownames(pheno1)<-as.matrix(pheno[,1])

## Apply a logistic regression to select signficant clinical covariates to 
#include in the model. Samples with missing values were deleted.
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

#####mapping the KEGG pathway genes and GEO genes#####
#mapping genes to the KEGG pathway
genename<-as.matrix(genename)
int<-intersect(proexp1[,1],genename[,2])
rownames(proexp1)<-proexp1[,1]
proexp2<-proexp1[,-1]
#There are 64 genes mapped to the pathway.
expressiongene<-proexp2[int,]
#The rows represent genes which need to transposed. 
expressiongene1<-t(expressiongene[,-c(188:192)])
#samples with missing values were deleted
expressiongene2<-expressiongene1[-indexna,]

w=seq(0,1,by=0.2)
nse=2000

###The speed of the mk package is limited by the HDI inference step which fits a 
#high-dimensional model when testing each variable. This is the step that takes the most time ### 
hsa05223<-mk(expressiongene2,state1,cov,w,nse,family="binomial",method="lasso",
standardize=TRUE,adjust=TRUE)
#hsa05223 is a p-value vector which is obtained after the Cauchy p-value 
#combination and multiple test adjustment.

#drawing a heatmap of -log10 p-values of main and interaction effects of genes
#deleting the covariates and do the -log10 transformation
pval<- -log10(as.matrix(hsa05223[-c(1,2),]))
#setting up the p-value matrix for heatmap
#letting the genes array like which are in the heatmap
s<-seq(ncol(expressiongene2)-1,1,by=-1)
lis<-list()
for(o in 1:length(s)){
lis[[o]]<-c(rep(o,s[o])) 
}
unlis<-unlist(lis)
P_value<-matrix(0,ncol(expressiongene2),ncol(expressiongene2))
P_value1<-matrix(0,ncol(expressiongene2),ncol(expressiongene2))

for (r in 1:length(lis)){
  clue1<-which(unlis==r)
#putting the main effects on the diagonal of the heatmap matrix
  P_value[r,r]<-pval[,1][c(1:(ncol(expressiongene2)-1))][r]
  P_value[ncol(expressiongene2),ncol(expressiongene2)]<-pval[,1][ncol(expressiongene2)]
#putting the interaction effects on the other places of heatmap matrix
  P_value[r,c((r+1):ncol(expressiongene2))]<-pval[,1][-c(1:ncol(expressiongene2))][clue1]
  P_value1[r,c((r+1):ncol(expressiongene2))]<-pval[,1][-c(1:ncol(expressiongene2))][clue1]
}

  rownames(P_value)<-colnames(expressiongene2)
  colnames(P_value)<-colnames(expressiongene2)
  rownames(P_value1)<-colnames(expressiongene2)
  colnames(P_value1)<-colnames(expressiongene2)

P_valuet<-t(P_value1)
P_value2<-P_valuet+P_value
library(reshape2)
library(ggplot2)
P_value2<-as.data.frame(P_value2)
dev.new()
P_value2$Genes<-rownames(P_value2)
P_value3<-melt(P_value2,id.vars=c("Genes"),value.name="pvalue",)

#drawing the heatmap of main and interaction effects
ggplot(P_value3,aes(Genes,variable))+ylab("Genes")+geom_tile(aes(fill=pvalue),
color="white",size=0.01)+theme_classic()+theme(legend.position="top",
legend.justification="right",legend.margin=margin(0,0,0,0),legend.box.margin=margin(0,0,-10,0),
legend.title = element_text(size = 8),legend.text = element_text(size = 6),
legend.key.size =unit(10,"pt"),
axis.text.x=element_text(color="black",angle=90,hjust=1,vjust=1,size=6),
axis.text.y=element_text(color="black",angle=0,hjust=1,vjust=1,size=6))+
scale_fill_gradient(low="light grey",high="red",breaks=seq(0,2.5,0.5),
name=expression(paste("-log"[10],"(p-value)")))+
geom_abline(intercept = 0,slope = 1,color="white")
  
}


% R documentation directory.
\keyword{mixed kernel }% use one of  RShowDoc("KEYWORDS")
\keyword{interaction }% __ONLY ONE__ keyword per line
