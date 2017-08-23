#Diversity calculations
rm(list=ls())
options(stringsAsFactors=FALSE)

FilePath <- "/Users/michellemiron/Dropbox (CUMC)/For Michelle/D233/Pooled 1+2/D233-BM-CD4-TEM-67-6ng.txt"

##TO DO
  #Add Gini coefficient (Gini coefficient, commonly used as a measure of income inequality in economics, was used to assess the inequality of clonotype distribution within a repertoire, source:Thomas PG, Handel A, Doherty PC, La Gruta NL. Ecological analysis of antigen-specific CTL repertoires defines the relationship between naive and immune T-cell populations. Proc Natl Acad Sci U S A. 2013;110:1839-44. doi: 10.1073/pnas.1222149110.)
    

  #Functions to calcualte different diversity metrics:
        #calcC : Compute Clonality
        #calcH: Compute Shannon Entropy
        #calcSI: Compute Simpson's index
        #calcr20 = Compute r20 index
        #calcr50 = Compute r50 index
        #calc JS = Compute Jenson-Shannon divergence
        #calc KL = Compute Kullback-Leibler divergence
            
      ##See this paper for more information on these diveristy metrics:
              #http://www.nature.com/nature/journal/v547/n7661/full/nature22383.html
      #Definitions of these diversity metrics and why you would use them:
            
            ##Shannon Entropy:Shannon's entropy quantifies the uncertainty in predicting the sequence identity of a random sequence from a dataset.1
            ##Clonality: To allow for comparisons between samples differing in the total number of reads, entropy was normalized by division of log 2 of the number of unique productive sequences. Clonality is the reciprocal of normalized Shannon's entropy with values ranging from 0 (most diverse) to 1 (least diverse).   
#References
#1 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4458014/
#2 

# Functions for TCR Diversity

# Functions included:
# Normalize data, Shanon, and JSD
# Functions to add: R20

normalize <- function(data) {
  nc = ncol(data)
  for (i in 1:nc) {
    data[,i] = data[,i] / sum(data[,i])
  }
  return(data)
}

shannon.entropy <- function(p)
{
  if (min(p) < 0 || sum(p) <= 0)
    return(NA)
  p.norm <- p[p>0]/sum(p)
  
  -sum(log2(p.norm)*p.norm)
}

jensen_shannon <- function(p, q){
  #input is columns of matrix with frequency of clones
  ## JSD = H(0.5 *(p +q)) - 0.5H(p) - 0.5H(q)
  # H(X) = \sum x_i * log2(x_i)
  #  p = p[p >0 & q >0]
  #  q = q[p>0 & q>0]
  p = p / sum(p)
  q = q / sum(q)
  Hj = shannon.entropy(0.5 *(p+q)) 
  Hp = shannon.entropy(p) 
  Hq = shannon.entropy(q)
  
  jsd = Hj - 0.5*(Hp+Hq)
  jsd = sqrt(jsd)
  #  cat(Hj, Hp, Hq, jsd, "\n")
  return(jsd)
}


counts2hist<-function(c)
{
  c=c[c>0]
  counts=c
  c=c/sum(c)
  h<-table(c)
  h2=table(counts)
  h<-cbind(names(h),h)
  h2<-cbind(names(h2),h2) 
  mode(h) <- "numeric"
  mode(h2) <- "numeric"
  h<-as.data.frame(h)
  h2<-as.data.frame(h2)
  names(h)<-c("fq","Y")
  names(h2)<-c("counts","Y")
  h$counts=h2$counts
  return(h)
}

getslope<-function(vals)
{
  vals=vals[vals$counts>0,]
  if(length(which(vals$Y==1))==1){
    vals2=vals
  }else{
    if(abs(diff(log(vals$fq[sort(which(vals$Y==1))[1:2]])))<1.5){
      ind=sort(which(vals$Y==1))[2]
    }
    else{
      ind=sort(which(vals$Y==1))[1]
    }
    vals2=vals[1:ind,]
  }
  fq=vals2$fq
  fit=lm(log(vals2$Y)~log(fq))
  slope=coefficients(fit)[[2]]
  return(slope)
}

calcH<-function(vals){
  vals=vals[vals>0]
  fq=vals/sum(vals)
  H=-sum(fq*log2(fq))
  return(H)
}

calcC<-function(vals){
  vals=vals[vals>0]
  H=calcH(vals)
  Hmax=log2(length(vals))
  C=1-H/Hmax
  return(C)
}

calcSI<-function(vals){
  vals=vals[vals>0]
  fq=vals/sum(vals)
  si=sum(fq^2)
  return(si)
}

calcr20 = function(X){
  X=sort(X,decreasing=T)
  X=X[X>0]
  CX=cumsum(X)
  num=length(which(CX/sum(X)<=0.2))
  den=length(X)
  return(num/den)
}

calcr50 = function(X){
  X=sort(X,decreasing=T)
  X=X[X>0]
  CX=cumsum(X)
  num=length(which(CX/sum(X)<=0.5))
  den=length(X)
  return(num/den)
}

# For probability distributions p,q
# Entropy can be obtained from calcH in the attached script.
calcJS<-function(p,q){
  p=p/sum(p)
  q=q/sum(q)
  M=0.5*(p+q)
  JS=entropy(M)-0.5*(entropy(p)+entropy(q))
  
  return(JS)
}

calcKL<-function(p,q)
{
  # normalize
  p=p/sum(p)
  q=q/sum(q)
  
  # compute KL
  Hp=p*log2(p)
  Hpq=p*log2(q)
  f=Hp-Hpq
  f[is.nan(f)]=0
  
  return(sum(f))
}



