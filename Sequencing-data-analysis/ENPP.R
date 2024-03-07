
#################################
#################################
### R script for application of ENPP 
#################################
### example data set generation 
### 1000 markers (assumed to be SNP markers) in total : 10 markers are significant, 990 markers are not significant 
### 5000 individuals 
### 10,000 permutation rounds for pruning step, and 39,999 permutations for unpruned markers 
### total number of permutations for unpruned markers is 50,000(=10,000+39,999+1(=original data analysis))  
#################################
#################################

library(mvtnorm)
set.seed(123)

sigma1<- matrix(0.2, nrow=11, ncol=11); diag(sigma1)<- 1 ; sigma1<-sigma1*10 
sig1<- rmvnorm(5000, mean=rep(100,11), sigma=sigma1)
y1<-sig1[,1]   #### response variable 

mat1<- matrix(NA, nrow=5000, ncol=10)

### marker generation with maf=0.1 and Hardy-Weinberg principle for significant markers 

for(i in 1:10)
 { 
   rank1<-rank(sig1[,i+1])
   mat1[which(rank1<=4050),i]<-0      
   mat1[which(rank1 > 4050 & rank1 <= 4950),i]<-1  
   mat1[which(rank1 >= 4950),i]<-2   
 }

### marker generation with mafs from a uniform distribution for insignificant markers 

mat2<- matrix(NA, nrow=5000, ncol=990)
for(i in 1:990)
 { 
   mat2[,i]<-rbinom(n=5000,size=2,prob=runif(1,min=0.05, max=0.50))
 }

mat_tot<-cbind(mat1,mat2) #### explanatory variables in total 
res1<-matrix(NA,nrow=1000, ncol=4)  #### linear model results matrix 
colnames(res1)<-c("estimate","se","t-value","p-value") 

for(i in 1:1000)
{ 
 res1[i,]<- summary(lm(y1~mat_tot[,i]))$coef[2,] 
} 

res1[1:20,]  ### first 20 results of 1,000 markers with parametric approaches 

#############################
#############################
### setting p-adjust, p-prun, and c-prun values, introduced in the article. 
#############################
#############################

n1=10^3
p1 = 0.05
p1adj = p1/n1   ##### Bonferroni type significance level 
p1prun <- p1adj 
coun1<-c(rep(NA,10^4))  ### permutation round is assumed to be 10,000 for pruning process 

for(i in 1:(10^4))
{ 
  p1val<- as.numeric(as.character(1-pbinom(0:min(i-1,20) ,i,p1adj)))
  coun1[i]<-which(signif(p1val,10)<=signif(p1prun,10))[1]
}

table(coun1)  #### coun1 is cprun for each permutation round 

#############################
#############################
### pruning step with 10,000 permutation rounds (assume that we are using only one thread) 
#############################
#############################

set.seed(123)

loca1<-c(1:1000) 
res2<-matrix(NA,nrow=1000,ncol=10^4)  
stat1<-abs(res1[,3]) 
sum1<-0 

for(i in 1:10000)
{ 
  y1n<-sample(y1) 
  res3<-matrix(NA, nrow=1000,ncol=4)
  for(j in  loca1)
   {
     res3[j,] <- summary(lm(y1n~mat_tot[,j]))$coef[2,] 
   }
  print(c("i",i,"regression done")) 
  loca2<-which(stat1<=abs(res3[,3]))   ### for case of having larger statistics than the original statistics.  
  loca3<-setdiff(loca1,loca2)   ### for case of having smaller statistics than the original statistics.  
  res2[loca2,i]<-1 
  res2[loca3,i]<-0  
  sum1<-sum1+res2[,i] 
  loca4<-which(sum1>=coun1[i]) 
  loca1<-setdiff(loca1,loca4) 
  print(c("i",i,"number of survived",length(loca1))) 
}

length(loca1) ; loca1   #### SNPs survived from the approach  
mat_tot2<-mat_tot[,loca1]  #### dataset of survived markers  (unpruned) 

#############################
#############################
### final evaluation step (39,999(=50,000-10,001) permutaion is additionally applied for each unpruned marker, to obtain final precise p-value) 
### parallel computing (multicore approach) can be applied 
#############################
#############################

set.seed(1234)

res4<-matrix(NA, nrow=39999,ncol=length(loca1))
for(i in 1:length(loca1))
 {
   for(j in 1: 39999)
    {
        y1n<-sample(y1) 
        res4[j,i]<- abs(summary(lm(y1n~mat_tot2[,i]))$coef[2,3])  
        if(j%%1000==0){print(c("i",i,"j",j))}  
    }
 }

res5<-matrix(NA,nrow=1000,ncol=2)  ### results of pruning step 
res6<-matrix(NA,nrow=length(loca1),ncol=2)  ### results of post-pruning step with 10 unpruned markers 

for(i in 1:1000)
 { 
   res5[i,1]<- sum(!is.na(res2[i,]))
   res5[i,2]<- sum(res2[i,],na.rm=T)
 }

for(i in loca1)
 { 
   res6[i,1]<- 39999
   res6[i,2]<- sum(stat1[loca1[i]]<=res4[,i])
 }

res7<-res5   #### combining permutation results of pruning step and post pruning step 
res7[loca1,]<-res5[loca1,]+res6 ; res7<- res7+1  

res_fin<- cbind(res1,res7,res7[,2]/res7[,1])   ### final results 
colnames(res_fin)<- c("ori_estimate","ori_se","ori_t-value","ori_p-value","permut_times","exceeding_times","permut_p-value")

################################################################








