### In this file, we include all the helper function for association test
## calc_score_test: the regular score test, including the regular variance of score statistic
## calc_robust_Var: the robust variance of score statistic (using prob(G_ij) for case)
## calc_scoreV: the vector of score statistic for J rare variants
## calc_scoreV_regVar(X,Y): the regular variance of score vector.
## calc_scoreV_RVS(X1,X2,Y): the robust variance of score vector (using p(G_ij) for case)
## calc_scoreV_RVS_regVar(X1,X2,Y): the robust variance of score vector using var(X(Y=1)) and var(X(Y=0))
## test_CAST, test_Calpha and test_Hotelling, the three test statistic for rare variants with input score S and variance Sigma
## calc_centralize_matrix, sample_boostrap used in the rare variant associaiton test.

# These packages should be downloaded from R-CRAN
library('MASS')
library('CompQuadForm')

#' clac_score_test calculates the score test for given expected value of genotype  for case and controls seperately (M1, M2)
#'
#' calc_score_test is called in \code{\link{test_regScore_perm}}
#' @param M1: the expected value of genotype for case E(G_ij|D_ij) on one snp (dimension: ncase*1)
#' @param M2: the expected value of genotype for control E(G_ij|D_ij) on one snp (dimension: ncont*1)
#' @return :the score test statistic calculated from M1 and M2
#' T=S^2/var(S), for a given j, S_j=sum_i((Y_i-bar(Y))E(G_ij|D_ij))=sum_case(1-bar(Y))E(G|D)-sum_cont(bar(Y)E(G|D))
calc_score_test = function(M1,M2){
  X = c(M1,M2)
  Y = c(rep(1,length(M1)),rep(0,length(M2)))
  p = length(M1)/length(X)
  q = 1 - p
  S = q*sum(M1)-p*sum(M2)
  vs = p*q*(length(X))*var(X)
  x = S^2/vs
  return (x)
}



#' clac_score_test calculates the score test for given expected value of genotype  for case and controls seperately (M1, M2)
#'
#' calc_score_test is called in \code{\link{test_regScore_perm}}
#' @param M1: the expected value of genotype for case E(G_ij|D_ij) on one snp (dimension: ncase*1)
#' @param M2: the expected value of genotype for control E(G_ij|D_ij) on one snp (dimension: ncont*1)
#' @return :the score test statistic calculated from M1 and M2 and empirical variance estimation seperated for case and controls.
#' T=S^2/var(S), for a given j, S_j=sum_i((Y_i-bar(Y))E(G_ij|D_ij))=sum_case(1-bar(Y))E(G|D)-sum_cont(bar(Y)E(G|D))
calc_score_test1 = function(M1,M2){
  X = c(M1,M2)
  ncase <- length(M1)
  ncont <- length(M2)
  Y = c(rep(1,ncase),rep(0,ncont))

  p = ncase/length(X)
  q = 1 - p
  S = q*sum(M1)-p*sum(M2)
  vs = p*ncont*(q*var(X[Y==1])+p*var(X[Y==0]))
  x = S^2/vs
  return (x)
}


#'  the robust variance estimation for case
#'
#'   Variance of genotypes Var(X)=E(X^2)-E(X)^2
#' @param P: the genotype frequency, calcualted from EM algorithm
#' @return the Robust variance of E(Gij|Dij)
#' called by \code{test_RVS_asy},\code{test_RVS_bootstrap}, and \code{calc_ScoreV_RVS}
calc_robust_Var = function(P){
  Sq = 4*P[3] + P[2]
  Sm = 2*P[3] + P[2]
  S = Sq - Sm^2
  return (S)
}






#' calc_ScoreV gets score vector S for nsnp number of rare variants
#'
#' it is called in \code{test_rare_regScore}, \code{test_rare_RVS} and \code{test_rare_RVS_regVar}
#'@param X: the expected conditional genotype probability
#'@param Y: the phenotype
#' @return the score statistic vector S for nsnp rare variant
calc_ScoreV = function(X,Y){ # NEED to check!!!
  L = length(X[1,])
  S = NULL
  for (i in 1:L){
    Yn = Y[!is.na(X[,i])]
    Xn = X[!is.na(X[,i]),i]
    xca = Xn[Yn==1]
    xco = Xn[Yn==0]
    my = mean(Yn)
    s = sum(xca,na.rm=T)*(1-my) - my*sum(xco,na.rm=T)
    S = c(S,s)
  }
  S = t(S)
  return (S)
}



#' calc_ScoreV_regVar: Get Variance of vector S by regular method
#'
#'  it is called in \code{test_rare_regScore}
#' @param X: the expected genotype probability
#' @param Y: the phenotype
#' @return the variance of score statistic vector S by regular method
calc_ScoreV_regVar = function(X,Y){

  X1 = X[Y==1,]
  X2 = X[Y==0,]
  l1 = length(X1[,1])
  l2 = length(X2[,1])
  J = length(X1[1,])

  a =colSums(is.na(X1),na.rm=T)
  b =colSums(is.na(X2),na.rm=T)

  ncase = rep(l1,J) - a
  ncont = rep(l2,J) - b
  nn = ncase+ncont
  L= length(Y)
  Xm = cov(X,use="pairwise.complete.obs")
  vs  = sqrt(diag(ncase*ncont/nn^2))
  diag_S  = vs%*%Xm%*%vs*L
  return (diag_S)
}

#' calc_ScoreV_RVS_regVar: Get Variance of vector S by RVS using regular expected prob
#'
#' it is called in \code{test_rare_RVS_regVar}
#'@param X1: the expected conditional genotype probability for case
#'@param X2: the expected conditional genotype probability for control
#'@param Y: the phenotype
#'@return the variance of score statistic vector S by RVS using regular expected prob
calc_ScoreV_RVS_regVar = function(X1,X2,Y){
  l1 = length((X1[,1]))
  l2 = length(X2[,1])
  J = length(X1[1,])

  a =colSums(is.na(X1),na.rm=T)
  b =colSums(is.na(X2),na.rm=T)

  ncase = rep(l1,J) - a
  ncont = rep(l2,J) - b
  nn = ncase+ncont

  Yhat = mean(Y)
  L= length(Y)
  p = l2/l1
  q = l1/l2
  Xm1 = cov(X1,use="pairwise.complete.obs")*(l1-1)
  Xm2 = cov(X2,use="pairwise.complete.obs")*(l2-1)

  vs  = sqrt(diag(ncase*ncont/nn^2))
  diag_S  = vs%*%(p*Xm1+ q*Xm2)%*%vs
  return (diag_S)
}


#' helper6: Get Variance of vector S by RVS using robust variance
#'
#'  it is called in \code{test_rare_RVS}
#'@param X1: the expected conditional genotype probability for case
#'@param X2: the expected conditional genotype probability for control
#'@param Y: the phenotype
#'@param P - genotype frequency for J variants P(G=0), P(G=1) and P(G=2), J by 3 matrix
#'
#' @return the variance of score statistic vector S by  RVS using robust variance
calc_ScoreV_RVS = function(X1,X2,P,Y){
  L = length(P[,1])
  V = rep(NA,length(L))
  for (i in 1:L){
    V[i] = calc_robust_Var(P[i,])
  }
  V = diag(sqrt(as.numeric(V)))
  l1 = length(X1[,1])
  l2 = length(X2[,1])
  Yhat = mean(Y)
  L= length(Y)
  p = l2/l1
  q = l1/l2
  Sgcase = t(V)%*%cor(X1,use="pairwise.complete.obs")%*%V*l1
  Xm2 = calc_centralize_matrix(X2)
  Xm2[is.na(Xm2)]<-0  ## if not set NA to 0, t(Xm2)%*%Xm2 will have too many NAs
  vs  = var(Y)
  diag_S  = vs*(p*Sgcase + q*t(Xm2)%*%Xm2)
  return (diag_S)
}

#
#' helper7: the linear test for nsnp rare variant
#'
#' it is called in \code{test_rare_regScore}, \code{test_rare_RVS} and \code{test_rare_RVS_regVar}
#' @param S: vector of length nsnp, the Score for each variant
#' @param Sigma: matrix Sigma is the variance matrix of vector S.
#' @param w: the weight for CAST test
#' @param nsnp: number of snp
#' @return p-values from CAST test for J variants (double)
#' @export
test_CAST= function(S,Sigma,w,nsnp){
  T = (w)%*%t(S)
  Sg = t(w)%*%Sigma%*%w
  linear = as.numeric(T)/sqrt(as.numeric(Sg))
  p = 2*(1-pnorm(abs(linear)))
  return (p)
}


#' helper8: the quadratic test for nsnp rare variant (Calpha)
#'
#' it is called in \code{test_rare_regScore}, \code{test_rare_RVS} and \code{test_rare_RVS_regVar}
#' @param S: vector of length nsnp, the Score for each variant
#' @param Sigma: matrix Sigma is the variance matrix of vector S.
#' @param nsnp: number of snp
#'@return p-values from Calpha test for J variants (double)

test_Calpha=function(S,Sigma){
  quad  =  (S)%*%t(S)
  lambda = Re(eigen(Sigma)$values)  ## the real part of the eigenvalue of Sigma will be the weight
  qp = davies(quad,lambda)   ### in CompQuadForm library, to calculate P[Q>quad] using davies's method
  return (abs(qp$Qq))
}

#'  helper9:  Hotelling test
#' Input:'
#' @param S: vector of length nsnp, the Score for each variant
#' @param Sigma: matrix Sigma is the variance matrix of vector S.
#' @param nsnp: number of snp

#' @return p-values from Hotelling test for J variants (double)
test_Hotelling = function(S,Sigma){
  dd  =  try(ginv(Sigma), silent = TRUE)
  if ( class(dd) == 'try-error'){
    cat('Inverse_error','\n')
    Sigma = diag(diag(Sigma))
  }
  X =ginv(Sigma)
  quad = S%*%X%*%t(S)
  rank = qr(X)$rank
  p = 1-pchisq(quad,rank)
  return (p)
}

#' helper10: centralize matrix X by column
#'
#'  called by \code{\link{sample_bootstrap}} and \code{cal_ScoreV_RVS}
#' @param X: the genotype likelihood matrix, one column for each snp
#' @return a volumn-centeralized matrix for the given input matrix X
calc_centralize_matrix = function(X){
  l = length(X[1,])
  mX = NULL
  for (i in 1:l){
    mm = mean(X[,i],na.rm=T)
    mv = X[,i]-mm
    mX = cbind(mX,mv)
  }
  return (mX)
}


#
#' helper11: Bootstrap resampling
#'
#' it calls \code{calc_centralized_matrix} and called by \code{\link{test_rare_RVS}}, \code{\link{test_rare_RVS_regVar}}.
#'
#' @param X1: genotype likelihood matrix for case
#' @param X2: genotype likelihood matrix for control
#' @return X: bootstrap for genotype likelihood matrix for case and control
sample_bootstrap = function(X1,X2){
  case = length(X1[,1])
  cont = length(X2[,1])
  X1 = calc_centralize_matrix(X1)
  X2 = calc_centralize_matrix(X2)
  ca = sample(1:case,case,replace=TRUE)
  co = sample(1:cont,cont,replace=TRUE)
  Xca = as.matrix(X1[ca,])
  Xco = as.matrix(X2[co,])
  X = rbind(Xca,Xco)
  #return (X)
  return(list('Xca'=Xca,'Xco'=Xco))
}

##
## helper12: to find the index of SNPs that are monormophic in case/controls
#' Input
#' @param subgroup: a matrix for conditional expected genotype or genotype calls for case or controls. dimension nsample*nsnp
#' @param nsnp: number of snps
#' @param type: specify it is conditional expected genotype or genotype call itself
#' @return  a vector indicate the snp columns that are monormophic.
calc_homo<-function(subgroup,nsnp,type)
{  hom=NULL;
   if(type %in%c('conditional','expected','exp'))
   {
     vrep <- 1:nsnp
   }else{
     vrep <- 7:(nsnp+6)
   }
   for (i in vrep)
   {
    tmp1=length(levels(as.factor(round(subgroup[!is.na(subgroup[,i]),i],1))))
    if(tmp1==1) hom=c(hom,i)
   }
    return(hom)
  }


  ##
## helper13: to remove the SNPs that are monormophic in case/controls
#' Input
#' @param case: a matrix for conditional expected genotype or genotype calls for case. Dimension ncase*nsnp
#' @param cont:  a matrix for conditional expected genotype or genotype calls for controls. Dimension ncont*nsnp
#' @param snp: a matrix for snps, including chr, position, expected sample frequency if the case and controls are conditional expected genotype
#' @return  a vector indicate the snp columns that are monormophic.
#' @export
rm_hom <-function(case,cont,snp,type='expected')
{
  if(!type%in%c('conditional','expected','exp','genotype','hard','calls','call'))
  {
    cat('Wrong data type for remove monormophic snps! need be expected or genotype!\n')
    return(NULL)
}
  hom1<-calc_homo(case,nrow(snp),type)
  hom2<-calc_homo(cont,nrow(snp),type)
  hom<-union(hom1,hom2)
  case=case[,-hom]
  cont=cont[,-hom]
  if(type %in% c('conditional','expected','exp'))
  {
    snp1=snp[-hom,3:5]
    return(list('case'=case,'cont'=cont,'snp'=snp,'P.rare'=snp1,'hom'=hom))
  }else{
    hom1<-hom-6
    snp1=snp[-hom1,]
    return(list('case'=case,'cont'=cont,'snp'=snp1,'hom'=hom1))
  }


}

