### in this file, include the simplified program for parallel computing, also similar to the setting of lm

#' RVS using asymptotic distribution for score test statistic, the variance of score using regular formula for Var_case(E(G_ij|D_ij))
#' when method = Regular (also called likelihood method in the paper), Var_case(G_ij)  is used for Var_case(E(G_ij|D_ij)) when method = RVS.
#'
#' Note, both function scaled the variance in equation (1) in Appendix A by dividing N_case*Ncontrol/N
#' also note the test statistics in anova gives us Rao=s^2/var(s^2), where RVS=s^2/robust var(s^2)
#' so RVS=Rao *var(s^2)/robust var(s^2), s=sum(y_j-mean(y))E(G_ij|D_ij), therefore var(s)//var(E(Gij|Dij)), var(X) in code.
#'
#' @param Y: a vector with the phenotype for n individuals (=ncase+ncont), first ncase (y=1) then ncont (y=0)  (integar)
#' @param X: a vector of genotype for a snp, first ncase and then ncont
#' @param P (vector with length 3  (double)) - estimated genotype frequency for a variant P(G=0), P(G=1) and P(G=2)
#' @param method, has value RVS or regular to indicate the estimation of variance of the score statistic. using the population frequency P
#' to estimate variance for case if method=RVS, otherwise use the expected genotype probability for case.
#' @return   p-values for the snp
#' @export
RVS_asy = function(Y,X,P,method='RVS'){
  if(!method %in% c('RVS','Rvs','rvs','regular','reg','Regular')) {
    cat('Wrong input for test method in RVS_asy, should be RVS or regular!\n')
    return(NULL)
  }
  ## run anlysis for non-NA X only
  Y = Y[!is.na(X)]
  X = X[!is.na(X)]
  a = (glm(Y~1,family='binomial'))
  b = glm(Y~X,family='binomial')
  p = length(X[Y==1])/length(X)
  q = 1 - p
  if(method %in% c('RVS','rvs','Rvs')) {
    v = q*as.numeric(calc_robust_Var(P)) + p*var(X[Y==0],na.rm=TRUE)
  }else{
    v = q*var(X[Y==1]) + p*var(X[Y==0])  ## regular estimation for variance on case at var(X[Y==1]) from test_RVS_asy
  }
  x = anova(a,b,test='Rao')
  x_rvs=x$Rao[2]*var(X)/v  ## var(S_inRao)=#case*#cont/#sample*var(X), var(RVS)=#case*#cont/#sample*(q*var(Xcase)+p*var(Xcont))
  cc = 1-pchisq(x_rvs,1)
  res<-list(coefficients=coef(summary(b))[2,1],
       score=x$Rao[2],
       p_Rao=x[2,6],
       RVS=x_rvs,
       p_rvs=cc)
  return (res)
}


#' use RVS to test associaton by bootrtrap, given phenotype (Y), expected values of genotypes for case and controls (X)
#' estimated genotype frequency (P) and number of times of bootstrap (nboot)
#' @param Y: a vector with the phenotype for n individuals (=ncase+ncont), first ncase (y=1) then ncont (y=0)  (integar)
#' @param X: a vector of genotype for a snp, first ncase and then ncont
#' @param P (vector with length 3  (double)) - estimated genotype frequency for a variant P(G=0), P(G=1) and P(G=2)
#' @param nboot: number of bootstrap
#' @param method, has value RVS or regular to indicate the estimation of variance of the score statistic. using the population frequency P
#' to estimate variance for case if method=RVS, otherwise use the expected genotype probability for case.
#' @return   p-values for the snp
#' @export
RVS_btrap = function(Y,X,P,nboot,method='RVS'){
  if(! method %in% c('RVS','Rvs','rvs', 'regular', 'reg','REG','Reg')) {
    cat('Wrong input for test method in RVS_btrap, should be RVS or regular!\n')
    return(NULL)
  }

  Y = Y[!is.na(X)]
  X = X[!is.na(X)]
  ncase1 = sum(Y==1)
  ncont1 = sum(Y==0)
  p = length(X[Y==1])/length(X)
  q = 1 - p
  if(method %in% c('RVS','Rvs','rvs')) {
    vcase = as.numeric(calc_robust_Var(P))
  } else{
    vcase = var(X[Y==1])
  }
  vcont = var(X[Y==0])
  Tobs = (mean(X[Y==1]) - mean(X[Y==0]))/sqrt(vcase/ncase1+vcont/ncont1)
  X1 = X[Y==1] - mean(X[Y==1])
  X2 = X[Y==0] - mean(X[Y==0])
  C = NULL
  for (j in 1:nboot){
    Xca = sample(X1[],ncase1,replace=TRUE)
    Xco = sample(X2[],ncont1,replace=TRUE)
    vcase = var(Xca)
    vcont = var(Xco)
    C =c(C,(mean(Xca) - mean(Xco))/sqrt(vcase/ncase1+vcont/ncont1))
  }
  cc = sum(abs(C)>=abs(Tobs))/(nboot+1)
  return(cc)
}


#' Use regular Score/Trend test for association for given phenotype (Y) and expected values of genotypes for case and controls (X).
#' it is evaluation by asymptotic distribution
#'
#' @param Y: a vector with the phenotype for n individuals (=ncase+ncont), first ncase (y=1) then ncont (y=0)  (integar)
#' @param X: a vector of genotype for a snp, first ncase and then ncont
#' @return p-values the variant (double)
#' @export
regScore_Rao = function(Y,X){
  ## run anlysis for non-NA X only
  Y = Y[!is.na(X)]
  X = X[!is.na(X)]
  a = (glm(Y~1,family='binomial'))
  b = glm(Y~X,family='binomial')
  x = anova(a,b,test='Rao')
  cc = 1-pchisq(x$Rao[2],1)
  res<-list(est=coef(summary(b))[2,1],
            Rao=x$Rao[2],
            p_Rao=x[2,6])
  return (res)
}

#' Use regular Score/Trend test for association for given phenotype (Y) and expected values of genotypes for case and controls (X)
#' it is evaluation by permutation distribution
#' @param Y: a vector with the phenotype for n individuals (=ncase+ncont), first ncase (y=1) then ncont (y=0)  (integar)
#' @param X: a vector of genotype for a snp, first ncase and then ncont
#' @param nperm: number of permutations
#' @return p-values the variant (double)
#' @export
regScore_perm = function(Y,X,nperm){
    Y = Y[!is.na(X)]
    X = X[!is.na(X)]
    ncase = sum(Y==1)  ## number of case with X!=NA
    ncont = sum(Y==0)  ## number of cont with X!=NA
    test = calc_score_test(X[1:ncase],X[(ncase+1):(ncase+ncont)])
    C = NULL
    for (j in 1:nperm){
      k = sample(length(X))
      X = X[k]
      a = calc_score_test(X[1:ncase],X[(ncase+1):(ncase+ncont)])
      C = c(C,a)
    }
    cc = sum(C>=test)/(nperm+1)
  return (cc)
}


# RVS analysis with rare variants
# Robust Variance Estimate
# Get p-values from RVS with CAST and C-alpha (Resampling: Bootstrap, Variance Estimate: Robust)
# paper mentioned the permutation does not work when have external controls, so use centered X to bootrstap
#' It calls functions \code{calc_ScoreV},\code{calc_ScoreV_RVS},\code{calc_ScoreV_likely},\code{test_CAST}, \code{test_Calpha},\code{sample_bootstrap}
# Input values matrix of expected values of genotypes given sequence data
#' @param Y - phenotype value, 1-cases and 0-controls
#' @param X - matrix of conditional expected values, it has dimensions n by J, J - number of variants, n - sample size
#' @param nboot - number of bootstrap
#' @param P - genotype frequency for J variants P(G=0), P(G=1) and P(G=2), J by 3 matrix
#' @return two p-values for CAST and C-alpha
#' @export
RVS_rare= function(Y,X,P,nboot,method='RVS'){
  if(! method %in% c('RVS','Rvs','rvs', 'regular', 'reg','REG','Reg')) {
    cat('Wrong input for test method in RVS_rare, should be RVS or regular!\n')
    return(NULL)
  }

  X1 = as.matrix(X[Y==1,])
  X2 = as.matrix(X[Y==0,])
  S = calc_ScoreV(X,Y)
  if(method%in%c('RVS','Rvs','rvs'))
  {     Sigma = calc_ScoreV_RVS(X,Y,P)
  }else{
        Sigma = calc_ScoreV_likely(X,Y)
  }
  w = rep(1,length(X[1,]))
  SLobs = test_CAST(S,Sigma,w)
  SQobs = test_Calpha(S,Sigma)
  Q = NULL
  L = NULL

  for (i in 1:nboot){
    Xs = sample_bootstrap(X1,X2)
    #Xa = as.matrix(Xs[Y==1,])
    #Xb = as.matrix(Xs[Y==0,])
    #  X = rbind(Xa,Xb)
    S = calc_ScoreV(Xs,Y)
    Sigma = calc_ScoreV_likely(Xs,Y)
    w = rep(1,length(X[1,]))
    L = c(L,test_CAST(S,Sigma,w))
    Q = c(Q,test_Calpha(S,Sigma))
  }
  pl = sum(L<=SLobs)/(nboot+1)
  pQ = sum(Q<=SQobs)/(nboot+1)
  preturn=data.frame(c(pl,pQ))
  colnames(preturn) <- c('p.CAST','p.Calpha')
  return (preturn)
}


# regScore test for J joint rare variants
# Get p-values from regScore with CAST and C-alpha (Resampling: permutation)
# paper mentioned the permutation does not work when have external controls, so use centered X to bootrstap
#' It calls functions \code{calc_ScoreV},\code{calc_ScoreV_RVS},\code{calc_ScoreV_likely},\code{test_CAST}, \code{test_Calpha},\code{sample_bootstrap}
# Input values matrix of expected values of genotypes given sequence data
#' @param Y - phenotype value, 1-cases and 0-controls
#' @param X - matrix of conditional expected values/genotype calls for J variants, it has dimensions n by J, J - number of variants, n - sample size
#' @param nperm - number of bootstrap
#' @return two p-values for CAST and C-alpha
#' @export
regScore_rare = function(Y,X,nperm){
  S = calc_ScoreV(X,Y)
  Sigma = calc_ScoreV_regVar(X,Y)
  w = rep(1,length(X[1,]))
  SLobs = test_CAST(S,Sigma,w)
  SQobs = test_Calpha(S,Sigma)
  Q = NULL
  L = NULL

  for (i in 1:nperm){
    k = sample(length(Y))
    YY = Y[k]
    S = calc_ScoreV(X,YY)
    Sigma = calc_ScoreV_regVar(X,YY)
    w = rep(1,length(X[1,]))
    L = c(L,test_CAST(S,Sigma,w))
    Q = c(Q,test_Calpha(S,Sigma))
  }
  pl = sum(L<=SLobs)/(nperm+1)
  pQ = sum(Q<=SQobs)/(nperm+1)
  return (c(pl,pQ))
}

