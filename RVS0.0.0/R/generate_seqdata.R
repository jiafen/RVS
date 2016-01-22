# Here include two functions, generate sequence data under null and alternative hypothesis
# Note to Andriy: The seq_call() function that is called from  generate_seqdata_alt(), doesn't have input args: gen_freq and R.
# I commented these two args out for now.  Also, gen_vect is an empty arg; see the seq_call() code.  --Ted, Agreed by Jiafen


#' Generate sequence data under null hypothesis, it calls functions \code{seq_call},\code{calc_EG_general},\code{calc_EM},\code{calc_pobs_ndepth}
#' \code{calc_Mr}
#'
#' @param nsnp, number of variants (integer) >0
#' @param ncase, number of cases (integer) >0
#' @param ncont, number of controls (integer) >0
#' @param mdcase, average read depth in cases (double) >0
#' @param sddcase, standard deviation of read depth in cases  (double) >0
#' @param mdcont, average read depth in controls (double) >0
#' @param sddcont, standard deviation of read depth in controls (double) >0
#' @param em, average error rate, probability that sequence call is wrong (double) >0
#' @param esd, standard deviation of error rate (double) >0
#' @param mmaf, a vector of MAF for nsnp variants, HWE is assumed >0
#'
#' @return MM - conditional expected value E(Gij|Dij)
#' @return P -estimated genotype frequencies P(G=0), P(G=1), P(G=2)
#' @return G - true genotype from which sequence data generated
#' @export
generate_seqdata_null <- function(nsnp, ncase, ncont, mdcase, sddcase, mdcont, sddcont, em, esd, mmaf){
  MAF = NULL
  A1 = 'C'
  A2 = 'T'
  MM = NULL
  G = NULL
  P = NULL
  for (i in 1:nsnp){
    if (i %% 10 == 1) {cat('i=',i,'\n')}
    v = NULL
    erv = NULL
    rdv = NULL
    gv = NULL
    maf = mmaf[i]
    MAF = c(MAF,maf)
    for (j in 1:ncase){
      rd = round(sddcase*rnorm(1) + mdcase)
      if (rd <= 0){rd=1}
      error = esd*rnorm(rd) + em
      k = rbinom(1,2,maf)
      gv = c(gv,k)
      gen_vect = c('CC','CT','TT')
      if (k==2){genotype=c('C','C')}
      if (k==1){genotype=c('C','T')}
      if (k==0){genotype=c('T','T')}


      a = seq_call(genotype, error, ndepth=rd)
      v = c(v,a)
      erv = c(erv,error)
      rdv = c(rdv,rd)
    }


    for (j in 1:ncont){
      rd = round(sddcont*rnorm(1) + mdcont)
      if (rd <= 0){rd=1}
      error = esd*rnorm(rd) + em
      k = rbinom(1,2,maf)
      gv = c(gv,k)
      gen_vect = c('CC','CT','TT')
      if (k==2){genotype=c('C','C')}
      if (k==1){genotype=c('C','T')}
      if (k==0){genotype=c('T','T')}
      a = seq_call(genotype, error,ndepth=rd)
      v = c(v,a)
      erv = c(erv,error)
      rdv = c(rdv,rd)
    }


    G = cbind(G,gv)
    t_ndepth=length(erv)
    M  = calc_Mr(erv,v,t_ndepth,rdv)
    Mm = calc_pobs_ndepth(erv,v,t_ndepth)
    p = calc_EM(M)
#    cat('i=',i,',p=',p,',sum_p=',sum(p),'\n')
    p[p<0]=0
####
    #p = c(0.9^2,2*0.9*0.1,0.1^2)
    P = rbind(P,p)
#     if ((sum(p< (-0.00001))>0) || (sum(p>1.00001)<0)){
#       cat('Problem in seqdata generating under Null Hypothesis!!!')
#       cat('for the', i,' th SNPs!\n')
#       cat('maf=', maf, 'k=',k,'\n')
#       return (NULL)
#     }
    EG = calc_EG_general(Mm,p,rdv)
    MM = cbind(MM,EG)
  }
  return (list(MM=MM,P=P,G=G))
}

#
#' Generate sequence data under alternative hypothesis

#' it calls functions \code{seq_call},\code{calc_EG_general},\code{calc_EG_Var}, \code{calc_EM}, \code{calc_pobs_ndepth}
#' \code{calc_Mr}
#' @param mmafCa MAF for nsnp variants in cases, HWE is assumed >0
#' @param mmafCo MAF for nsnp variants in controls, HWE is assumed >0
#' @param all the other params are describe in \code{generate_seqdata_null}, here MAF for cases/controls are different.
#' @return For alt, outputs same object with additional:
#'
#' @return MM - conditional expected value E(Gij|Dij)
#' @return P -estimated genotype frequencies P(G=0), P(G=1), P(G=2)
#' @return G - true genotype from which sequence data generated
#' @return MM2-  variance  of the conditional expected value E(Gij|Dij)
#generate_seqdata_alt <- function(nsnp,ncase,ncont,mdcase,sddcase,mdcont,sddcont,em,esd,mmafCa,mmafCo,R){
#' @export
generate_seqdata_alt <- function ( nsnp, ncase, ncont, mdcase, sddcase, mdcont, sddcont, em, esd, mmafCa, mmafCo) {
  MAFCa = NULL
  MAFCo = NULL
  A1 = 'C'
  A2 = 'T'
  MM = NULL
  MM2 = NULL
  G = NULL
  P = NULL
  for (i in 1:nsnp){
    if (i %% 10 == 1) {cat('i=',i,'\n')}
    v = NULL
    erv = NULL
    rdv = NULL
    gv = NULL
    maf = mmafCa[i]
    MAFCa = c(MAFCa,maf)
#    cat('mafca=',maf,'\n')
    for (j in 1:ncase){
      rd = round(sddcase*rnorm(1) + mdcase)
      if (rd <= 0){rd=1}
      error = esd*rnorm(rd) + em
      k = rbinom(1,2,maf)
      gv = c(gv,k)
      gen_vect = c('CC','CT','TT')
      if (k==2){genotype=c('C','C')}
      if (k==1){genotype=c('C','T')}
      if (k==0){genotype=c('T','T')}
      # a = seq_call(genotype, error, gen_vect, gen_freq,R)
      a=seq_call(genotype, error, ndepth=rd)
      v = c(v,a)
      erv = c(erv,error)
      rdv = c(rdv,rd)
    }

    maf = mmafCo[i]
    MAFCo=c(MAFCo,maf)
#    cat('mafco=',maf,'\')
    for (j in 1:ncont){
      rd = round(sddcont*rnorm(1) + mdcont)
      if (rd <= 0){rd=1}
      error = esd*rnorm(rd) + em
      k = rbinom(1,2,maf)
      gv = c(gv,k)
      gen_vect = c('CC','CT','TT')
      if (k==2){genotype=c('C','C')}
      if (k==1){genotype=c('C','T')}
      if (k==0){genotype=c('T','T')}
      # a = seq_call(genotype, error, gen_vect, gen_freq,R)
      a = seq_call(genotype, error, rd)
      v = c(v,a)
      erv = c(erv,error)
      rdv = c(rdv,rd)
    }

    G = cbind(G,gv)

    t_ndepth=length(erv)
    M  = calc_Mr(erv,v,t_ndepth,rdv)
    Mm = calc_pobs_ndepth(erv,v,t_ndepth)
    p = calc_EM(M)
    #p = c(0.9^2,2*0.9*0.1,0.1^2)
   # cat('i=',i,',p=',p,',sum_p=',sum(p),'\n')
    p[p<0]=0
    P = rbind(P,p)
#     if ((sum(p< (-0.00001))>0) || (sum(p>1.00001)<0)){
#       cat('Problem in seqdata generating under alternative Hypothesis!!!')
#       cat('for the', i,' th SNPs!\n')
#       return (NULL)
#     }

    EG = calc_EG_general(Mm,p,rdv)
    EG2 = calc_EG_Var(Mm,p,rdv)
    MM = cbind(MM,EG)
    MM2 = cbind(MM2,EG2)
  }
  return (list(MM=MM,P=P,G=G,MM2=MM2))
}

