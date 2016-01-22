#' seq_call produces the base call from given a genotype and the error scores at one locus.
#'
#' It is called from the \code{generate_seqdata_null} and \code{generate_seqdata_alt} functions.
#' @param genotype: vector of two alleles for a given genotype, ie. c('C', 'T')
#' @param error: vector of base call errors: Phen scores (double).  Length of this vector is the read depth.
#' @param ndepth: the read depth (integer)
#' @return vect_reads: vector of row sequence reads ('A','G','T','C') for single variant. The length of this vector is the read depth.
#' @export
seq_call <- function(genotype, error, ndepth){
  if(length(error) != ndepth) {
    cat ('dimension error in seq_call!\n')
    return(NULL)
  }
  vect_reads = NULL
  all_alleles = c('A','G','T','C')
  for (i in 1:ndepth){
    s = rbinom(1,1,0.5)+1
    g = genotype[s]
    ng = all_alleles[all_alleles!=g]
    e3 = error[i]/3
    e23 = e3*2
    e = e3*3
    r = runif(1)

    if (r <= e3) {
      vect_reads <- c(vect_reads,ng[1])
    }else if (r > e3 & r <= e23) {
      vect_reads <- c(vect_reads,ng[2])
    }else if (r > e23 & r <= e) {
      vect_reads <- c(vect_reads,ng[3])
    }else{
      vect_reads <- c(vect_reads,g)
    }
  }
  return (vect_reads)
}


## helper3
#' calc_psingle computes the genotype likelihood  P(D=g|G=G1G2,e) for a single base read given the sequencing error rate and true genotype.
#'
#' calc_psingle is called from the \code{calc_pobs_ndepth} and \code{calc_Mr} functions.
#' please see Appendix B for how to calculate, P(D=g|G=G1G2)=1/2*P(D=g|G1)+1/2*P(D=g|G2) and  p(D=g|G1=g)=1-e, p(D=other letter|G1=g)=e/3.
#' @param sing_read: single read, ie. ('A')
#' @param A1: allele 1 in the true genotype, eg. (A1-'A',A2-'C')
#' @param A2:  allele 2 in the true genotype, ie. (A1-'A',A2-'C')
#' @param error: sequencing error (double)
#' @return p - likelihood of observing sing_read given genotype (A1A2) and error rate
#' @seealso \code{\link{calc_pobs_ndepth}} or \code{\link{calc_Mr}}
##

calc_psingle = function(sing_read,A1,A2,error){
  if (sing_read == A1){
    p1 = 1-error
  }else{
    p1 = error/3
  }

  if (sing_read == A2){
    p2 = 1-error
  }else{
    p2 = error/3
  }
  p = 1/2*p1 + 1/2*p2

  return (p)
}

#' calc_pobs_ndepth:  Get likelihood P(Di|G) matrix given a vector of error rates and a vector of sequence reads.
#'
#' This function calculates a likelihood matrix given a vector of error rates and a vector of sequence reads.
#'  It is called from the \code{generate_seqdata_null} and \code{generate_seqdata_alt} functions.
#' s  The output of this function  are used as inputs to \code{calc_EG}.
#' it use function \code{calc_psingle} to calculate the likelihood for each read and one possible genotype (TT, CT or CC)
#' Note to Andriy:  Why is LG hardcoded to LG=c('TT','CT','CC')  ??  Guess it's to denote the AA, Aa, aa.
#'
#' @param Error: vector of base call error rates; Uses the same \code{error} input in the \code{seq_call} function
#'               but concatenate all \code{error} values over all samples.
#' @param vect_row_reads: vector of row sequence reads ('A','G','T','C') outputted from function \code{seq_call},
#'                        but concatenated over all samples in same fashion as \code{Error}.
#' @param ndepth: integer, read depth, should be length(error)=length(vect_row_reads)=ndepth.
#' @return M matrix of P(read|AA), P(read|Aa) and P(read|aa), where row dimensions equal \code{length(Error)}=length(vect_row_reads), by 3 columns.
#' @seealso \code{\link{seq_call}}
calc_pobs_ndepth= function(Error,vect_row_reads,ndepth=1){
  L = length(Error)
  L2=length(vect_row_reads)
  if(L!=ndepth || L2!=ndepth) {
    cat("Lengths mismatch between error, reads and read-depth in cal_pobs_ndepth!");
    return(NULL)
  }

  M = NULL
  LG  =c('TT','CT','CC')

  for (i in 1:L){
    m = NULL
    for  (j in 1:3){
      G = LG[j]
      A1 = substring(G,1,1)
      A2 = substring(G,2,2)
      m = c(m,calc_psingle(vect_row_reads[i],A1,A2,Error[i]))
    }
    M = rbind(M,m)
  }
  return (M)
}


#' Andriy's original get_Mr, Get likelihood P(Di|G) for more general case,
#'
#' It is used in \code{generate_seqdata_null} and \code{generate_seqdata_alt}
#' @describeIn calc_pobs_ndepth with additional param \code{rdv}
#' @param rdv: vector of read depths concatenated for all samples (integer).
#' @param Error: vector of base call error rates; Uses the same \code{error} input in the \code{seq_call} function
#'               but concatenate all \code{error} values over all samples.
#' @param vect_row_reads: vector of row sequence reads ('A','G','T','C') outputted from function \code{seq_call},
#'                        but concatenated over all samples in same fashion as \code{Error}.
#' @param ndepth: integer, read depth, length(error)=length(vect_row_reads)=ndepth.
#' @return M matrix of P(read|AA), P(read|Aa) and P(read|aa), where row dimensions equal \code{length(Error)}=length(vect_row_reads), by 3 columns.
#' Looks like getting likelihood over all reads cumulatively? Andriy, please explain.

calc_Mr = function(Error,vect_row_reads,ndepth=1,rdv){
  L = length(Error)
  L2=length(vect_row_reads)
  if(L!=ndepth || L2!=ndepth) {cat("Lengths mismatch between error, reads and read-depth in calc_Mr!"); return(NULL)}

  L = length(rdv)
  M = NULL
  LG  =c('TT','CT','CC')
  S = 0
  for (i in 1:L){
    m = NULL
    for  (j in 1:3){
      LL = 1
      for (kk in 1:rdv[i]){
        G = LG[j]
        A1 = substring(G,1,1)
        A2 = substring(G,2,2)
        LL = LL*calc_psingle(vect_row_reads[S+kk],A1,A2,Error[S+kk])
      }
      m =c(m,LL)
    }
    S =S+rdv[i]
    M = rbind(M,m)
  }
  return (M)
}



## helper5
#' calc_EG calculates the conditional expected value of the genotype given the genotype likelihoods and frequencies E(Gij|Dij)
#'
#' It is called from the \code{generate_seqdata_null} and \code{generate_seqdata_alt} functions.
#' talk to Andriy, he said in the kk loop, he did it for more general case, we only need rdv=1 for our case.
#' so I renamed Andriy's get_EG to calc_EG_general, delete the kk loop in calc_EG function.
#'
#' @param M: genotype likelihoods AA, Aa, aa, matrix nsample by 3 (double), genotype likelihood for one locus in VCF;
#' uses output from \code{calc_pobs_ndepth} function for simulation data and output from \code{getgenexp} or \code{getMAF} for VCF input
#' @param p: genotype frequencies AA, Aa, aa (double) for each SNP, vector of length 3; output from \code{calc_EM} function.
#' @param rdv: read depth (vector of integers) for all samples, should be a vector of 1 with length=nsample.
#' @return EG: a vector containing conditional expectation values. Size should be equal to length(rdv).
#'
#' E(Gij|Dij)=sum_(g=0)^2 gP(G_ij=g|D_ij)
#' P(G_ij=g|D_ij)=P(D_ij|G_ij=g)*P(G_ij=g)/P(D_ij),
#' where P(D_ij|G_ij=g) (input M) are from the VCF file or \code{calc_pobs_ndepth}, p(G_ij=g) (input p) are from \code{calc_EM}.
#' all values are rescaled by P(D_ij), (see m below, but it is not probability, so m/sum(m) gives us sth. like prob.)

calc_EG <- function(M, p, rdv) {
  LL = length(rdv)

  EG = NULL
  g = c(0,1,2)
  for (i in 1:LL){
    m = NULL

    ### m=(M[i,1]*p[1],M[i,2]*p[2],M[i,3]*p[3])
    for  (j in 1:3){
      L=M[i,j]
      m = c(m,L*p[j])
    }
    pm = sum(m/sum(m)*g)
    EG = c(EG,pm)
  }
  return (EG)
}

#' Andriy's original get_EG function, more general case.
#'
#'  It is used in \code{generate_seqdata_null} and \code{generate_seqdata_alt}
#'  The more general case for \code{calc_EG}
#'
#' @param M: genotype likelihoods AA, Aa, aa, matrix sum(rdv) by 3 (double);
#' uses output from \code{calc_pobs_ndepth} function for simulation data and output from \code{getgenexp} or \code{getMAF} for VCF input
#' @param p: genotype frequencies AA, Aa, aa (double); output from \code{calc_EM} function.
#' @param rdv: read depth (vector of integers) for all samples, should be a vector of length=nsample.
#' @return EG: a vector containing conditional expectation values. Size should be equal to length(rdv).
#'
calc_EG_general= function(M,p,rdv){
  LL = length(rdv)
  S = 0
  EG = NULL
  g = c(0,1,2)
  for (i in 1:LL){
    m = NULL
    for  (j in 1:3){
      L = 1
      for (kk in 1:rdv[i]){

        L = L*M[S + kk,j]
      }
      m = c(m,L*p[j])

    }
    S = S + rdv[i]
    pm = sum(m/sum(m)*g)
    EG = c(EG,pm)
  }
  return (EG)
}

## helper 6
#' This function calculates the variance of conditional expected value of the genotype given the genotype likelihoods and frequencies Var(E(Gij|Dij)).
#'
#'  It is called from the \code{generate_seqdata_alt} functions.
#' talk to Andriy, he said in the kk loop, he did it for more general case, we only need rdv=1 for our case.
#' so I renamed Andriy's get_EG2 to \code{calc_EG_var_general}, and delete the kk loop in this function.
#' @describeIn calc_EG also, see \code{calc_EG},
#' @param M: genotype likelihoods AA, Aa, aa, matrix sum(rdv) by 3 (double);
#' #' uses output from \code{calc_pobs_ndepth} function for simulation data and output from \code{getgenexp} or \code{getMAF} for VCF input
#' @param p: genotype frequencies AA, Aa, aa (double); output from \code{calc_EM} function.
#' @param rdv: read depth (vector of integers) for all samples
#' @return the variance of E(G_ij|D_ij)
#'
#' ##, var(X)=E(X^2)-E(X)^2, E(X) is calculated in calc_EG
#'
calc_EG_Var = function(M,p,rdv){
  LL = length(rdv)
  EG2 = NULL
  g = c(0,1,2)
  for (i in 1:LL){
    m = NULL
    for  (j in 1:3){
      L=M[i,j]
      m = c(m,L*p[j])
    }

    pm = sum(m/sum(m)*g^2) - sum(m/sum(m)*g)^2
    EG2 = c(EG2,pm)
  }
  return (EG2)
}

#' Andriy's original cal_EG2
#'
#' @describeIn calc_EG_Var also, see \code{calc_EG_Var},
#' @param M: genotype likelihoods AA, Aa, aa, matrix sum(rdv) by 3 (double);
#' #' uses output from \code{calc_pobs_ndepth} function for simulation data and output from \code{getgenexp} or \code{getMAF} for VCF input
#' @param p: genotype frequencies AA, Aa, aa (double); output from \code{calc_EM} function.
#' @param rdv: read depth (vector of integers) for all samples
#' @return the variance of E(G_ij|D_ij)
calc_EG_Var_general = function(M,p,rdv){
  L = length(rdv)
  S = 0
  EG = NULL
  g = c(0,1,2)
  for (i in 1:L){
    m = NULL
    for  (j in 1:3){
      L = 1
      for (kk in 1:rdv[i]){

        L = L*M[S + kk,j]
      }
      m = c(m,L*p[j])
    }
    S = S + rdv[i]
    pm = sum(m/sum(m)*g^2) - sum(m/sum(m)*g)^2
    EG = c(EG,pm)
  }
  return (EG)
}

#' EM calculation: Given a n by 3 matrix (M) containing likelihoods, this function uses the EM algorithm to estimate genotype frequencies.
#'
#' This part is the new part of the analysis where we can use reads instead of hard calls.
#' For this to work we need to estimate genotype frequencies.
#' Data consists of matrix n by 3
#' p = P(G=0)
#' q = P(G=1)
#'
#' @param M matrix of likelihoods derived from PL data.
#' @return a vector of three values: p0, q0, 1-p0-q0
calc_EM <- function(M){
  p_0 = 0.15
  q_0 = 0.15
  q_n = 1
  p_n = 0
  d_n = 0
  k = 0
  while ((p_n - p_0)^2 + (q_n - q_0)^2>0.000001){
    d_0 = 1-p_0 - q_0
    v = c(p_0,q_0,d_0)
    p_D = M%*%(v)
    E_p = M[,1]*p_0/p_D
    E_q = M[,2]*q_0/p_D
    p_n= p_0
    q_n = q_0
    d_n = 1-q_0-p_0
    p_0 = sum(E_p)/length(E_p)
    q_0 = sum(E_q)/length(E_q)
    k = k+1
    if (k==1000){
 #     cat('hi','\n')
      return ( c(p_0, q_0, 1-p_0-q_0) )
    }
  }
  return ( c(p_0, q_0, 1-p_0-q_0) )
}




