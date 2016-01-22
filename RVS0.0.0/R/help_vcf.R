#' Seprates the case and control IDs.
#'
#' This function read the header of the vcf file (rows starting with # in vcf) and split the case and controls according to the caseID_file.
#' @param vcf_file, The 'address/filename' of the vcf file.
#' @param caseID_file, The 'address/filename' of the file including case IDs.
#' @param IDline, The number of the line in the vcf file whch contains the ID's, the provided vcf has IDline=128.
#' @param headnum, The number of columns in the headers before the IDs starts in the vcf file, default=9.
#' @return a list which contains the position (column) numbers which correspond to the cases and controls in the row of the VCF with IDs, but columns are shifted by the headnum columns.
#' @examples getIDs("../data/1g113low_1g56exomehigh_filtered.hg19.chr11_1000snps.vcf","../data/bam.txt",128,9)
#' @export
get_IDs <- function(vcf_file,caseID_file,IDline=128,headnum=9){
  filen<-vcf_file
  filecon = file(filen, open='r')
  on.exit(close(filecon))

  tt2 = readLines(filecon, n=IDline) ## read the first IDlines of the file

  S = unlist(strsplit(tt2[IDline],'\t')) #S is now a list of all the ID's including the names of columns before IDs start
  Sl = length(S)
  h1 = headnum+1
  S=S[h1:Sl] #Have to get rid of the first few components which are not the IDs

  CaseIDs <- read.table(caseID_file) #CASE ID's
  L = length(S)
  case <- NULL
  cont <- NULL
  #IF theres an identical ID like the Exome ID in the vcf, can call it a case, otherwise call it a control
  case=which(S%in%CaseIDs$V1)
  controls=1:L
  controls=controls[-case]
  output <- list("cases"=case,"controls"=controls)
  return(output)
}


#' Filters SNPs based on specific criteria:
#'
#' This function read the whole vcf file, filter out SNPs and return the genotype likelihood in VCF file for case and control separately.
#' it calls \link{\code{get_Likelihood_vcf}} to convert the vcf format into phred score.
#' @details the keeping criteria are: SNPs that PASS filters, do not have duplicates, have only one reference  allele and at most 50% missing rate.
#' @param vcf_file, The address/filename' of the vcf file.
#' @param case, The cases location in the whole samples, output from the getIDs function.
#' @param cont, The controls location in the whole samples, output from the getIDs function.
#' @param nsnps, total number of snps in the vcf file. default=1000
#' @param nread, each time the read.table will read nread number of rows in file, default=300
#' @param headnum, default =9, the number of columns before each individual's data
#' @param missing_th, the missing rate for filter out a snp, default is set=0.5
#' @return The function returns 6 matrices and 2 vectors:
#' @return 6 matrices:case00, case01, case11, cont00, cont01 and cont11 which contain the genotype likelihoods from VCF
#' @return 2 vectros: chr and position (Loc) for the snps kept for further analysis.
#'  @examples filterSNPs('1g113low_1g56exomehigh_filtered.hg19.chr11_1000snps.vcf',IDs$cases,IDs$controls,)
filter_SNPs<- function(vcf_file,case,cont,nsnps=1000, nread=300, headnum=9,missing_th=0.5){
  filen = vcf_file  # ../data/1g113low_1g56exomehigh_filtered.hg19.chr11_1000snps.vcf
  filecon = file(filen, open='r')
  on.exit(close(filecon))

  # Genotype likelihoods for cases
  # A0M -  major homozygous
  # A1M -  minor heterozygous
  # A2M -  minor homozygous
  A0M <- NULL
  A1M <- NULL
  A2M <- NULL
  #
  # Genotype likelihoods for controls
  # B0M -  major homozygous
  # B1M -  minor heterozygous
  # B2M -  minor homozygous
  B0M <- NULL
  B1M <- NULL
  B2M <- NULL
  Loc <- NULL  ## to store the location of returned snps
  Chr <- NULL

  mod <- nsnps %% nread
  if (mod ==0) { nloop = nsnps/nread
  s <- TRUE
  }else {
    nloop <- floor(nsnps/nread)+1
    s <- FALSE
  }

  loop <- 0
  while ( loop< nloop) {
    loop <- loop+1

    # Start reading vcf file
    # when filecon has tens of thousands rows, it is hard to read them all in one time.
    if(s == TRUE | loop < nloop) { FF = read.table(filecon, nrows = nread, sep='\t')
    }
    if (s == FALSE & loop ==nloop) {
      FF = read.table(filecon, nrows = mod, sep='\t')
    }
    if ( class(FF) == 'try-error'){break}
    if ( length(FF) == 0){break}

    chr<-as.character(FF$V1)
    loc<-as.character(FF$V2)
    ref<-as.character(FF$V4)
    alt<-as.character(FF$V5)
    filter<-as.character(FF$V7)


    #
    # select SNPs have no duplication (uniq_loc), pass the filter (pass_loc)
    # have unique reference allele or alternative allele (uniq_ref_loc and uniq_alt_loc)
    ### see the intersection afterwards
    #

    ## find the location of SNPs have no duplication
    dup=which(duplicated(loc))
    unique_loc=which(!loc%in%loc[dup])   ### location of the unique location

    ## find the location of SNPs that pass the criteria in GATK algorithm
    pass_loc=which(filter%in%'PASS')
    ## find the location of SNPs that only have one reference allele or alternative allele
    uniq_ref_loc=which(nchar(ref)==1)
    uniq_alt_loc=which(nchar(alt)==1)


    ## find the location of SNPs that meet all the above criteria
    t11=intersect(unique_loc,pass_loc)
    t12=intersect(t11,uniq_ref_loc)
    t13=intersect(t12,uniq_alt_loc)

    ## SNPs pass the above criteria
    Fpass1<-FF[t13,]

    #### get the likelihood for these first filtered SNPs
    AA=get_likelihood_vcf(Fpass1,headnum)  ## need to ask Andriy why he has AA, and BB and SameI and SameKK

    L1 = AA$SNPs[,'loc']
    chr= AA$SNPs[,'chr']

    A00=AA$L_matrix$L00[case,]
    A01=AA$L_matrix$L01[case,]
    A02=AA$L_matrix$L11[case,]

    B00=AA$L_matrix$L00[cont,]
    B01=AA$L_matrix$L01[cont,]
    B02=AA$L_matrix$L11[cont,]

    #### select variants that with missing rate smaller than the threshold given by TH for case and controls separately .

    TH = missing_th

    case1=apply(A00,2,is.na)  ## if A00 is missing for one row, A01, A02 are missing as well
    case2=apply(case1,2,sum)  ## count how many people have missing probability for each location
    mis_case=case2/nrow(A00) ## missing rate for each SNP in case

    cont1=apply(B00,2,is.na)
    cont2=apply(cont1,2,sum)
    mis_cont=cont2/nrow(B00)

    col_keep= (mis_case<TH) & (mis_cont<TH)

    A0M = cbind(A0M, A00[,col_keep])
    A1M = cbind(A1M, A01[,col_keep])
    A2M = cbind(A2M, A02[,col_keep])

    B0M = cbind(B0M, B00[,col_keep])
    B1M = cbind(B1M, B01[,col_keep])
    B2M = cbind(B2M, B02[,col_keep])

    Loc = c(Loc,L1[col_keep])  ## SNPs will be analyzed further
    Chr = c(Chr,chr[col_keep])

    cat('while loop:',loop,'\n')
    cat('case dimension: ', dim(A0M),'\n')
  }

  rm('A00','A01','A02','B00','B01','B02','AA','L1')
  snps<-data.frame(cbind(Chr,Loc))
  colnames(snps)=c('chr','loc')
  return(list("SNPs"=snps,"case00"=A0M,"case01"=A1M,"case11"=A2M,"cont00"=B0M,"cont01"=B1M,"cont11"=B2M))
}

#' Generates the expected probabilities of the genotypes E(G_ij|D_ij).
#'
#' @param A0M, Genotype lkelihood for major homozygous of the cases.  output of filterSNPs, extract from VCF files
#' @param A1M, Genotype lkelihood for minor heterozygous of the cases.
#' @param A2M, Genotype lkelihood for minor homozygous of the cases.
#' @param B0M, Genotype lkelihood for major homozygous of the controls.
#' @param B1M, Genotype lkelihood for minor heterozygous of the controls.
#' @param B2M, Genotype lkelihood for minor homozygous of the controls.
#' @param chr, chromosome of the snp, output from the  function filterSNPs
#' @param Loc, location of the snp, Output from the function filterSNPs.
#' @return list including: 2 matrix (MM and P) and 1 data frame SNPs
#' @return MM -  Matrix contains expected genotypes for each person at each loci, row - individual, column - SNP
#' @return  P -  Matrix for genotype frequency in all the sample at each loci, 3 columns contains P(G=0), P(G=1) and P(G=2), row - SNP
#' @return SNPs - a data frame contains chr and loc of the returned SNPs
#' @examples
#' geneexp <- getgenexp(SNPs$A0M,SNPs$A1M,SNPs$A2M,SNPs$B0M,SNPs$B1M,SNPs$B2M,SNPs$L11)
getgenexp <- function(A0M,A1M,A2M,B0M,B1M,B2M,chr,Loc){

  A0 = A0M
  A1 = A1M
  A2 = A2M

  B0 = B0M
  B1 = B1M
  B2 = B2M
  loc = Loc

  MM = NULL
  P = NULL

  #number of snps
  Lsnp = ncol(B0)

  for (i in 1:Lsnp){
    A = cbind(A0[,i],A1[,i],A2[,i])
    B = cbind(B0[,i],B1[,i],B2[,i])
    M = rbind(A,B)

    EG = rep(NA,nrow(M))
    kk0 = !is.na(M[,1])
    M = M[kk0,]
    ## length(M[,1]) is the number of sample, here read depth for each sample =1
    rdv = rep(1,nrow(M))
    p = calc_EM(M)
    #p = c(0.9^2,2*0.9*0.1,0.1^2)

    if ((sum(p< (-0.00001))>0) || (sum(p>1.00001)>0)){
      cat('Unreasonable genotype frequency by EM algorithm!!!')
      return (NULL)
    }
    t1<-which(p<(-0.00001)); p[t1]=0;
    t1<-which(p>1.00001); p[t1]=1;
    P = rbind(P,p)
    EG[kk0] = calc_EG(M,p,rdv)
    MM = cbind(MM,EG)
  }
  t1<-1:nrow(P)
  rownames(P)<-paste('snp',t1,sep='')
  colnames(P)=c('P00','P01','P11')

  colnames(MM)=rownames(P)
  t1<-1:nrow(MM)
  rownames(MM)<-paste('sample',t1,sep='')
  SNP<-data.frame(cbind(chr,loc))
  return(list("SNPs"=SNP,"pop_frq"=P,"exp_cond_prob"=MM))
}


#' Calculates the MAF
#'
#' @description
#' Calculates the MAF and returns matrices which are to be inputed into a test function of the rvs package.
#'
#' @param P Output from the getgenexp function which contains the probabilities of the genotypes.
#' @param MM Output matrix from the getgenexp function.
#' @param ncase, number of cases
#' @param method, need value of common or rare
#' @return The function returns matrix P which contains the probabilities of the genotypes and
#'  matrices M1 and M2, which can all be inputed into a test function of the rvs package.
#' @examples
#' getMAF(geneexp$P,geneexp$MM)
getMAF <- function(P,MM,chr,loc,ncase,method='common'){
  if(!method %in%c('common','com','Com','Common','Rare','rare')){
    cat('Wrong method!, need be common or rare only!\n');
    return(NULL)
  }
  maf = P%*%c(0,1,2)/2
  maf[maf>0.5] = 1-maf[maf>0.5]

  if(method %in% c('common','com','Com','Common'))
  {  kk=which(maf>0.05)
  }
  if(method %in%c('Rare', 'rare'))
  {
    kk=which(maf<0.05)
  }
  MM = MM[,kk]
  P = P[kk,]
  loc1=loc[kk]
  chr1=chr[kk]
  maf=maf[kk]
  M1 = MM[1:ncase,] # cases
  M2 = MM[(ncase+1):nrow(MM),] #controls
  #  SNPs<-cbind(chr,loc)
  #  head(SNPs)
  SNPs<-cbind(chr1,loc1)
  SNPs<-cbind(SNPs,P,maf)
  return(list("SNPs"=SNPs,"cases"=M1,"controls"=M2))
}


#' Gets Likelihoods for all possible genotypes given the set of alleles defined in the REF and ALT fields
#' This function fetches the genotype likelihoods from the VCF file (wo header) for the analysis and
#' must contain a 9th column (ie. 'FORMAT' field) consisting of (GT:AD:DP:GQ:PL).
#' The first step is to read the vcf file as a matrix and outputs an R object with 6 lists including 5 vectors and 1 list including 3 matrix:
#'  5 vectors: chr, position,ref allele, alt allele, filter indicator for 'PASS'
#' 1 list: likelihood of L(00|ref,alt), L(01|ref,alt), L(11|ref,alt), which are three matrix with row for sample and col for variants
#'
#' Input matrix:
#' 2nd column (POS)    - location
#' 4th column (REF)    - reference allele
#' 5th column (ALT)    - alternative allele
#' 7th column (FILTER) - vcr filter value
#' 10th column and onwards are 'data' columns in format (GT:AD:DP:GQ:PL)
#'
#' GL = Phred-scaled likelihoods for genotypes as defined in the VCF specification format of 1000 Genome project,
#' likelihoods for all possible genotypes given the set of alleles defined in the REF and ALT fields.
#' Phred-scaled likelihoods: -log10(L(00|ref,alt)), -log10(L(01|ref,alt)), -log10(L(11|ref,alt)).
#' PL : the phred-scaled genotype likelihoods rounded to the closest integer (and otherwise defined precisely as the GL field) (Integers)
#'
#' Example: 0/0:2,0:4:6:0,6,42
#'
#' @param M The vcf file without the first headers in a matrix (nsnps rows and nhead column +each column for one individuals).
#' @param nhead, the number of columns before the data for each individual
#' @return Outputs an R object with:\cr
#' 1. chromosome of variant
#' 2. POS of variant location (vector)\cr
#' 3. REF allele (vector)\cr
#' 4. ALT allele (vector)\cr
#' 5. FILTER indicator for 'PASS' (vector)\cr
#' 6. likelihoods: L(D|0), L(D|1), L(D|2). List of 3 matrix.
#'
#' L0 - matrix of P(D|0), rows are individuals and columns are SNPS\cr
#' L1 - matrix of P(D|1), rows are individuals and columns are SNPS\cr
#' L2 - matrix of P(D|2), rows are individuals and columns are SNPS\cr
get_likelihood_vcf <- function(M, nhead ){
  ref = as.vector(M[,4])
  alt = as.vector(M[,5])
  loc = as.vector(M[,2])
  filter = as.vector(M[,7])
  chr = as.vector(M[,1])

  L00 = NULL
  L01 = NULL
  L11 = NULL
  Lsnp = nrow(M)
  Lind = ncol(M)-nhead

  for (i in 1:Lsnp){
    l0=NULL
    l1=NULL
    l2=NULL
    for (j in 1:Lind){
      temp<-M[i,j+nhead]
      if(temp =='./.') {
        p_loci=c(NA,NA,NA)
      }else {
        p_loci=get_lsingle_vcf(temp)
      }
      l0=c(l0,p_loci[1])
      l1=c(l1,p_loci[2])
      l2=c(l2,p_loci[3])
    }


    L00 = cbind(L00,l0)
    L01 = cbind(L01,l1)
    L11 = cbind(L11,l2)
  }
  SNP=cbind(chr,loc,ref,alt,filter)
  return ( list('SNPs'=SNP, "L_matrix"=list( "L00"=L00, "L01"=L01, "L11"=L11) ) )
}


#' Parse Phred-scaled likelihoods
#'
#' This function parses the data column (GT:AD:DP:GQ:PL) to get the PL values \code{-10log10(x)},
#' and converts them to a vector of likelihoods: L0,L1,L2. It's called within \code{get_likelihood_vcf}.
#'
#' @param M:  a single data unit of format: GT:AD:DP:GQ:PL.  Example: get_l('0/0:2,0:4:6:0,6,42')
#' @keywords called by \code{get_likelihood_vcf}
#' @return a vector of three values: l0, l1, l2
#' @seealso \code{\link{get_likelihood_vcf}} for parent function.
get_lsingle_vcf <- function( Mij ){
  a = unlist(strsplit(as.vector(Mij),':'))
  ## the original code in Andiry:
  if(a[5]=='.') {return(c(NA,NA,NA)) }
  #if ((length(a)==1 & a =='./.' )|| (length(a)>1 & a[5] == '.')) { return (c(NA,NA,NA)) }
  l = as.numeric(unlist(strsplit(a[5],',')))
  l0 = l[1]
  l1 = l[2]
  l2 = l[3]
  Pl0 = 10^(-l0/10)
  Pl1 = 10^(-l1/10)
  Pl2 = 10^(-l2/10)
  return ( c(Pl0,Pl1,Pl2) )
}


#' the main function to call the VCF file and return the  expected probabilities of the genotypes E(G_ij|D_ij) from VCF file
#'
#' @param file, the VCF file
#' @param case_ID_file, the caseID files
#' @param IDline, the row number of your VCF which contains all sample ID
#' @param nhead, the column number before the sample ID
#' @param missing_th, the missing rate cut off for a SNP to be filtered, if missing>missing_th in case or controls, the SNP will be removed
#' @param nsnp, should be the number of SNPs in your VCF file, but if VCF is too large, you can choose only read in the first nsnp SNPs
#' @param index, if it is common,  SNPs with maf<MAF will be removed, if it is rare, SNPs with maf>MAF will be removed
#' @param MAF, the MAF to be used to filter out SNP,
#'
#' @return matrix of E(G_ij|D_ij)
#' @export
vcf_process<-function(file='example/1g113low_1g56exomehigh_filtered.hg19.chr11_1000snps.vcf',
                      caseID_file='example/caseID.txt',
                      IDline=128, nhead=9, missing_th=0.2, nsnp=1000, nread=300, MAF=0.05)
{
  ## seperate the case and controls, return their column locations
  IDs <- get_IDs(file,caseID_file,IDline,nhead);

  nsample=length(IDs$cases)+length(IDs$controls)
  cat('This VCF includes ',nsample,' samples with ', length(IDs$cases), 'cases!\n')
  cat('cases location in all samples:\n')
  cat(IDs$cases)
  cat('\n\n')

  ## filter out SNPs
  SNPs <- filter_SNPs(file, IDs$cases, IDs$controls, nsnps=nsnp, nread=nread, missing_th=missing_th);

  nsnp.left=nrow(SNPs$SNPs)
  cat(nsnp.left, 'SNPs out of ', nsnp, 'SNPs are kept!\n\n')

  ## getgenexp1 has second filter on missingness of samples, they have the same result if the missing_th2=missing in filter_SNPs.
  ##  exp_prob <- getgenexp1(SNPs$case00, SNPs$case01, SNPs$case11, SNPs$cont00, SNPs$cont01, SNPs$cont11, SNPs$chr, SNPs$Loc, missing_th2=0.5)
  exp_prob<-getgenexp(SNPs$case00, SNPs$case01, SNPs$case11, SNPs$cont00, SNPs$cont01, SNPs$cont11, SNPs$SNPs[,'chr'], SNPs$SNPs[,'loc'])

  common<-getMAF(exp_prob$pop_frq,exp_prob$exp_cond_prob,exp_prob$SNPs[,'chr'],exp_prob$SNPs[,'loc'],ncase=length(IDs$cases),method='common')
  rare<-getMAF(exp_prob$pop_frq,exp_prob$exp_cond_prob,exp_prob$SNPs[,'chr'],exp_prob$SNPs[,'loc'],ncase=length(IDs$cases),method='rare')

  save(common,file='Common_variants')
  save(rare,file='Rare_variants')

  ncommon=nrow(common$SNPs)
  nrare=nrow(rare$SNPs)
  cat('There are ', ncommon, 'common variants and ', nrare, 'rare variants!\n\n')
  cat('Common variants info and expected probabilities are saved in Rdata Common_variants!\n')
  cat('Rare variants info and expected probabilities are saved in Rdata Rare_variants!\n')
  return(list('common'=common,'rare'=rare))
}

