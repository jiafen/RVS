
#' raw_reorder should be called if both case and controls are in the same dataset.
#'
#' It reorganize the cases' and controls' genotype, and return genotype with the first ncase rows are for case followed by ncont controls
#' @param path1: input for additive coding genotype files
#' @param path2: the file including the list of IIDs for cases
#' @param file3: the file name you want your processed genotype and phenotype to be named.
#' @return a list including- geno: the genotype data with the first ncase rows for cases, genotype data starting from 7th column, ncase: Num. of case
#' @export
raw_reorder<-function(path1,path2,file3)
{
raw<-read.table(path1,header=T,stringsAsFactors=F)
case<-read.table(path2,header=F,stringsAsFactors=F)

nsample=nrow(raw)
geno.case=subset(raw,IID%in%case$V1)
geno.cont=subset(raw,!IID%in%case$V1)
ncase<-nrow(case)
ncont<-nsample-ncase
Y=c(rep(1,ncase),rep(0,ncont))

nsnp=ncol(raw)-6
geno=rbind(geno.case,geno.cont)
geno.only=geno[,-(1:6)]
save(geno,geno.only,Y,ncase,ncont,file=file3)
return( list('geno'=geno,'Y'=Y,'ncase'=ncase))
}


#' subset_ind will subset the cases' or controls' genotype from a bigger dataset, and return the genotype for the wanted individuals.
#' @param path1: input for genotype files, ped file is prefered
#' @param path2: the file including the list of IIDs for individuals to be kept
#' @param file3: the file name you want your processed genotype and phenotype to be named.
#' @return a list including- geno: the genotype data with the first ncase rows for cases, genotype data starting from 7th column, ncase: Num. of case
subset_ind<-function(path1,path2)
{
  raw<-read.table(path1,header=F,stringsAsFactors=F)
  keep<-read.table(path2,header=F,stringsAsFactors=F)

  geno1=subset(raw,IID%in%keep$V1)
  n1<-nrow(geno1)

  #save(geno,geno.only,Y,ncase,ncont,file=file3)
  return( list('geno1'=geno1,'n1'=n1))
}


# Function #1
#' common_coord is to used to get intersect between sets
#' @param cord_set1: the coordinate of the snps in file1
#' @param cord_set2: the coordinate of the snps in file2
#' @return com_cord: a list include two vectors, one for the location of common snps in file1 and the other for file2.
#' @export
common_coord = function(cord_set1,cord_set2){
  set1 = NULL
  set2 = NULL
  L = length(cord_set1)
  for (i in 1:L){
    k = which(cord_set2 == cord_set1[i])
    if (length(k)==1){
      set1 = c(set1,i)
      set2 = c(set2,k)
    }
  }
  com_cord = list(cord1 = set1,cord2 = set2)
  return (com_cord)
}

# Function #2
#' subset_snp is to extract the genotype columns for chosen snps from a ped format file.
#'
#' This can also be done by software plink with option --extact snplist.txt, but snpnames are needed in snplist.txt to use plink.
#' @param M, the genotype file in ped format,Note, each snp have two columns in the ped file and the first 6 column need to be included.
#' @param set, the return of common_coord, the vector of integers representing the location of the common snps in the two sets.
#' @return Msub, the genotype file for the snps common in two sets only, each snp has two columns.
#' @export
subset_snp = function(M, set){ ## Note, each snp have two columns in the ped file

  set_k = 2*set + 5 ## 2*(set-1)+7, the start column Num of the genotype for the set-th snp
  set_k1 = set_k + 1  ## 2*set+6, the end column Num of the genotype for the set-th snp

  L = length(set_k)
  all_colns = NULL #1:6
  for (i in 1:L){
    all_colns  = c(all_colns,set_k[i],set_k1[i])
  }
  Msb = M[,all_colns]
  return (Msb)
}


#' snp_true_2T is used to change all TRUE genotype into T (if all the individual has genotype 'T' in one strand, R automatically change them into 'TRUE')
#'
#' @param M, the genotype matrix for one sample
#' @return M, the genotype matrix after change TRUE into T
#' @export
snp_true_2T <-function(M)
{
  kk=which(M[1,] == 'TRUE')
  for(i in kk){ M[,i]=factor('T')}
  return(M)
}

# Function #3
#'  ped_2_additive is to Combine two ped files for case and controls and transform them into additive coding
#'
#'   This function remove snps that are missing in one sample, monomorphic in whole sample, monomorphic in each sample with diff alleles, non biallelic
#'  alternatively, we can use plink merge the two files and generate the additive coding with option --recodeA
#'  plink will also check the strand +/- and the major/minor alleles.
#' @param M1, the genotype data for the first sample (normally case)
#' @param M2, the genotype data for the secondsample (normally controls)
#' @param ncols, the number of column before the genotype data start
#' @return A list including: Ad_all: the additive coding for all snps left, ncase: the first ncase rows are for cases;
#'  minor: the minor allele for each snp; rmed, the index of the removed snps from the input
#' @export
ped_2_additive = function(M1,M2,ncols){
  L = (length(M1[1,])-ncols)/2
  A1 = NULL
  A2 = NULL
  Minor=NULL
  rmed=NULL
  for (i in 1:L){
    c1=2*i-1+ncols
    c2=2*i+ncols
#    c1=1
#    c2=2
    samp1=M1[,c(c1,c2)]
    samp2=M2[,c(c1,c2)]
    l1 = levels(unlist(c(list(M1[,c1],M1[,c2]))))
    l2 = levels(unlist(c(list(M2[,c1],M2[,c2]))))

    set1=l1[l1 != '0']
    set2=l2[l2 != '0']
    d1 = length(set1)
    d2 = length(set2)

    if(d1 ==0 && d2 > 0)
    { cat('snp',i,' is missing in the first sample and removed! \n')
      rmed=c(rmed,i)
    }else if(d1>0 && d2 == 0){
      cat('snp',i,' is missing in the second sample and removed! \n')
      rmed=c(rmed,i)
    }else if(d1 == 1 && d2 == 1 && set1 == set2){
      cat('snp',i,' is monomorphic in the whole samples and removed! \n')
      rmed=c(rmed,i)
    }else if(d1 == 1 && d2 == 1 && set1 != set2){
      cat('snp',i,' is monomorphic in each sample but with different calls and removed! \n')
      rmed=c(rmed,i)
    }else if (d1 + d2>=3 && d1 + d2 <=4 && length(unique(c(set1,set2))) > 2)
    {cat('snp',i,' is biallelic in each sample but has more than 2 alleles in whole sample! removed! \n')
      rmed=c(rmed,i)
    }else if ( d1>2 || d2>2){cat('snp ',i,' has more than 2 alleles and removed! \n')
        rmed=c(rmed,i)
    }else{
     ma=minor_allele(M1[,c(c1,c2)],M2[,c(c1,c2)])
     A1=cbind(A1,convert(samp1,ma))
     A2=cbind(A2,convert(samp2,ma))
     Minor=c(as.character(Minor),as.character(ma))
    }
  }
    Mall=rbind(A1,A2)
    ncase=nrow(A1)
     addi = list(Ad_all=Mall,ncase=ncase,Minor=Minor,rmed=rmed)
     return (addi)
}

#' minor_allele is used to find the minor allele in the whole sample
#' @param M1: the 2 column of a snp from sample1
#' @param M2: the 2 column of a snp from sample2
#' @return minor: the minor allele from the sample1+sample2.
minor_allele <- function(M1, M2){
  c1 =unlist(c(list(M1[,1],M1[,2],M2[,1],M2[,2])))
  c1=c1[!c1%in%'0']
  t1=data.frame(table(c1))
  tt=which(t1[,1]%in%'0')
  if(length(tt>0)) {t1=t1[-tt,] }
  t1=t1[order(t1$Freq),]
  minor=(t1[1,1])
  return(minor)
}


#' convert two columns of a snp from ped file into additive coding
#'
#' @param M: 2 columns of a snp from ped file
#' @param ref, the first allele for the snp
#' @return Ad, a vector representing additive coding with length nsample for that snp.
convert = function(M,ref){
  L = length(M[,1])
  Ad = rep(0,L)
  for (i in 1:L){
    a  = M[i,]
    k = which(a=='0')
    if (length(k)>0){
      Ad[i] = NA
    }else{
      k = which(a==as.character(ref))
      Ad[i] = length(k)
    }
  }
  return (Ad)
}


#' filter_out filter out snps by missing rate 1-p
#'
#' It removes all the SNPs that missing rate are smaller than 1-p in either sample.
#' @param M1,M2 the additive coding of genotype for case and controls seperately
#' @param p, the keeping rate cut off, only keep snps with missing rate <(1-p) in both samples
#' @return list, M, including row combined genotype  and k-index for whether a snp pass the missingness cut off
#' @export
filter_out = function(M1,M2,p){
  L = length(M1[1,])
  J1 = length(M1[,1])
  J2 = length(M2[,1])
  k = rep(FALSE,L)
  for (i in 1:L){
    ll1 = sum(is.na(M1[,i]))
    ll2 = sum(is.na(M2[,i]))
    if ((ll1<(1-p)*J1) & (ll2<(1-p)*J2)){
        k[i] = TRUE
    }
  }
  tt=which(k>0)
#  cat('filter_keep=',tt)
  M = rbind(M1[,k],M2[,k])
  return (list(M=M, k =k))
}

#' get_maf is used to calculate the MAF
#'
#' @param M the additive coding of genotype for the whole sample
#' @return maf based on the non-missing data
#' @export
get_maf = function(M){
  L = length(M[1,])
  maf = rep(0,L)
  for (i in 1:L){
    LL = sum(!is.na(M[,i]))
    X = sum(M[!is.na(M[,i]),i])
    maf[i] = X/(2*LL)
  }
  return (maf)
}

#' get_minor is used to reset the maf if maf>0.5
#'
#' This is not needed if your genotype coding is based on the Num of minor allele
#' @param M the additive coding of genotype for the whole sample
#' @return maf based on the non-missing data
get_minor = function(M){
  L = length(M[1,])
  maf = get_maf(M)
  for (i in 1:L){
    if (maf[i]>0.5){
      a = abs(M[,i] - 2)
      M[,i] = a
    }
  }
  return (M)
}


#' combine_twogeo is used to combine two sets of ped pairs for two samples
#'
#' @param file1, the file name for the first sample before the extension. Eg. file1=sample1 if the file names are sample1.ped/map
#' @param file2,  the file name for the second sample before the extension. Eg. file1=sample2 if the file names are sample2.ped/map
#' @param keep1, if only a subset of sample1 is needed in the final analyis, keep1=1, ow, keep1=NA
#' @param keep2, if only a subset of sample2 is needed in the final analyis, keep2=1, ow, keep2=NA
#' @return a list include data sets 'common' and 'rare' for common  or rare variants and 'snp.Nomiss' for the SNPs kept in common and rare.
#' @export
combine_twogeno<-function(file1,file2,keep1=NA,keep2=NA)
{
  ## find the common SNPs in two datasets
  map1=read.table(paste0(file1,'.map'))
  map2=read.table(paste0(file2,'.map'))
  snpsub=common_coord(map1$V4,map2$V4)
  
  #### ped file to take the common variants in both datasets
  ped1=read.table(paste0(file1,'.ped'))
  ped2=read.table(paste0(file2,'.ped'))
  
  ped1.com=subset_snp(ped1,snpsub$cord1)
  ped2.com=subset_snp(ped2,snpsub$cord2)
  ## convert the subsetted ped file into additive coding
  ped1.com=snp_true_2T(ped1.com)
  ped2.com=snp_true_2T(ped2.com)
  
  pedall<-ped_2_additive(ped1.com,ped2.com,0)
  snp.all<-map1[snpsub$cord1,][-pedall$rmed,]
  
#  save(snp.all,file='genotype_allsample_ped2additive')
#  cat('keep1=',keep1,'keep2=',keep2,'\n')
  if(!is.na(keep1))
  {   tmp1<-read.table(paste0(file1,'_keep.txt'))
  index1<-which(ped1$V2%in%tmp1$V2)
  ncase <- length(index1)
  }
  if(!is.na(keep2))
  {   tmp2<-read.table(paste0(file2,'_keep.txt'))
  index2<-which(ped2$V2%in%tmp2$V2)
  }
  if(exists("index1") && exists("index2")) {
    ped.final <- pedall$Ad_all[c(index1,index2+pedall$ncase),]
  }else if (exists("index1") && (!exists("index2"))){
    rm=(1:pedall$ncase)[-index1]
    ped.final <- pedall$Ad_all[-rm,]
  }else{
    ped.final <- pedall$Ad_all[c(1:pedall$ncase,index2+pedall$ncase),]
  }
  if(!exists("ncase")) ncase<-pedall$ncase

  nall<-nrow(ped.final)
 # cat('ncase=',ncase,'nall=',nall,'\n')
  geno2test=filter_out(ped.final[1:ncase,],ped.final[(ncase+1):nall,],0.8)
  snp.Nomiss=snp.all[geno2test$k,]
  Y=c(rep(1,ncase),rep(0,nall-ncase))
  
  maf=get_maf(geno2test$M)
  common=geno2test$M[,maf>=0.05]
  rare=geno2test$M[,maf<0.05]
  snp.Nomiss=cbind(snp.Nomiss,maf)
  final=list('common'=common,'rare'=rare,'snp.Nomiss'=snp.Nomiss)
  cat('\nThere are', ncol(common), 'common SNPs and ',ncol(rare),' rare SNPs!\n')
  cat('All data are saved in binary data genotype_data, can be reached by final$common, final$rare and final$SNPs after loaded!\n') 
  save(final,file='genotype_data')
  return(final)
}  



