getKS<-function(s,l) {
  if(class(s)=='character'&class(l)=='character'){
    V <- match(s,l)
  }
  if(class(s)=='numeric'&class(l)=='character'){
    V <- match(names(s),l)
  }
  if(class(s)=='character'&class(l)=='numeric'){
    V <- match(s,names(l))
  }
  if(class(s)=='numeric'&class(l)=='numeric'){
    V <- match(names(s),names(l))
  }
  V <- V[!is.na(V)]
  V <- sort(V)
  n <- length(l)
  t <- length(V)
  j <- 1:t
  a <- j/t - V/n
  a <- max(a)
  b <- V/n - (j - 1)/t
  b <- max(b)
  if (a > b) {
    ks = a
  }
  else {
    ks = -b
  }
  return(ks)
}

CalculateSES<-function (labels.list, correl.vector = NULL) {
  tag.indicator <- labels.list
  no.tag.indicator <- 1 - tag.indicator
  N <- length(labels.list)
  Nh <- length(labels.list[labels.list == 1])
  Nm <- N - Nh
  correl.vector <- abs(correl.vector)
  sum.correl.tag <- sum(correl.vector[tag.indicator == 1])
  norm.tag <- 1/sum.correl.tag
  norm.no.tag <- 1/Nm
  RES <- cumsum(tag.indicator * correl.vector * norm.tag -
                  no.tag.indicator * norm.no.tag)
  max.ES <- max(RES)
  min.ES <- min(RES)
  ES <- signif(ifelse(max.ES > -min.ES, max.ES, min.ES), digits = 5)
  return(ES)
}


#'@title Identifying the optimal drugs
#'@description Function "optimaldrugs" used to identify the optimal drugs for specific disease.
#'@param SubPathscore A dataframe with three columns which are "SubPathID","Weighted-ES","Pvalue" (The result of function "getSubPathscore").
#'@param Drug_Pscore_matrix A matrix with n rows and m columns. n is the number of subpathways and m is the number of all drugs. The values in this matrix is weighted enrichmentscore of subpathways in every drug. The users could obtain this matrix from our example data.
#'@param nperm Number of random permutations (default: 1000).
#'@param cut There are two ways to select up-regulated and down-regulated subpathways. The up-regulated subpathways (down-regulated subpathways) is the top (bottom) subpathways of list in descending order of weighted enrichmentscore (weighted-ES) when cut="top". When cut="p",up-regulated subpathways and down-regulated subpathways is screened based on ES and pvalue in the results of GSEA.
#'@param topcut When cut="top", topcut represents the number of selected up-regulated subpathways or down-regulated subpathways.The topcut defaults to 20.
#'@param pcut When cut="p", pcut represents the threshold of statistical significance level for screen subpathways. The pcut defaults to 0.05.
#'@param weight A boolean value determines the method for calculating the drug-disease association score of the drug. "weight=FALSE"(default): Similar to "CMap" (Lamb et al., 2006), no weight is needed. "weight=TRUE": KS random walk statistic with individualized subpathway activity score as weight was used to calculate the drug-disease reverse association score.
#'@return A dataframe with four columns which are "Drug"(drug names),"KS"(final drug-disease assciation score),"pvalue"(statistical significance),"FDR"(statistical significance after adjust).
#'@importFrom stats p.adjust
#'@usage optimaldrugs(SubPathscore,Drug_Pscore_matrix,nperm=1000,cut='top',
#'                    topcut=20,pcut=0.05,weight=FALSE)
#'@export
#'@examples
#'##Obtain input data
#'#Weighted enrichmentscore of subpathways this function need were stored
#'#in packet "DRviaSPCNData". "DRviaSPCNData" has been uploaded to the
#'#github repository.Users can download and install through "install_github"
#'#function and set parameter url="hanjunwei-lab/DRviaSPCNData".
#'#After installing and loading package "DRviaSPCNData",
#'#users can use the following command to get the data.
#'#DrugPscoreMatrix<-Getlist('DrugPscoreMatrix')
#'DE2SubPathresult<-GetExample("DE2SubPathresult")
#'GEP<-GetExample("GEP")
#'label<-GetExample("label")
#'\donttest{
#'SubPathscore<-getSubpathscore(DE2SubPathresult=DE2SubPathresult,inexpData=GEP,Label=label)
#'#Run the function
#'Opdrugresult<-optimaldrugs(SubPathscore=SubPathscore,Drug_Pscore_matrix=DrugPscoreMatrix,
#'                           nperm=1000,cut='p',topcut=20,pcut=0.01,weight=FALSE)}


optimaldrugs<-function(SubPathscore,Drug_Pscore_matrix,nperm=1000,cut='top',topcut=20,
                       pcut=0.05,weight=FALSE){
  havestats <- PackageLoaded("stats")
  if (havestats == FALSE) {
    stop("The 'stats' library, should be loaded first")
  }
  if(cut=='top'){
    SubPathscore<-SubPathscore[order(SubPathscore[,2],decreasing = TRUE),]
    uppath<-SubPathscore[,1][1:topcut]
    downpath<-SubPathscore[,1][(length(SubPathscore[,1])-(topcut-1)):length(SubPathscore[,1])]
  }
  if(cut=='p'){
    uppath<-SubPathscore[SubPathscore[,2]>0&SubPathscore[,3]<pcut,1]
    downpath<-SubPathscore[SubPathscore[,2]<0&SubPathscore[,3]<pcut,1]
  }
  if(weight==FALSE){
    ks_up <- apply(Drug_Pscore_matrix, 2,function(y) {
      names(y) <- rownames(Drug_Pscore_matrix)
      y <- sort(y,decreasing = TRUE)
      V <- match(uppath, names(y))
      V <- V[!is.na(V)]
      V <- sort(V)
      n <- length(y)
      t <- length(V)
      j <- 1:t
      a <- j/t - V/n
      a <- max(a)
      b <- V/n - (j - 1)/t
      b <- max(b)
      if (a > b) {
        ks_up = a
      }
      else {
        ks_up = -b
      }
      return(ks_up)
    })
    ks_down <- apply(Drug_Pscore_matrix, 2,function(y) {
      names(y) <- rownames(Drug_Pscore_matrix)
      y <- sort(y,decreasing = TRUE)
      V <- match(downpath, names(y))
      V <- V[!is.na(V)]
      V <- sort(V)
      n <- length(y)
      t <- length(V)
      j <- 1:t
      a <- j/t - V/n
      a <- max(a)
      b <- V/n - (j - 1)/t
      b <- max(b)
      if (a > b) {
        ks_down = a
      }
      else {
        ks_down = -b
      }
      return(ks_down)
    })
  }
  if(weight==TRUE){
    ks_up <- apply(Drug_Pscore_matrix, 2,function(y) {
      names(y) <- rownames(Drug_Pscore_matrix)
      y <- sort(y,decreasing = T)
      P.rank<-names(y)
      tag.indicator <- sign(match(P.rank,
                                  uppath, nomatch = 0))
      up.ES <- CalculateSES(tag.indicator, correl.vector = y)
      return(up.ES)
    })
    ks_down <- apply(Drug_Pscore_matrix, 2,function(y) {
      names(y) <- rownames(Drug_Pscore_matrix)
      y <- sort(y,decreasing = TRUE)
      P.rank<-names(y)
      tag.indicator <- sign(match(P.rank,
                                  downpath, nomatch = 0))
      down.ES <- CalculateSES(tag.indicator, correl.vector = y)
      return(down.ES)
    })
  }
  KS<-c()
  for (i in 1:length(ks_up)) {
    if(ks_up[i]*ks_down[i]>0){
      k<-0
    }else{
      k<-ks_up[i]-ks_down[i]
    }
    KS<-c(KS,k)
  }
  KS<-data.frame(KS,ks_up)

  #names(KS)<-colnames(Drug_Pscore_matrix)
  ##标准化
  MA<-max(KS[,1])
  MI<-min(KS[,1])
  for(i in 1:length(KS[,1])){

    if(KS[,1][i]>0){
      KS[,1][i]<-KS[,1][i]/MA
    }
    if(KS[,1][i]<0){
      KS[,1][i]<--(KS[,1][i]/MI)
    }
  }
  KS<-KS[order(KS[,1],KS[,2],decreasing = TRUE),]
  KSorder<-KS[,1]
  names(KSorder)<-rownames(KS)
  ###按药名提取集合，富集到KS上，计算初始ks值
  drugname<-c()
  for (i in 1:length(KSorder)) {
    d<-unlist(strsplit(names(KSorder[i]),'_'))[1]
    drugname<-c(drugname,d)

  }
  uniquename<-unique(drugname)

  druglist<-list()
  for (i in 1:length(uniquename)) {
    druglist[[i]]<-KSorder[which(drugname==uniquename[i])]
    names(druglist)[i]<-uniquename[i]
  }

  ks_list <- lapply(druglist, getKS,KSorder)##计算初始ks值

  ##扰动
  permlist<-list()
  for (i in 1:length(druglist)) {
    permmatrix<-matrix(NA,length(druglist[[i]]),nperm)
    for (j in 1:nperm) {
      permmatrix[,j]<-sample(names(KSorder),length(druglist[[i]]),replace = FALSE)
    }
    permlist[[i]]<-permmatrix
  }

  pml<-lapply(permlist,function(x){
    apply(x, 2,getKS,KSorder)
  })##得到每个药扰动后的ks

  ##算p值
  pvalue<-c()
  for(i in 1:length(ks_list)){
    if(ks_list[[i]]>0){
      p<-length(which(pml[[i]]>ks_list[[i]]))/nperm
    }

    if(ks_list[[i]]<0){
      p<-length(which(pml[[i]]<ks_list[[i]]))/nperm
    }
    p<-length(which(pml[[i]]<ks_list[[i]]))/nperm
    pvalue<-c(pvalue,p)
  }
  fdr<-p.adjust(pvalue,'fdr',length(pvalue))
  result<-data.frame(matrix(unlist(ks_list),
                            nrow = length(ks_list),byrow = TRUE),
                     stringsAsFactors = FALSE)

  result<-cbind(names(ks_list),result)
  result<-cbind(result,pvalue)
  result<-cbind(result,fdr)
  result<-result[order(result[,3],decreasing = F),]
  colnames(result)<-c('Drug','KS','pvalue','FDR')
  return(result)
}
