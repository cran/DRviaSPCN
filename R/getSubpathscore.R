#'@title Calculating weighted enrichmentscore of subpathways
#'@description Function "getSubpathscore" used to calculate enrichmentscore weighted by centrality score of subpathways.
#'@param DE2SubPathresult A dataframe with seven columns those are subpath ID, subpath name, subpath size, genes in subpath, centralscore (eigenvector centrality), Pvalue and FDR (The result of function "DE2SubPath").
#'@param inexpData A gene expression profile of interest (rows are genes, columns are samples).The data in the expression profile is best not be log2 converted.
#'@param Label A character vector consist of "0" and "1" which represent sample class in gene expression profile. "0" means normal sample and "1" means disease sample.
#'@return A dataframe with three columns which are "SubPathID","Weighted-ES","Pvalue".
#'@importFrom clusterProfiler GSEA
#'@usage getSubpathscore(DE2SubPathresult,inexpData,Label)
#'@export
#'@examples
#'#Load depend package
#'library(clusterProfiler)
#'#Get input data (The "DE2SubPathresult" is the result of function "DE2SubPath")
#'DE2SubPathresult<-GetExample('DE2SubPathresult')
#'GEP<-GetExample('GEP')
#'label<-GetExample('label')
#'#Run the function
#'\donttest{SubPathscore<-getSubpathscore(DE2SubPathresult=DE2SubPathresult,inexpData=GEP,Label=label)}


getSubpathscore<-function(DE2SubPathresult,inexpData,Label){
  haveclusterProfiler <- PackageLoaded("clusterProfiler")
  if (haveclusterProfiler == FALSE) {
    stop("The 'clusterProfiler' library, should be loaded first")
  }
  DE<-getDEscore(inexpData,Label)
  pathset<-data.frame(path='1',gene='1')
  for (i in 1:length(DE2SubPathresult[,'SubPathID'])) {
    n<-DE2SubPathresult[,'SubPathID'][i]
    g<-unlist(strsplit(DE2SubPathresult[i,'Gene'],','))
    ps<-data.frame(path=rep(n,length(g)),gene=g)
    pathset<-rbind(pathset,ps)
  }
  pathset<-pathset[-1,]
  pathset[,1]<-as.character(pathset[,1])
  pathset[,2]<-as.character(pathset[,2])
  diseasegene<-DE[,1]
  names(diseasegene)<-rownames(DE)
  diseasegene<-sort(diseasegene,decreasing = T)


  diseaseES<-GSEA(diseasegene, TERM2GENE =pathset,exponent=1,minGSSize=0,
                  pAdjustMethod = "fdr", pvalueCutoff = 1)

  diseaseES<-diseaseES@result[,c(1,4,6)]

  diseasecentral<-as.numeric(DE2SubPathresult[,'Centralscore'])
  names(diseasecentral)<-as.character(DE2SubPathresult[,'SubPathID'])

  diseasecentral<-diseasecentral[match(rownames(diseaseES),names(diseasecentral))]
  a<-(sum(diseasecentral^2))^0.5
  b<-1+(diseasecentral/a)
  pathscore<-b*diseaseES[,2]
  result<-cbind(diseaseES[,1],pathscore)
  result1<-cbind(result,diseaseES[,3])
  result1<-as.data.frame(result1)
  result1[,1]<-as.character(result1[,1])
  result1[,2]<-as.numeric(as.character(result1[,2]))
  result1[,3]<-as.numeric(as.character(result1[,3]))
  colnames(result1)<-c('SubPathID','Weighted-ES','Pvalue')
  return(result1)
}
