## ---- include = FALSE---------------------------------------------------------
library(DRviaSPCN)

## ----eval=FALSE---------------------------------------------------------------
#  ### Download DRviaSPCNData package from GitHub.
#  library(devtools)
#  install_github("hanjunwei-lab/DRviaSPCNData",force = TRUE)
#  library(DRviaSPCNData)
#  ### Get weighted-ES of subpathways.
#  DrugPscoreMatrix<-Getlist("DrugPscoreMatrix")
#  ## Get pvalue of subpathways centrality score.
#  DrugPvalueMatrix<-Getlist("DrugPvalueMatrix")

## ----message=FALSE------------------------------------------------------------
###Load depend package
library(igraph)
###Obtain input data
GEP<-GetExample('GEP')# Get the gene expression profile.
label<-GetExample('label')# Get the sample class label.
SubPathwayInfo<-GetExample('SubPathwayInfo')# Get the subpathway data
GoInfo<-GetExample('GoInfo')# Get the biology process data
Jaccardscore<-GetExample('Jaccardscore')# Get the jaccardscore matrix
GoSubPconGene<-GetExample('GoSubPconGene')# Get shared genes matrix

## ----eval=FALSE---------------------------------------------------------------
#  ###Run the function
#  DE2SubPathresult<-DE2SubPath(inexpData=GEP,Label=label,
#               Subpathway=SubPathwayInfo,Go=GoInfo,Jaccard=Jaccardscore,
#               Go_SubPath_gene=GoSubPconGene,perm=FALSE)
#  
#  DE2SubPathresult_P<-DE2SubPath(inexpData=GEP,Label=label,
#               Subpathway=SubPathwayInfo,Go=GoInfo,Jaccard=Jaccardscore,
#               Go_SubPath_gene=GoSubPconGene,perm=TRUE)
#  

## ----echo=FALSE---------------------------------------------------------------
###Get the result of this function
DE2SubPathresult<-GetExample('DE2SubPathresult')
DE2SubPathresult_P<-GetExample('DE2SubPathresult_P')

## -----------------------------------------------------------------------------
###view first ten subpathways result without random permutations
DE2SubPathresult[1:10,c(1,3,5,6,7)]
###view first ten subpathways result with random permutations
DE2SubPathresult_P[1:10,c(1,3,5,6,7)]


## ----message=FALSE------------------------------------------------------------

###Load depend package
library(clusterProfiler)

###Run the function
SubPathscore<-getSubpathscore(DE2SubPathresult=DE2SubPathresult,
                              inexpData=GEP,Label=label)

###view first ten subpathways result
head(SubPathscore,10)

## ----eval=FALSE---------------------------------------------------------------
#  
#  ###Run the function
#  Opdrugresult<-optimaldrugs(SubPathscore=SubPathscore,
#                 Drug_Pscore_matrix=DrugPscoreMatrix,nperm=1000,cut='p',
#                 topcut=20,pcut=0.01,weight=FALSE)
#  

## ----echo=FALSE---------------------------------------------------------------
###Get the result of this function
Opdrugresult<-GetExample('Opdrugresult')

## -----------------------------------------------------------------------------
###view first ten drugs result
head(Opdrugresult,10)

## ----message=FALSE,fig.width=7,fig.height=5-----------------------------------
###load depend package
library(igraph)
###plot network graph of the subpathway "00020_4"
plotSPW("00020_4")


## ----results='hide',message=FALSE,fig.width=7,fig.height=5--------------------
###Load depend package
library(ChemmineR)
library(rvest)
###Obtain molecular formula and visualize it.
Mole_formula<-getMolecularFm(drugname ="methotrexate")
plot(Mole_formula)


## ----message=FALSE,results='hide',fig.width=7,fig.height=5--------------------

###Load depend package
library(GSVA)
library(pheatmap)

###Run the function
Disease2SPWheatmap(DE2SubPathresult_P,exp=GEP,Label=label,pcut=0.05
                   ,bk=c(-2,2),cluster.rows=FALSE,cluster.cols=FALSE,
                   show.rownames=TRUE,show.colnames=FALSE,
                   col=c("navy","firebrick3"),cell.width=NA,
                   cell.height=NA,scale="row",fontsize=7,
                   fontsize.row=9,fontsize.col=10)

## ----eval=FALSE---------------------------------------------------------------
#  ###Load depend package
#  library(GSVA)
#  library(pheatmap)
#  ###Run the function
#  heatmap.list<-Drug2SPWheatmap(drugname="methotrexate",
#                Drug_Pvalue_matrix=DrugPvalueMatrix,exp=GEP,
#                Label=label,pcut=0.05,bk=c(-2,2),cluster.rows=FALSE,
#                cluster.cols=FALSE,show.rownames=TRUE,
#                show.colnames=FALSE,col=c("navy","firebrick3"),
#                cell.width=NA,cell.height=NA,scale="row",
#                fontsize=6,fontsize.row=9,fontsize.col=10)
#  
#  ###view the result
#  heatmap.list[[1]]
#  dev.off()

## ----echo=FALSE---------------------------------------------------------------
knitr::include_graphics("../inst/heat.png")

