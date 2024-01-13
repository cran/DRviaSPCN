## ----include = FALSE----------------------------------------------------------
library(DRviaSPCN)

## ----eval=FALSE---------------------------------------------------------------
#  ### Download DRviaSPCNData package from GitHub
#  
#  library(devtools)
#  install_github("hanjunwei-lab/DRviaSPCNData",force = TRUE)
#  library(DRviaSPCNData)
#  ### Get weighted-ES of subpathways
#  DrugSPESCMatrix<-GetData("DrugSPESCMatrix")
#  ### Get p-value of subpathways centrality score
#  DrugSPPvalueMatrix<-GetData("DrugSPPvalueMatrix")

## ----message=FALSE------------------------------------------------------------
###Load depend package
library(igraph)
###Obtain input data
GEP<-GetExample('GEP')# Get the gene expression profile
Slabel<-GetExample('Slabel')# Get the sample class label

## ----eval=FALSE---------------------------------------------------------------
#  ###Run the function
#  CentralityScoreResult<-CalCentralityScore(ExpData=GEP,Label=Slabel,nperm = 1000)
#  
#  

## ----echo=FALSE---------------------------------------------------------------
###Get the result of this function
CentralityScoreResult<-GetExample('CentralityScoreResult')

## -----------------------------------------------------------------------------
###view first ten subpathways result
CentralityScoreResult[1:10,c(1,3,5,6,7)]



## ----eval=FALSE---------------------------------------------------------------
#  
#  ###Run the function
#  Opdrugresult<-Optimaldrugs(ExpData=GEP,Label=Slabel,DrugSPESC=DrugSPESCMatrix,
#                CentralityScore=CentralityScoreResult,nperm=1000,topcut=10,
#                pcut=0.01,weight=FALSE)
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


## ----echo=FALSE---------------------------------------------------------------
###Get the result of this function
sdf<-GetExample('methotrexate')

## ----eval=FALSE---------------------------------------------------------------
#  #Get the drug sdf data from DRviaSPCNData package
#  sdf<-DRviaSPCNData::GetData("sdfSET")

## ----results='hide',message=FALSE,fig.width=7,fig.height=5--------------------
getMolecularFm(drugname ="methotrexate",sdfSET=sdf)


## ----message=FALSE,results='hide',fig.width=7,fig.height=5--------------------
###Load depend package
library(GSVA)
library(pheatmap)

###Run the function
Disease2SPheatmap(CentralityScore=CentralityScoreResult,ExpData=GEP,Label=Slabel,pcut=0.05,
                   bk=c(-2,2),cluster.rows=FALSE,cluster.cols=FALSE,
                   show.rownames=TRUE,show.colnames=FALSE,
                   col=c("navy","firebrick3"),cell.width=NA,
                   cell.height=NA,scale="row",fontsize=7,
                   fontsize.row=9,fontsize.col=10)

## ----eval=FALSE---------------------------------------------------------------
#  ###Load depend package
#  library(GSVA)
#  library(pheatmap)
#  ###Run the function
#  Drug2SPheatmap(drugname="methotrexate_HL60_6_8.8e-06",
#                DrugSPPvalue=DrugSPPvalueMatrix,ExpData=GEP,
#                Label=Slabel,pcut=0.05,bk=c(-2,2),cluster.rows=FALSE,
#                cluster.cols=FALSE,show.rownames=TRUE,
#                show.colnames=FALSE,col=c("navy","firebrick3"),
#                cell.width=NA,cell.height=NA,scale="row",
#                fontsize=6,fontsize.row=9,fontsize.col=10)
#  
#  

## ----echo=FALSE---------------------------------------------------------------
knitr::include_graphics("../inst/heat.png")

