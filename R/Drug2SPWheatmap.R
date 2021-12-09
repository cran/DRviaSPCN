##' @title Plot a heat map of the subpathways activity regulated by drugs
##' @description The "Drug2SPWheatmap" function plots a heat map of the subpathways that are regulated by a specific drug.
##' @param drugname A character which represent interest drug name.
##' @param Drug_Pvalue_matrix A matrix which colunms represent drugs and rows respresent subpathways. Values in this matrix is the pvalue of subpathways centrality score regulated by drugs.
##' @param exp A gene expression profile of interest.
##' @param Label A character vector consist of "0" and "1" which represent sample class in gene expression profile. "0" means normal sample and "1" means disease sample.
##' @param pcut A numeric value which represent threshold. Subpathways with p-value less than this threshold will be screened out and visualized.
##' @param bk A numeric vector that covers the range of values. Users could adjust color depth through this parameter.
##' @param cluster.rows Boolean values determining if rows should be clustered or hclust object.
##' @param cluster.cols Boolean values determining if columns should be clustered or hclust object.
##' @param show.rownames Boolean specifying if row names are be shown.
##' @param show.colnames Boolean specifying if column names are be shown.
##' @param col Vector of colors used in heatmap.
##' @param cell.width Individual cell width in points. If left as NA, then the values depend on the size of plotting window.
##' @param cell.height Individual cell height in points. If left as NA, then the values depend on the size of plotting window.
##' @param scale Character indicating if the values should be centered and scaled in either the row direction or the column direction, or none. Corresponding values are "row", "column" and "none".
##' @param fontsize Base fontsize for the plot (default: 10).
##' @param fontsize.row Fontsize for rownames (default: 10).
##' @param fontsize.col Fontsize for colnames (default: 10).
##' @return A list where the elements are heat maps of drugs acting on different cell lines,duration and concentrations.
##' @importFrom GSVA gsva
##' @importFrom pheatmap pheatmap
##' @importFrom grDevices colorRampPalette
##' @usage Drug2SPWheatmap(drugname="",Drug_Pvalue_matrix,exp,Label,pcut=0.05,
##'          bk=c(-2,2),cluster.rows=FALSE,cluster.cols=FALSE,show.rownames=TRUE,
##'          show.colnames=FALSE,col=c("navy","firebrick3"),
##'          cell.width=NA,cell.height=NA,scale="row",fontsize=7,
##'          fontsize.row=9,fontsize.col=10)
##' @export
##' @examples
##' #Load depend package
##' library(GSVA)
##' library(pheatmap)
##' ##Obtain input data
##' #Statistic significance of subpathways centrality score
##' #this function need were stored in packet "DRviaSPCNData".
##' #"DRviaSPCNData" has been uploaded to the github repository.
##' #Users can download and install through "install_github" function and
##' #set parameter url="hanjunwei-lab/DRviaSPCNData".
##' #After installing and loading package "DRviaSPCNData",
##' #users can use the following command to get the data:
##' #DrugPvalueMatrix<-Getlist('DrugPvalueMatrix')
##' GEP<-GetExample('GEP')
##' label<-GetExample('label')
##'
##' #Run the function
##' \donttest{heatmap.list<-Drug2SPWheatmap(drugname = "methotrexate",
##'                  Drug_Pvalue_matrix=DrugPvalueMatrix,
##'                  exp=GEP,Label=label,pcut=0.05,bk=c(-2,2),
##'                  cluster.rows=FALSE,cluster.cols=FALSE,show.rownames=TRUE,
##'                  show.colnames=FALSE,col=c("navy","firebrick3"),
##'                  cell.width=NA,cell.height=NA,scale="row",fontsize=7,
##'                  fontsize.row=9,fontsize.col=10)
##' #view the result
##' heatmap.list[[1]]
##' dev.off()
##' heatmap.list[[2]]
##' dev.off()
##' heatmap.list[[3]]
##' dev.off()}



Drug2SPWheatmap<-function(drugname="",Drug_Pvalue_matrix,exp,Label,pcut=0.05,bk=c(-2,2),
                             cluster.rows=FALSE,cluster.cols=FALSE,show.rownames=TRUE,
                             show.colnames=FALSE,col=c("navy","firebrick3"),
                             cell.width=NA,cell.height=NA,scale="row",fontsize=7,
                             fontsize.row=9,fontsize.col=10){
  haveGSVA <- PackageLoaded("GSVA")
  havepheatmap <- PackageLoaded("pheatmap")
  if (haveGSVA == FALSE) {
    stop("The 'GSVA' library, should be loaded first")
  }
  if (havepheatmap == FALSE) {
    stop("The 'pheatmap' library, should be loaded first")
  }

  sp<-GetExample('SubPathwayInfo')

  d.name<-c()
  for (i in 1:length(colnames(Drug_Pvalue_matrix))) {
    z<-unlist(strsplit(colnames(Drug_Pvalue_matrix)[i],split = '_'))[1]
    d.name<-c(d.name,z)
  }

  Drug_Pscore_matrix1<-Drug_Pvalue_matrix[,which(d.name==drugname)]
  pscore.list<-list()
  for (i in 1:length(Drug_Pscore_matrix1[1,])) {
    pscore.list[[i]]<-Drug_Pscore_matrix1[,i]
  }
  names(pscore.list)<-colnames(Drug_Pscore_matrix1)

  spm.list<-lapply(pscore.list,function(p){

    p<-names(which(p<pcut))

    sp<-sp[sp[,1]%in%p,]

    pl<-list()
    for (i in 1:length(sp[,1])) {
      pl[[i]]<-strsplit(sp[i,4],',')[[1]]
    }

    names(pl)<-sp[,1]

    #library(GSVA)
    spw_matrix = gsva(as.matrix(exp), pl, method = "ssgsea",
                      kcdf = "Gaussian", abs.ranking = TRUE,min.sz=2)
    return(spw_matrix)
  })

  heat.list<-list()
  for (i in 1:length(spm.list)) {

    spw_matrix<-spm.list[[i]]
    colnames(spw_matrix)[which(Label=='0')]='normal'
    colnames(spw_matrix)[which(Label=='1')]='disease'
    spw_matrix<-spw_matrix[,order(colnames(spw_matrix))]
    colann=data.frame(Sample=factor(rep(names(table(colnames(spw_matrix))),
                                        table(colnames(spw_matrix)))))
    samples<-paste(colnames(spw_matrix),1:length(spw_matrix[1,]))
    colnames(spw_matrix)<-samples
    rownames(colann)<-samples
    colnames(colann)<-'Group'

    bk1<-c(seq(bk[1],-0.1,by=0.1),seq(0,bk[2],by=0.1))
    #library(pheatmap)
    heat.map<-pheatmap(spw_matrix,
                       scale = scale,
                       breaks = bk1,
                       color = colorRampPalette(c(col[1], "white", col[2]))(50),
                       cluster_rows=cluster.rows,cluster_cols=cluster.cols,
                       annotation_col =colann,
                       show_rownames=show.rownames,show_colnames=show.colnames,
                       cellwidth = cell.width,cellheight = cell.height,
                       fontsize = fontsize,fontsize_row = fontsize.row,
                       fontsize_col = fontsize.col,
                       main = paste("Heatmap of the activities of significant subpathways regulated by the",
                                    names(spm.list)[i])

    )

    heat.list[[i]]<-heat.map
  }
  names(heat.list)<-names(spm.list)
  return(heat.list)
}
