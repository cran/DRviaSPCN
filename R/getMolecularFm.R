#' @title Plot chemical molecular formula  of drugs
#' @description The function "getMolecularFm" outputs the chemical molecular formula  of a drug or compound . The results can be visualized by the "plot" function.
#' @param drugid A character string of DrugBank ID.
#' @param drugname A character string of drug name.
#' @param sdfSET Sdf data of drug structure.
#' @param main An overall title for the chemical structure graph.
#' @param sub A sub title for the chemical structure graph.
#' @return Chemical molecular formula  of the drug or compound.
#' @importFrom sp plot
#' @importFrom ChemmineR sdf2ap
#' @usage getMolecularFm(drugid = NULL, drugname = NULL,sdfSET,main = "", sub = "")
#' @export
#' @examples
#'#"sdfSET" has been uploaded to the
#'#github repository.Users can download and install through "install_github"
#'#function and set parameter url="hanjunwei-lab/DRviaSPCNData".
#'#After installing and loading package "DRviaSPCNData",
#'#users can use the following command to get the data.
#' # Obtain molecular formula and visualize it.
#' \donttest{
#' #Get the sdf data of drug structure from DRviaSPCNData package
#' #library("Chemminer")
#' #sdf<-GetData('sdfSET')
#' #Run the function
#' #Mole_formula<-getMolecularFm(drugname ="methotrexate",sdfSET=sdf)
#' #plot(Mole_formula)}



getMolecularFm<-function(drugid = NULL, drugname = NULL,sdfSET,main = "", sub = "") {
  dn<-GetExample('dn')

  if(is.null(drugid)==TRUE){
    drugname <- unlist(strsplit(drugname, "\\("))[1]
    drugname1 <- tolower(drugname)
    drugid<-dn[which(dn$Drug.name==drugname1),1]
  }
  if(length(drugid)==0){
    print("Sorry, the drug ID or name you provided is not in our records.")
  }

  sdfset<-sdfSET[which(sdfSET@ID==drugid)]
  if(length(sdfset)==0){
    stop("Sorry, the drug ID or name you provided is not in our records.")
  }
  if (main == "") {
    sdfset@ID <- drugid
  }
  else {
    sdfset@ID <- main
  }
  ap <- sdf2ap(sdfset[[1]])
  sp::plot(sdfset)

}




