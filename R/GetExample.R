#' @title Get example data
#' @description This function is used to achieve exxample data.
#' @param exampleData  A character, should be one of"GEP","label","SubPathwayInfo","GoInfo","GoSubPconGene",
#' "Jaccardscore","SubPathwaymapdata","Drugs_CID","DE2SubPathresult","DE2SubPathresult_P",
#' "Opdrugresult" and "heatmap.list".
#' @return example data
#' @usage GetExample(exampleData)
#' @export


GetExample<-function(exampleData){
  if(!exists("envData")) {
    utils::data("envData",package="Disease2Drug")
  }

  if (exampleData=="GEP")
  {
    dataset<- get("GEP",envir=envData)
    return(dataset)
  }
  if (exampleData=="label")
  {
    dataset<- get("label",envir=envData)
    return(dataset)
  }
  if (exampleData=="SubPathwayInfo")
  {
    dataset<- get("SubPathwayInfo",envir=envData)
    return(dataset)
  }
  if (exampleData=="GoInfo")
  {
    dataset<- get("GoInfo",envir=envData)
    return(dataset)
  }
  if (exampleData=="Jaccardscore")
  {
    dataset<- get("Jaccardscore",envir=envData)
    return(dataset)
  }
  if (exampleData=="GoSubPconGene")
  {
    dataset<- get("GoSubPconGene",envir=envData)
    return(dataset)
  }
  if (exampleData=="SubPathwaymapdata")
  {
    dataset<- get("SubPathwaymapdata",envir=envData)
    return(dataset)
  }
  if (exampleData=="Drugs_CID")
  {
    dataset<- get("Drugs_CID",envir=envData)
    return(dataset)
  }
  if (exampleData=="DE2SubPathresult")
  {
    dataset<- get("DE2SubPathresult",envir=envData)
    return(dataset)
  }
  if (exampleData=="DE2SubPathresult_P")
  {
    dataset<- get("DE2SubPathresult_P",envir=envData)
    return(dataset)
  }
  if (exampleData=="Opdrugresult")
  {
    dataset<- get("Opdrugresult",envir=envData)
    return(dataset)
  }
  if (exampleData=="heatmap.list")
  {
    dataset<- get("heatmap.list",envir=envData)
    return(dataset)
  }
}
