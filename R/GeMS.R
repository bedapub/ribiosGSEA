#' @note TODO GeMS examples are suspended because the host cannot be reached outside of Roche network
NULL

#' GeMS base URL
#' To set GeMS base URL in your environment, use `GeMS_BASE_URL=value` in your
#' "~/.Renviron" file
#' @export
GeMS_BASE_URL <- Sys.getenv("GeMS_BASE_URL")
if(GeMS_BASE_URL=="") GeMS_BASE_URL <- "https://pred-gems.roche.com/api"

#' GeMS insert URL
#' @export
GeMS_INSERT_URL <- paste(GeMS_BASE_URL, "/insert", sep="")

#' GeMS remove URL
#' @export
GeMS_REMOVE_URL <- paste(GeMS_BASE_URL, "/remove", sep="")

#' GeMS genesets retrieval URL
#' @export
GeMS_GENESETS_URL <- paste(GeMS_BASE_URL, "/genesets", sep="")

## test URLs, GeMS API running as pod container on the biocomp server

#' GeMS URL for testing
#' @export
GeMS_TEST_URL <- "http://rkalvbiocomp.kau.roche.com:1234/api"

#' GeMS geneset retrieval URL for testing
#' @export
GeMS_TEST_GENESETS_URL <- paste(GeMS_TEST_URL, "/genesets", sep="")


#' Test whether GeMS is reachable
#' @return Logical value
#' @examples 
#' \dontrun{
#'   ## isGeMSReachble()
#' }
#' @export
isGeMSReachable <- function() {
   !httr::http_error(httr::GET(GeMS_GENESETS_URL)) 
}

#' Send a list as JSON query to an URL and fetch the response
#' 
#' @param url The destination URL
#' @param body A list to be sent to the URL, which will be encoded in the JSON format internally
#' 
#' @return The response from the webserver
#' 
#' @examples 
#' \dontrun{
#'    ## getJsonResponse(GeMS_GENESETS_URL, list(user=ribiosUtils::whoami()))
#' }
#' @export
getJsonResponse <- function(url, body) {
  response <- httr::POST(url, body=body, encode='json')
  ## TODO: in case the returned information is not a JSON, throw an error
  returnJSON <- jsonlite::fromJSON(httr::content(response, 'text', encoding='UTF-8'))
  return(returnJSON$response)
}

#' Construct message body to insert into GeMS
#' @param gmtList A \code{GmtList} object defined in the \code{BioQC} package
#' @param geneFormat Integer index of gene format. 0 stands for official human gene symbol
#' @param source Character, source of the gene set
#' @param taxID Integer, NCBI taxonomy ID of the species.
#' @param user The user name
#' @param subtype Subtype of the geneset
#' @param domain Domain of the geneset
#' 
#' @return A list with three items: \code{headers}, \code{parsed}, and \code{params}
#' 
#' @examples 
#' testList <- list(list(name="GS_A", desc=NULL, genes=c("MAPK14", "JAK1", "EGFR")),
#'   list(name="GS_B", desc="gene set B", genes=c("ABCA1", "DDR1", "DDR2")),
#'   list(name="GS_C", desc="gene set C", genes=NULL))
#' testGmt <- BioQC::GmtList(testList)
#' insertGmtListToGeMSBody(testGmt, geneFormat=0, source="Test")
#' @export
insertGmtListToGeMSBody <- function(gmtList,
                                    geneFormat=0,
                                    source="PubMed",
                                    taxID=9606,
                                    user=ribiosUtils::whoami(),
                                    subtype="",
                                    domain="") {
  parsed <- lapply(gmtList, function(x) {
    name <- x["name"]
    desc <- x["desc"]
    type <- ""
    genes <- x["genes"]
    res <- unname(unlist(c(name, desc, type, genes)))
  })
  names(parsed) <- NULL
  gemsPars <- list(gf=geneFormat,
                   so=source,
                   ti=taxID,
                   us=user)
  dataList <- list(headers=c("setName", "desc", "type", "genes"),
                   parsed=parsed,
                   params=gemsPars)
  return(dataList)
}

#' Insert a GmtList object to GeMS
#' @param gmtList A \code{GmtList} object defined in the \code{BioQC} package
#' @param geneFormat Integer index of gene format. 0 stands for official human gene symbol
#' @param source Character, source of the gene set
#' @param taxID Integer, NCBI taxonomy ID of the species.
#' @param user The user name
#' @param subtype Subtype of the geneset
#' @param domain Domain of the geneset
#' 
#' @return Response code or error message returned by the GeMS API. A value of \code{200} indicates a successful insertion.
#' 
#' @seealso \code{\link{removeFromGeMS}}
#' @examples 
#' \dontrun{
#'   testList <- list(list(name="GS_A", desc=NULL, genes=c("MAPK14", "JAK1", "EGFR")),
#'     list(name="GS_B", desc="gene set B", genes=c("ABCA1", "DDR1", "DDR2")),
#'     list(name="GS_C", desc="gene set C", genes=NULL))
#'   testGmt <- BioQC::GmtList(testList)
#'   ## insertGmtListToGeMS(testGmt, geneFormat=0, source="Test")
#'   ## removeFromGeMS(setName=c("GS_A", "GS_B", "GS_C"), source="Test")
#' }
#' @export
insertGmtListToGeMS <- function(gmtList,
                                geneFormat=0,
                                source="PubMed",
                                taxID=9606,
                                user=ribiosUtils::whoami(),
                                subtype="",
                                domain="") {
  dataList <- insertGmtListToGeMSBody(gmtList=gmtList,
                                  geneFormat=geneFormat,
                                  source=source,
                                  taxID=taxID,
                                  user=user,
                                  subtype=subtype,
                                  domain=domain)

  return(getJsonResponse(GeMS_INSERT_URL, dataList))
}

#' Message body to remove one or gene sets of the same source and user from GeMS
#' 
#' @param setName A vector of character strings, defining set names to be renamed. They must all have the same \code{source}, \code{user}, and \code{subtype}
#' @param source Character string, source of the gene set(s)
#' @param user Character string, user name
#' @param subtype Character string, subtype of the gene set(s)
#' 
#' @return A list of genesets to be removed, to be sent as message body
#' 
#' @examples 
#' removeFromGeMSBody(setName=c("GS_A", "GS_B", "GS_C"), source="Test")
#' @export
removeFromGeMSBody <- function(setName="",
                               source="", 
                               user=ribiosUtils::whoami(), 
                               subtype="") {
  toRemove <- lapply(setName, function(sname) {
    list(setName=sname,
         source=source,
         user=user,
         subtype=subtype)
  })
  res <- list(genesets=toRemove)
  return(res)
}

#' Rmove one or gene sets of the same source and user from GeMS
#' @param setName A vector of character strings, defining set names to be renamed. They must all have the same \code{source}, \code{user}, and \code{subtype}
#' @param source Character string, source of the gene set(s)
#' @param user Character string, user name
#' @param subtype Character string, subtype of the gene set(s)
#' 
#' @return Response code or error message returned by the GeMS API. A value of \code{200} indicates a successful insertion.
#' 
#' @seealso \code{\link{insertGmtListToGeMS}}
#' @examples 
#' \dontrun{
#'   testList <- list(list(name="GS_A", desc=NULL, genes=c("MAPK14", "JAK1", "EGFR")),
#'     list(name="GS_B", desc="gene set B", genes=c("ABCA1", "DDR1", "DDR2")),
#'     list(name="GS_C", desc="gene set C", genes=NULL))
#'   testGmt <- BioQC::GmtList(testList)
#'   ## insertGmtListToGeMS(testGmt, geneFormat=0, source="Test")
#'   ## removeFromGeMS(setName=c("GS_A", "GS_B", "GS_C"), source="Test")
#' }
#' @export
removeFromGeMS <- function(setName="", 
                           source="", 
                           user=ribiosUtils::whoami(), 
                           subtype="") {
  body <- removeFromGeMSBody(setName, source, user, subtype)
  return(getJsonResponse(GeMS_REMOVE_URL, body))
}

#' Get gene sets of a user from GeMS
#' @param user User name
#' @return A data.frame including following columns: 
#' \enumerate{
#'   \item setName
#'   \item desc
#'   \item domain
#'   \item source
#'   \item subtype
#' }
#' 
#' @examples 
#' \dontrun{
#' #### my gene-sets
#' ## getUserSetsFromGeMS()
#' #### from another user
#' ## getUserSetsFromGeMS("kanga6")
#' }
#' @export
getUserSetsFromGeMS <- function(user=ribiosUtils::whoami()) {
  fieldsOfInterest <- c("setName", "desc", "domain",
                         "source", "subtype")
  body <- list(user=user,
               returnParams=as.list(fieldsOfInterest))
  df <- getJsonResponse(GeMS_GENESETS_URL, body)
  if(is.list(df) && length(df)==0) {
    res <- NULL
  } else {
    res <- df[, fieldsOfInterest]
  }
  return(res)
}

#' Get one gene-set with its name
#' @param setName Character string
#' @return A list of two elements
#' \enumerate{
#'   \item name
#'   \item genes
#' }
#' @seealso \code{\link{getSetsWithNamesFromGeMS}}
#' @examples 
#' \dontrun{
#' getSetWithNameFromGeMS("Plasma_sc")
#' }
#' @export
getSetWithNameFromGeMS <- function(setName) {
  body <- list(setName=setName)
  response <- getJsonResponse(GeMS_GENESETS_URL, body)
  my_content <- unlist(lapply(response$genes[[1]], function(x) x[[1]][1]))
  res <- list(name=setName, 
              desc=response$desc[1],
              genes=my_content)
  return(res)
}

#' Get one or more gene-sets with their names
#' @param setNames Character strings
#' @return A GmtList object
#' @seealso \code{\link{getSetWithNameFromGeMS}}
#' @examples 
#' \dontrun{
#' getSetsWithNamesFromGeMS(c("Plasma_sc", "Bcell_l_Danaher17"))
#' }
#' @export
getSetsWithNamesFromGeMS <- function(setNames=NULL) {
  gmtList <- lapply(setNames, getSetWithNameFromGeMS)
  res <- BioQC::GmtList(gmtList)
  return(res)
}

#' Get gene-sets for application
#' @param property Character string, property to query
#' @param value Character string, property value
#' @return A GmtList object
#' @examples 
#' \dontrun{
#' getSetsWithPropertyFromGeMS("meta.application", "rtbeda_CIT")
#' }
#' @export
getSetsWithPropertyFromGeMS <- function(property="meta.application",
                                       value="") {
  body <- list(property=value)
  names(body) <- property
  response <- getJsonResponse(GeMS_GENESETS_URL, body)
  gmt_lists <- lapply(1:nrow(response), function(i) {
    list(name=response$setName[i],
         desc=response$desc[i],
         genes=sapply(response$genes[[i]], function(x) x[[1]][1]))
  })
  res <- GmtList(gmt_lists)
  return(res)
}

