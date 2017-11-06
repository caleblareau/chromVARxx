#' bagDeviations
#'
#' Given a deviations object from the cisBP
#' curtailed chromVARmotifs list, "bag" them
#' according to correlation where the most variable
#' goes first. Make sure to specify organism and
#' correlation. The lower the the correlation, 
#' the less bagging that occurs. Row data then
#' contains a comma-separated list of motifs that are
#' bagged.
#' 
#' @param object chromVARDeviations object from
#' chromVAR after 
#' @param cor = 0.7 minimum correlation between PWMS
#' @param organism = "human"; only other option is "mouse"
#' @import chromVAR
#' @import SummarizedExperiment
#' @export
#' @author Caleb Lareau
#' @examples
#'
#' rdsA<-paste0(system.file('rds',package='chromVARxx'),'/dev_humanSimple.rds')
#' object <- readRDS(rdsA)
#' bagged <- bagDeviations(object, cor = 0.3, organism = "human")
#' 
setGeneric("bagDeviations",
   function(object, cor, organism) standardGeneric("bagDeviations"))

#' @describeIn bagDeviations 
#' @export
setMethod("bagDeviations", c(object = "chromVARDeviations", cor = "numeric", organism = "character"),
          function(object, cor = 0.7, organism = "human"){
              
    stopifnot(organism %in% c("human", "mouse"))
    stopifnot(cor < 1 & cor > 0.15)
    
    # Get variability + TFs
    vb <- computeVariability(object, bootstrap_error = FALSE)[,"variability"]
    tfs <- rownames(object)
    names(vb) <- tfs
        
    vb <- vb[order(vb, decreasing=TRUE)]
    TFnames <- names(vb)
    
    # Import correlation based on PWMS for the organism
    rdsA <- paste0(system.file('cisBP_correlation',package='chromVARxx'),'/', 
                 organism, '_cisBP_correlationMatrix.rds')
    
    cormat <- readRDS(rdsA)
    
    stopifnot(all(TFnames %in% levels(cormat$TF1)))
    
    i <- 1
    TFgroups <- list()
    while(length(TFnames) != 0){

        # Pull out correlated hit
        tfcur <- TFnames[1]
        boo <- (cormat$TF1 == tfcur | cormat$TF2 == tfcur) & cormat$Pearson >= cor
        hits <- cormat[boo, ]
        tfhits <- unique(c(as.character(hits$TF1), as.character(hits$TF2), tfcur))
        
        #Update lists
        TFnames <- TFnames[!(TFnames %in% tfhits)]
        listHits <- list(tfhits)
        newname <- c(names(TFgroups), paste0("Group", as.character(i), "--", tfcur))
        TFgroups[i] <- listHits
        cormat <- cormat[cormat$TF1 %in% TFnames & cormat$TF2 %in% TFnames,]
        i <- i + 1
    }
    sentinalTFs <- sapply(TFgroups, function(x) x[1])
    objReduced <- object[sentinalTFs,]
    rowData(objReduced)$groups <- TFgroups
    return(objReduced)
})