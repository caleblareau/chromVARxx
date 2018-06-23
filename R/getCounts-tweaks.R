#' get counts for large single cell data
#'
#' If bulk, use the standard getCounts function in
#' chromVAR. Otherwise, if you have a barcoded
#' sample and a particular sam tag that distinguishes
#' cells, this function should be much faster. Assumes
#' that duplicate reads are already removed and other
#' artifacts (e.g. proper pairing) are also already
#' handled upstream. 
#' 
#' @param bamfile Valid string for filepath to bam file
#' @param peaks A GRanges object of accessibility peaks.
#' Use chromVAR::getPeaks to import a bed file for this.
#' @param barcodeTag The two-letter sam tag associated
#' with the barcoding of the different cells
#' @param mapqFilter Minimum mapping quality metric necessary
#' for read to be used in counts matrix.
#' 
#' @return A summarized Experiment roughly equivalent
#' to what chromVAR::getCoutns would provide
#' 
#' @import SummarizedExperiment
#' @import dplyr
#' @importFrom GenomicAlignments readGAlignments findOverlaps
#' @importFrom S4Vectors queryHits subjectHits 
#' @importFrom Rsamtools ScanBamParam
#' @export
#' @author Caleb Lareau
#' @examples
#'
#' bamfile <-paste0(system.file('raw',package='chromVARxx'),'/chr1.barcode.small.bam')
#' peakfile <- paste0(system.file('raw',package='chromVARxx'),'/chr1.peaks.small.bed')
#' peaks <- suppressWarnings(chromVAR::getPeaks(peakfile))
#' barcodeTag <- "CB"
#' SE <- getCountsByID(bamfile, peaks, barcodeTag, mapqFilter =0 )
#' 
setGeneric("getCountsByID",
           function(bamfile, peaks, barcodeTag, mapqFilter) standardGeneric("getCountsByID"))

#' @describeIn getCountsByID 
#' @export
setMethod("getCountsByID", c(bamfile = "character", peaks = "GenomicRanges", barcodeTag = "character", mapqFilter = "numeric"),
          function(bamfile, peaks, barcodeTag, mapqFilter = 0){
              
              # Import alignments and get overlaps with peaks
              GA <- readGAlignments(bamfile, use.names = TRUE, param = ScanBamParam(
                  tag = c(barcodeTag), mapqFilter = 0))
              ovPEAK <- findOverlaps(peaks, GA)
              
              # Determine universe of unique barcodes
              barcodes <- as.character(mcols(GA)[,barcodeTag])
              uniqueBarcodes <- unique(barcodes)
              id <- factor(barcodes, levels = uniqueBarcodes)
              
              # Assemble counts
              countdf <- data.frame(peaks = queryHits(ovPEAK),
                                    sample = as.numeric(id)[subjectHits(ovPEAK)],
                                    read = names(GA)[subjectHits(ovPEAK)]) %>%
                  distinct() %>%  # by filtering on distinct read / peak / sample trios, ensure that PE reads that overlap peak aren't double counted
                  select(-one_of("read")) %>% 
                  group_by(peaks,sample) %>% summarise(count = n()) %>% data.matrix()
              
              m <- Matrix::sparseMatrix(i = c(countdf[,1], length(peaks)),
                                        j = c(countdf[,2], length(uniqueBarcodes)),
                                        x = c(countdf[,3],0))
              colnames(m) <- uniqueBarcodes
              
              # Generate colData
              depth <- data.frame(
                  sample = as.numeric(id),
                  read = names(GA)) %>%
                  distinct() %>% group_by(sample) %>% summarise(depth = n()) %>% data.matrix()
              colData <- data.frame(
                sample = uniqueBarcodes,
                depth = depth[,2],
                FRIP = Matrix::colSums(m)/depth[,2]
              )
              
              # Make sure that the SE can be correctly constructed
              stopifnot(all(colData$sample == colnames(m)))
              
              # Make summarized Experiment
              SE <- SummarizedExperiment(
                  rowRanges = peaks, 
                  assays = list(counts = m),
                  colData = colData
              )
              return(SE)
          })
