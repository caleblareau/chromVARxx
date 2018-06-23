
context("bagging deviations works")

test_that("Bagging deviations reduces the feature space", {
    
    # Human
    rdsA<-paste0(system.file('rds',package='chromVARxx'),'/dev_humanSimple.rds')
    object <- readRDS(rdsA)
    bagged <- bagDeviations(object, cor = 0.3, organism = "human")
    expect_equal(dim(object)[1], 1764)
    expect_equal(dim(bagged)[1], 8)
    
    # Mouse
    rdsB<-paste0(system.file('rds',package='chromVARxx'),'/dev_mouseSimple.rds')
    object <- readRDS(rdsB)
    baggedm <- bagDeviations(object, cor = 0.3, organism = "mouse")
    expect_equal(dim(object)[1], 1346)
    expect_equal(dim(baggedm)[1], 5)
})

context("New import getCounts functions do what we think that they should")

test_that("getCountsByID creates a proper SummarizedExperiment object", {
    
    bamfile <-paste0(system.file('raw',package='chromVARxx'),'/chr1.barcode.small.bam')
    peaks <- chromVAR::getPeaks(paste0(system.file('raw',package='chromVARxx'),'/chr1.peaks.small.bed'))
    barcodeTag <- "CB"
    SE <- getCountsByID(bamfile, peaks, barcodeTag, mapqFilter=0)
    expect_equal(dim(SE)[1], 18427)
    expect_equal(dim(SE)[2], 2748)
})
