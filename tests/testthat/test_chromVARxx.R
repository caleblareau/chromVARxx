
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
