set.seed(1515)
D <- matrix(unlist(lapply(seq_len(3), function(m) rnorm(100))), nrow=100)

D_char <- D
D_char[, 2] <- ifelse(D[, 2]>0, "Pos", "Neg")


dt <- data.frame(D)
names(dt) <- paste0("X", seq_len(ncol(dt)))

dt_fact <- dt
dt_fact$X2 <- factor(ifelse(dt$X2>0, "Pos", "Neg"))
test_that("Argument check", {
    #matrix with at least one character variable
    expect_error(hull_sample(D_char), "X must be numeric")
    #data.frame with a factor variable
    expect_error(hull_sample(dt_fact), "X must be numeric")
    
    #specifying only a single exposure value with matrix
    expect_error(hull_sample(D[, 1]), "Exposure is not multivariate")
    #specifying only a single exposure value with data.frame
    expect_error(hull_sample(dt[, 1]), "Exposure is not multivariate")
    
    #specifying trim_hull but no quantile
    expect_error(hull_sample(D, trim_hull=TRUE), "trim_hull set to TRUE but trim_quantile not specified.")
    #improperly specifying trim_quantile value
    expect_error(hull_sample(D, trim_hull=TRUE, trim_quantile=0.01), "trim_quantile must be between [0.5, 1]", fixed=TRUE)
})

test_that("Performance check", {
    #matrix of dimension 2
    mat2 <- hull_sample(D[, 1:2])
    expect_named(mat2, c("hpts_vs", "grid_pts", "X"))
    
    #data.frame of dimension 2
    df2 <- hull_sample(dt[, 1:2])
    expect_named(df2, c("hpts_vs", "grid_pts", "X"))
    
    #matrix of dimension 3
    mat3 <- hull_sample(D)
    expect_named(mat3, c("hpts_vs", "grid_pts", "X"))
    expect_null(mat3$grid_pts)
    
    #data.frame of dimension 3
    df3 <- hull_sample(dt)
    expect_named(df3, c("hpts_vs", "grid_pts", "X"))
    expect_null(df3$grid_pts)
    
    #proper trimming check
    df3_trim <- hull_sample(dt, trim_hull=TRUE, trim_quantile=0.99)
    expect_named(df3_trim, c("hpts_vs", "grid_pts", "X"))
    expect_null(df3_trim$grid_pts)
    expect_true(all(apply(dt, 2, quantile, 0.99) >= apply(df3_trim$X, 2, max)))
    expect_true(all(apply(dt, 2, quantile, 1- 0.99) <= apply(df3_trim$X, 2, min)))
})