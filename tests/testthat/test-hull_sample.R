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
})

test_that("Performance check", {
    #matrix of dimension 2
    mat2 <- hull_sample(D[, 1:2])
    expect_named(mat2, c("hpts_vs", "grid_pts"))
    
    #data.frame of dimension 2
    df2 <- hull_sample(dt[, 1:2])
    expect_named(df2, c("hpts_vs", "grid_pts"))
    
})