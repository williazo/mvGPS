#generating data for testing
set.seed(1515)
D <- matrix(unlist(lapply(seq_len(3), function(m) rnorm(100))), nrow=100)
C1 <- matrix(unlist(lapply(seq_len(3), function(m) rnorm(100))), nrow=100)
C2 <- matrix(unlist(lapply(seq_len(2), function(m) rnorm(100))), nrow=100)
C3 <- matrix(unlist(lapply(seq_len(1), function(m) rnorm(100))), nrow=100)
C <- list(C1, C2, C3)

dt <- data.frame(D, C1, C2, C3)
names(dt) <- paste0("X", seq_len(ncol(dt)))

test_that("Argument formatting", {
    expect_error(mvGPS(D[, 1], C), "Exposure must be multivariate. See details to ensure formula is properly specified")
    
    #improper length of arguments in C
    C_drop <- C
    C_drop[[1]] <- C_drop[[1]][-1, ]
    expect_error(mvGPS(D, C_drop), "Each matrix in C must have same number of observations, n, as D")
    
    #improper length of list in C
    C_drop <- C
    C_drop <- C_drop[-1]
    expect_error(mvGPS(D, C_drop), "Set of confounders not equal to number of exposures, m.")
    
    #common confounders specified without common argument
    expect_error(mvGPS(D, C1), "common=FALSE, C must be list of length m")
    
    #list of confounders specified with common argument
    expect_error(mvGPS(D, C, common=TRUE), "common=TRUE, expecting C to be single matrix of common confounders")
})

test_that("Expected output", {
    out <- mvGPS(D, C)
    expect_equal(length(out), 2)
    expect_named(out, c("score", "wts"))
    expect_equal(unname(lapply(out, length)), rep(list(100), 2))
    
    # checking that you can specify data.frames as well as matrices
    out_dt <- mvGPS(dt[, 1:3], C=list(dt[, 4:6], dt[, 7:8], dt[, 9]), common=FALSE)
    expect_equal(out$wts, out_dt$wts)
    
    out <- mvGPS(D, C1, common=TRUE)
    expect_equal(length(out), 2)
    expect_named(out, c("score", "wts"))
    
    #checking to make sure the trimming is occurring properly
    out_trim <- mvGPS(D, C1, common=TRUE, trim_w=TRUE)
    expect_true(quantile(out$wts, 0.99)==max(out_trim$wts))
    expect_true(quantile(out$wts, 1-0.99)==min(out_trim$wts))
})