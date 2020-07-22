set.seed(1515)
D <- matrix(unlist(lapply(seq_len(3), function(m) rnorm(100))), nrow=100)
colnames(D) <- paste0("D", seq_len(3))
C1 <- matrix(unlist(lapply(seq_len(3), function(m) rnorm(100))), nrow=100)
colnames(C1) <- paste0("C1_", seq_len(3))
C2 <- matrix(unlist(lapply(seq_len(2), function(m) rnorm(100))), nrow=100)
colnames(C2) <- paste0("C2_", seq_len(2))
C3 <- matrix(unlist(lapply(seq_len(1), function(m) rnorm(100))), nrow=100)
colnames(C3) <- paste0("C3_", seq_len(1))
C <- list(C1, C2, C3)
Y <- rnorm(100)

dt <- data.frame(Y, D, C1, C2, C3)

dt_bal <- bal(c("mvGPS", "PS", "entropy"), dt[, colnames(D)], C)
W <- dt_bal$W
test_that("Argument check", {
    #improper formula specified, character value
    expect_error(mvGPS.fit("Y~D1+D2+D3", data=dt), 
                 "formula improperly specified. See details")
    #improper formula specified, single sided
    expect_error(mvGPS.fit(~D1+D2+D3, data=dt), "formula must be 2-sided")
    
    #improper W length
    W_miss <- lapply(W, function(x) x[-1])
    expect_error(mvGPS.fit(formula=Y~D1+D2+D3, data=dt, W=W_miss),
                 "all weights in W must be same length as number of units in exposure and confounders")
})

test_that("Performance check", {
    rslt <- mvGPS.fit(formula=Y~D1+D2+D3, W=W, data=dt)
    expect_equal(length(rslt), length(W) + 1)
    expect_named(rslt, c(names(W), "unweighted"))
})