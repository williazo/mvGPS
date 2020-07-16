set.seed(1515)
D <- matrix(unlist(lapply(seq_len(3), function(m) rnorm(100))), nrow=100)
C1 <- matrix(unlist(lapply(seq_len(3), function(m) rnorm(100))), nrow=100)
C2 <- matrix(unlist(lapply(seq_len(2), function(m) rnorm(100))), nrow=100)
C3 <- matrix(unlist(lapply(seq_len(1), function(m) rnorm(100))), nrow=100)
C <- list(C1, C2, C3)
test_that("Argument check", {
    expect_error(bal(model_list="mvGPS", D[, 1], C), 
                 "Exposure must be multivariate. See details to ensure formula is properly specified")
    
    #unused model arguments
    expect_warning(bal(model_list=c("mvGPS", "kernel", "super"), D, C))
})