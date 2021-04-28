set.seed(1515)
D <- matrix(unlist(lapply(seq_len(3), function(m) rnorm(80))), nrow=80)
C1 <- matrix(unlist(lapply(seq_len(3), function(m) rnorm(80))), nrow=80)
C2 <- matrix(unlist(lapply(seq_len(2), function(m) rnorm(80))), nrow=80)
C3 <- matrix(unlist(lapply(seq_len(1), function(m) rnorm(80))), nrow=80)
C <- list(C1, C2, C3)

dt <- data.frame(D, C1, C2, C3)
names(dt) <- paste0("X", seq_len(ncol(dt)))
test_that("Argument check", {
    #unused model arguments warning
    expect_warning(bal(model_list=c("mvGPS", "kernel", "super"), D, C))
    
    #checking fuzzy matching for arguments
    out <- bal(model_list=c("mv", "ent", "PS"), D=D, C=C)
    expect_equal(out$models, c("mvGPS", "entropy", "PS"))
    expect_named(out, c("W", "cor_list", "bal_metrics", "ess", "models"))
    
    out_trim <- bal(model_list=c("mv", "ent", "PS"), D=D, C=C, trim_w=TRUE)
    
    #checking the min and max trimming
    expect_equal(unname(unlist(lapply(out$W, quantile, 0.99))),
                 unname(unlist(lapply(out_trim$W, max, 0.99))))

    expect_equal(unname(unlist(lapply(out$W, quantile, 1-0.99))),
                 unname(unlist(lapply(out_trim$W, min, 0.99))))
    
    #for each univariate method we expect them to have metric for each exposure,
    #i.e, 2 * 3 and we have one from mvGPS and one unweighted
    expect_equal(nrow(out$bal_metrics), 1 + 1 + 2 * 3)
    
    #if all_uni=FALSE then we should have only one weight for each method
    out_uni <- bal(model_list=c("mv", "ent", "PS"), D=D, C=C, all_uni=FALSE)
    expect_equal(nrow(out_uni$bal_metrics), 1 + 1 + 2 * 1)
    
    out_uni_trim <- bal(model_list=c("mv", "ent", "PS"), D=D, C=C, 
                        all_uni=FALSE, trim_w=TRUE)
    #checking the min and max trimming
    expect_equal(unname(unlist(lapply(out_uni$W, quantile, 0.99))),
                 unname(unlist(lapply(out_uni_trim$W, max))))
    
    expect_equal(unname(unlist(lapply(out_uni$W, quantile, 1-0.99))),
                 unname(unlist(lapply(out_uni_trim$W, min))))
    
    #there is an additional warning whenever GBM is used
    expect_warning(bal(model_list=c("GBM"), D=D, C=C))
    
    #checking passing options to weightit functions. here we expect the GBM function to not have an error
    expect_warning(bal(model_list=c("GBM"), D=D, C=C, stop.method="p.mean"), regexp=NA)
    
})
