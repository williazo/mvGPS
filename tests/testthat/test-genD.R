test_that("Argument check", {
    expect_error(gen_D("o", 200, 0.8, 1, 1, 3, rep(0, 3), 1, rep(0, 3), rep(0, 3)))
    
})
test_that("Performance check", {
    #checking to make sure seed is properly set
    u1<- gen_D(method="u", n=200, rho_cond=0.8, s_d1_cond=1, s_d2_cond=1, k=3, 
               C_mu=rep(0, 3), C_cov=0.2, C_var=1, 
               d1_beta=rep(0, 3), d2_beta=rep(0, 3), seed=1)
    u2 <- gen_D(method="u", n=200, rho_cond=0.8, s_d1_cond=1, s_d2_cond=1, k=3, 
                C_mu=rep(0, 3), C_cov=0.2, C_var=1, 
                d1_beta=rep(0, 3), d2_beta=rep(0, 3), seed=1)
    expect_equal(u1$D, u2$D)
    
    #checking that marginal correlation is properly calculated
    zero_rho <-gen_D(method="u", n=200, rho_cond=0, s_d1_cond=1, s_d2_cond=1, k=3, 
                     C_mu=rep(0, 3), C_cov=0, C_var=1, 
                     d1_beta=rep(0, 3), d2_beta=rep(0, 3), seed=1)
    expect_equal(zero_rho$rho, 0)
    
    high_rho <-gen_D(method="u", n=200, rho_cond=0.8, s_d1_cond=1, s_d2_cond=1, k=3, 
                     C_mu=rep(0, 3), C_cov=0, C_var=1, 
                     d1_beta=rep(0, 3), d2_beta=rep(0, 3), seed=1)
    expect_equal(high_rho$rho, 0.8)
    
    #checking names of output
    expect_named(u1, c("D", "C", "D_Sigma", "rho"))
    
    #checking vector normal method
    v1<- gen_D(method="v", n=200, rho_cond=0.8, s_d1_cond=1, s_d2_cond=1, k=3, 
               C_mu=rep(0, 3), C_cov=0.2, C_var=1, 
               d1_beta=rep(0, 3), d2_beta=rep(0, 3), seed=1)
    v2 <- gen_D(method="v", n=200, rho_cond=0.8, s_d1_cond=1, s_d2_cond=1, k=3, 
                C_mu=rep(0, 3), C_cov=0.2, C_var=1, 
                d1_beta=rep(0, 3), d2_beta=rep(0, 3), seed=1)
    expect_equal(v1$D, v2$D)
    
    #checking matrix normal method
    m1<- gen_D(method="m", n=200, rho_cond=0.8, s_d1_cond=1, s_d2_cond=1, k=3, 
               C_mu=rep(0, 3), C_cov=0.2, C_var=1, 
               d1_beta=rep(0, 3), d2_beta=rep(0, 3), seed=1)
    m2 <- gen_D(method="m", n=200, rho_cond=0.8, s_d1_cond=1, s_d2_cond=1, k=3, 
                C_mu=rep(0, 3), C_cov=0.2, C_var=1, 
                d1_beta=rep(0, 3), d2_beta=rep(0, 3), seed=1)
    expect_equal(m1$D, m2$D)
    
})