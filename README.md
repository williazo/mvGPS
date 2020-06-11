Multivariate Generalized Propensity Score (mvGPS)
================

The goal of this package is to expand currently available software to
estimate weights for multivariate continuous exposures. Weights are
formed assuming a multivariate normal distribution for the simultaneous
exposures.

# Installation

You can install `mvGPS` from GitHub using the following code:

``` r
install.packages("devtools")
devtools::install_github("williazo/mvGPS")
```

## Example

### Data Generating

To illustrate a simple setting where this multivariate generalized
propensity score would be useful, we can construct a simple directed
acyclic graph (DAG) with a bivariate exposure, **D**=(D<sub>1</sub>,
D<sub>2</sub>), confounded by a set of confounders **C** shown below.

<img src="README_files/figure-gfm/dag_draw-1.png" width="50%" style="display: block; margin: auto;" />

To generate this data we first draw *n*=200 samples from **C** assuming
a multivariate normal distribution with mean equal to zero, variance
equal to 1, and constant covariance of 0.3. For this example we assume
that there are 3 confounders

``` r
require(MASS)
require(matrixNormal)
set.seed(06102020)
n <- 200
C_mu <- rep(0, 3)
C_cov <- 0.3
#generating the covariance matrix of 
C_Sigma <- I(3) + ((J(3) - I(3)) * C_cov)
#drawing our observed samples
C <- MASS::mvrnorm(n, mu=C_mu, Sigma=C_Sigma)
```

Next we define our exposure as a linear function of our confounders. In
this case we assume C<sub>1</sub> and C<sub>2</sub> are associated with
D<sub>1</sub>, while C<sub>2</sub> and C<sub>3</sub> are associated with
D<sub>2</sub>. Explicitly these two equations are defined as

E\[D<sub>1</sub>|**C**\]=0.75C<sub>1</sub>+C<sub>2</sub>,

E\[D<sub>2</sub>|**C**\]=C<sub>2</sub>+0.75C<sub>3</sub>.

With this construction, the exposures have one confounder in common,
C<sub>2</sub>, and one independent confounder of equal effect size. We
assume that the conditional distribution of **D** given **C** is
bivariate normal with conditional correlation equal to 0.2 and
conditional variance equal to 2.

``` r
s_d1_cond <- 2
s_d2_cond <- 2
s_d_cond <- c(s_d1_cond, s_d2_cond)
rho_cond <- 0.2

d1_beta <- c(0.75, 1, 0) #exposure D1 effect of C1 and C2
d2_beta <- c(0, 1, 0.75) #exposure D2 effect of C2 and C3
d_beta <- cbind(d1_beta, d2_beta)
d_xbeta <- C %*% d_beta #constructing the conditional mean expression
d1_xbeta <- d_xbeta[, 1]
d2_xbeta <- d_xbeta[, 2]

#construction conditional covariance matrix
D_corr_cond <- I(2) + matrix(c(0, rho_cond, rho_cond, 0), nrow=2, ncol=2, byrow=TRUE)
D_Sigma_cond <- outer(s_d_cond, s_d_cond) * D_corr_cond

#drawing bivariate exposure
D1 <- rnorm(n, d1_xbeta, s_d1_cond)
D2 <- rnorm(n, d2_xbeta + (s_d2_cond / s_d1_cond) * rho_cond * (D1 - d1_xbeta), 
            sqrt((1 - rho_cond^2) * s_d2_cond^2))
D <- cbind(D1, D2)
```

By construction our marginal correlation of D is a function of
parameters from the distribution of **C**, coefficients of conditional
mean equations, and conditional covariance parameter. For the above
specification the true marginal correlation of exposure is equal to 0.4
and our observed marginal correlation is equal to 0.41.

Finally, we specify our outcome, Y, as a linear combination of the
confounders and exposure. The mean of the dose-response equation is
shown below,

E\[Y|**D**,
**C**\]=0.5C<sub>1</sub>+0.3C<sub>2</sub>+0.6C<sub>3</sub>+D<sub>1</sub>+D<sub>2</sub>.

Both exposures have treatment effect sizes equal to one. The standard
deviation of our outcome is set equal 4.

``` r
alpha <- c(0.5, 0.3, 0.6, 1, 1)
sd_Y <- 4
X <- cbind(C, D)
Y <- X%*%alpha + rnorm(n, sd=sd_Y)
```

### Generating Weights

With the data generated, we can now use our primary function `mvGPS()`
to estimate weights. These weights are constructed such that the
numerator is equal to the marginal density, with the denominator
corresponding to the conditional density, i.e., the multivariate
generalized propensity score.

<img src="https://latex.codecogs.com/gif.latex?w=\frac{f(\mathbf{D})}{f(\mathbf{D}\mid\mathbf{C})}" title="w=\frac{f(\mathbf{D})}{f(\mathbf{D}\mid\mathbf{C})}" />

In our case since the bivariate exposure is assumed to be bivariate
normal, we can break both the numerator and denominator into full
conditional densities knowing that each univariate conditional
expression will remain normally distributed.

<img src="https://latex.codecogs.com/gif.latex?w=\frac{f(D_2\mid&space;D_1)f(D_1)}{f(D_2\mid&space;D_1,&space;C_2,&space;C_3)f(D_1\mid&space;C_1,&space;C_2)}" title="w=\frac{f(D_2\mid D_1)f(D_1)}{f(D_2\mid D_1, C_2, C_3)f(D_1\mid C_1, C_2)}" />

Notice in the equation above, we are also able to specify the
confounding set for each exposure separately.

``` r
require(mvGPS)
w <- mvGPS(D=D, C=list(C[, 1:2], C[, 2:3]))
```

This vector w now can be used to test balance of confounders by
comparing weighted vs. unweighted correlations and to estimate the
treatment effects using weighted least squares regression.

### Balance Assessment

For continuous exposure(s) we can asses balance using several metrics
such as euclidean distance, maximum absolute correlation, and average
absolute correlation where correlation refers to the Pearson correlation
between exposure and covariate.

Below we use the function `mvGPSsim_bal()` to specify a set of potential
models to use for comparison. Possible models that are available
include: mvGPS, Entropy, CBPS, GBM, and PS. For methods other than mvGPS
which can only estimate univariate continuous exposure, each exposure is
fit separately so that weights are generated for both exposures.

``` r
require(cobalt)
require(kableExtra)
bal_results <- mvGPSsim_bal(model_list=c("mvGPS", "entropy", "CBPS", "PS", "GBM"), D, C=list(C[, 1:2], C[, 2:3]))
bal_summary <- bal_results$bal_metrics #contains overall summary statistics with respect to balance
knitr::kable(bal_summary, digits=4, row.names=FALSE, 
             col.names=c("Euc. Distance", "Max. Abs. Corr.", "Avg. Abs. Corr.", "Method")) %>%
    kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width=FALSE) %>%
    kableExtra::row_spec(grep("mvGPS", bal_summary$method), bold=TRUE, background="yellow")
```

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:right;">

Euc. Distance

</th>

<th style="text-align:right;">

Max. Abs. Corr.

</th>

<th style="text-align:right;">

Avg. Abs. Corr.

</th>

<th style="text-align:left;">

Method

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;font-weight: bold;background-color: yellow !important;">

0.3126

</td>

<td style="text-align:right;font-weight: bold;background-color: yellow !important;">

0.2499

</td>

<td style="text-align:right;font-weight: bold;background-color: yellow !important;">

0.1329

</td>

<td style="text-align:left;font-weight: bold;background-color: yellow !important;">

mvGPS

</td>

</tr>

<tr>

<td style="text-align:right;">

0.4279

</td>

<td style="text-align:right;">

0.3704

</td>

<td style="text-align:right;">

0.1462

</td>

<td style="text-align:left;">

entropy\_D1

</td>

</tr>

<tr>

<td style="text-align:right;">

0.5267

</td>

<td style="text-align:right;">

0.4095

</td>

<td style="text-align:right;">

0.1852

</td>

<td style="text-align:left;">

entropy\_D2

</td>

</tr>

<tr>

<td style="text-align:right;">

0.4706

</td>

<td style="text-align:right;">

0.4049

</td>

<td style="text-align:right;">

0.1612

</td>

<td style="text-align:left;">

CBPS\_D1

</td>

</tr>

<tr>

<td style="text-align:right;">

0.6016

</td>

<td style="text-align:right;">

0.4524

</td>

<td style="text-align:right;">

0.2122

</td>

<td style="text-align:left;">

CBPS\_D2

</td>

</tr>

<tr>

<td style="text-align:right;">

0.4863

</td>

<td style="text-align:right;">

0.4245

</td>

<td style="text-align:right;">

0.1854

</td>

<td style="text-align:left;">

PS\_D1

</td>

</tr>

<tr>

<td style="text-align:right;">

0.6753

</td>

<td style="text-align:right;">

0.4942

</td>

<td style="text-align:right;">

0.2970

</td>

<td style="text-align:left;">

PS\_D2

</td>

</tr>

<tr>

<td style="text-align:right;">

0.5976

</td>

<td style="text-align:right;">

0.4867

</td>

<td style="text-align:right;">

0.2511

</td>

<td style="text-align:left;">

GBM\_D1

</td>

</tr>

<tr>

<td style="text-align:right;">

0.7264

</td>

<td style="text-align:right;">

0.5429

</td>

<td style="text-align:right;">

0.2652

</td>

<td style="text-align:left;">

GBM\_D2

</td>

</tr>

<tr>

<td style="text-align:right;">

0.9062

</td>

<td style="text-align:right;">

0.5166

</td>

<td style="text-align:right;">

0.4515

</td>

<td style="text-align:left;">

unweighted

</td>

</tr>

</tbody>

</table>

We can see that our method `mvGPS` has the lowest balance metrics across
both exposure dimensions.

### Bias Reduction
