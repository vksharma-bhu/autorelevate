# autorelevate

<!-- badges: start -->
<!-- badges: end -->

`autorelevate` package implements the **autorelevated family** of probability
distributions. The autorelevated family is obtained by applying the autorelevation
transformation of Krakowski (1973) and Dileepkumar and
Sankaran (2022) to a baseline lifetime distribution. Ten baseline
distributions, namely Weibull, Lomax, Burr XII, Gompertz,
Log-Logistic, Chen, Exponentiated Exponential, Power Lindley,
Log-normal, and Gamma constitute the members of the family. The Weibull member ("Autorelevated Weibull")
is studied in detail by Dileep Kumar, Shabeer, and Sankaran
(2025), whose results this package implements and generalizes.

## Installation

```r
# install.packages("devtools")
devtools::install_github("vksharma-bhu/autorelevate")
```

## Example

```r
library(autorelevate)

# Simulate from an Autorelevated Weibull distribution
set.seed(1)
x <- rautorelevate(200, dist = "weibull", p1 = 0.5, p2 = 1.5)

# Fit by maximum likelihood
fit <- fit_autorelevate(x, dist = "weibull", method = "mle")
summary(fit)
plot(fit)

# Compare all ten baseline distributions by information criteria
compare_families(x)

# A real dataset is bundled with the package
data(bladder_cancer)
ttt_plot(bladder_cancer)
compare_families(bladder_cancer)
```

## Citation

See `citation("autorelevate")` to cite the package itself. The
underlying methodological references (Krakowski, 1973; Dileepkumar &
Sankaran, 2022; Dileep Kumar, Shabeer, & Sankaran, 2025; Sharma, Pal,
Bhardwaj, & Tyagi, 2026, submitted) are listed in
`?autorelevate-package` and throughout the function documentation.

## License

MIT

