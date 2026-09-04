# autorelevate

<!-- badges: start -->
<!-- badges: end -->

`autorelevate` implements the **autorelevated family** of probability
distributions: the family obtained by applying the autorelevation
(self-relevation) transform of Krakowski (1973) and Dileepkumar and
Sankaran (2022) to a baseline lifetime distribution. Ten baseline
families are supported: Weibull, Lomax, Burr XII, Gompertz,
Log-Logistic, Chen, Exponentiated Exponential, Power Lindley,
Log-normal, and Gamma. The Weibull member (the "Autorelevated Weibull",
ARW) is studied in detail by Dileep Kumar, Shabeer, and Sankaran
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

# Compare all ten baseline families by AIC
compare_families(x)

# A real dataset is bundled with the package
data(bladder_cancer)
ttt_plot(bladder_cancer)
compare_families(bladder_cancer)
```

## Citation

See `citation("autorelevate")` for the package citation and the
methodological references.

## License

MIT

