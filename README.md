
# icer

[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Travis-CI Build
Status](https://travis-ci.org/ojwatson/icer.png?branch=master)](https://travis-ci.org/ojwatson/icer)

## Motivation

`icer` was designed to assess for interference between *Plasmodium*
species to see if some species appeared in infections more often than
expected under the assumption of independent random infection events,
i.e. one species type is introduced from one infectious bite.

`icer` works out the most likely population frequency of each species
and the most likely number of infection types in individual by
comparison to the observed data. These are then used to test if our data
can be explained by this distribution by comparing the observed data to
the bootstrapped samples from the resultant multinomial distribution
describing the probability of being infected with each infection
compositon type.

The approach is somewhat generic and could be used for any coinfection
or cooccurence data in which the observed prevalence of events may not
be predictive of the frequency of the infecting entities.

-----

## Example

Let’s assume we have some data about occurence of *Plasmodium
falciparum*, *malariae* and *vivax* as follows:

``` r
data <- c("pf" = 902, "pm" = 14, "pv" = 44,
          "pf/pm" = 9, "pf/pv" = 29, "pm/pv" = 1, 
          "pf/pm/pv" = 1)
```

We can test to see if any of the coinfection types occur more or less
than best explained by the data and the independence model.

``` r
library(icer)
res <- cooccurence_test(data)
#> final  value 12.355387 
#> converged
#> No id variables; using all as measure variables
```

![](https://i.imgur.com/kewI3d7.png)<!-- -->

The resultant plot suggests that independent infection events explain
this data well, with the boostrapped estimates from the best fitting
model including the observed data well (5%-95% quantile in blue, median
in white and observed data in red).

We also can seee what the best estimate of the frequencies of each
species, the mean number of infections (different to number of different
species, i.e. you could have 6 infections but be Pf/Pv if you have 4 Pf
infections and 2 Pv infections).

``` r
str(res$params)
#> List of 3
#>  $ freqs   : Named num [1:3] 0.9225 0.0192 0.0583
#>   ..- attr(*, "names")= chr [1:3] "pf" "pm" "pv"
#>  $ lambda  : Named num 0.563
#>   ..- attr(*, "names")= chr ""
#>  $ multinom: Named num [1:7] 0.9008 0.0144 0.04413 0.00981 0.03006 ...
#>   ..- attr(*, "names")= chr [1:7] "pf" "pm" "pv" "pf/pm" ...
```

-----

## Extensions Needed

1.  Currently infections is escribed by a zero truncated poisson
    distribution and would be good to include a more dispersed zero
    truncated negative binomial.

2.  Tweak the independence assumption and explore the assumption that
    there could be some suppresion, i.e. if the first strain is Pf, the
    prob of Po is lower and fit this too in our model.

3.  Incorporate alternative explanations for “suppression”. One could
    simply be fine scale spatial heterogeneity in the prevalence of the
    species, i.e. what looks like suppression is just fine spatial
    patterns with different frequencies of each species, and as such
    your more likely to see mono species infection because of local
    hotspots of a particular species.
