
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
data <- c("pf" = 902, "pm" = 200, "pv" = 44,
          "pf/pm" = 9, "pf/pv" = 29, "pm/pv" = 1, 
          "pf/pm/pv" = 1)
```

We can test to see if any of the coinfection types occur more or less
than best explained by the data and the independence model.

``` r
library(icer)
#> Registered S3 methods overwritten by 'ggplot2':
#>   method         from 
#>   [.quosures     rlang
#>   c.quosures     rlang
#>   print.quosures rlang
res <- cooccurence_test(data)
#> iter   10 value 38.971224
#> final  value 38.971219 
#> converged
#> No id variables; using all as measure variables
```

![](https://i.imgur.com/l5kTqQt.png)<!-- -->

The resultant plot suggests that independent infection events do not
perfectly explain this data, with the predicted number of pf/pm being
too high, with the boostrapped estimates from the best fitting model
notincluding the observed data well (5%-95% quantile in blue, median in
white and observed data in red).

We can then change what type of model we use to describe the acquisition
of species. For example maybe there is some interference between
falciparum and malariae and the inverse for falciparu and vivax. We can
specify starting values for how the probability of acquiring an
additional strains depends of the first strain acquired. e.g. `k_12` is
the multiplication of the probability of the next infection being `pm`
given the first infection was `pf` and vice versa, i.e. next infection
being `pf` given the first infection was `pm`. (The numbers refer to the
species by alphabetic
order):

``` r
res_interf <- cooccurence_test(data, density_func = icer:::interference, 
                               k_12 = 0.5, k_13 = 2, k_23 = 1)
#> iter   10 value 15.580477
#> iter   20 value 13.643069
#> iter   30 value 13.516104
#> iter   40 value 13.188206
#> iter   50 value 13.092270
#> final  value 13.082618 
#> converged
#> No id variables; using all as measure variables
```

![](https://i.imgur.com/mIGrkTw.png)<!-- -->

That’s a lot better\! We can see what the best estimate of the
frequencies of each species, the mean number of infections (different to
number of different species, i.e. you could have 6 infections but be
Pf/Pv if you have 4 Pf infections and 2 Pv infections), and the values
for the interference.

``` r
res_interf$params
#> $params
#>           pf           pm           pv           mu         size 
#> 7.657616e-01 1.755459e-01 5.869248e-02 1.726855e+00 1.000043e+02 
#>         k_12         k_13         k_23 
#> 7.764434e-03 5.867866e-02 1.000000e-02 
#> 
#> $multinom
#>           pf           pm           pv        pf/pm        pf/pv 
#> 0.7604943863 0.1687349573 0.0367985499 0.0076748673 0.0245925814 
#>        pm/pv     pf/pm/pv 
#> 0.0012177302 0.0004869276
```

-----
