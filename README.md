
# SBMStudy

SBMStudy regroups tools to generate Stochastic Block Model (SBM) or
mixture of SBMs and compare multiples techniques of inference.

## Installation

You can install the development version of SBMStudy from
[GitHub](https://github.com/) with:

``` r
install.packages("devtools")
devtools::install_github("Myosis-R/SBMStudy")
```

## Dependencies

``` r
install.packages('purrr')
install.packages('gtools')
install.packages('dplyr')
install.packages('FixedPoint')
install.packages('parallel')
install.packages('doParallel')
install.packages('foreach')
install.packages('itertools')
devtools::install_github('GrossSBM/sbm')
install.packages('graphclust')
devtools::install_github('Chabert-Liddell/colSBM')
```

## Quick walk around SBMStudy

Most useful functions :

- gen_sbm : Compare inference algorithms on a given graph

- gen_mix_sbm : Compare inference algorithms on a collection of graphs

- gen_vfold : Compare inference algorithms on a given graph using
  cross-validation technique

- vem_f : implementation of variationnal expectation-maximization
  algorithm.

- exICL_f : implementation of exact integrated complete likelihood
  algorithm.

## Example 1 : gen_sbm

We analyse inference on a thousand of networks with fifty nodes. The
output of gen_sbm is a dataframe where each row is a new network and
each column is the result of a criterion on a given inference.

Criteria :

- heterogen (measure of the heterogeneity of interactions between
  classes)

- ARIs (Ajusted Rank Index between infered nodes memberships and truth)

- cl_diff (Number of infered classes minus truth)

- exICL (Only for algorithm using exICL, return exICL end value)

Models :

1.  package sbm

2.  package sbm with the constrain of 3 classes

3.  exICL with a random initialization

4.  exICL with another random initialization

5.  exICL with another random initialization

6.  package sbm using exICL in replacement of ICL in order to select
    best model

``` r
suppressPackageStartupMessages(library(tidyverse))

res = SBMStudy::gen_sbm(N=50,R=1000)
glimpse(res)
#> Rows: 1,000
#> Columns: 19
#> $ heterogen1 <dbl> 0.9116106, 0.5717669, 0.5118550, 1.9131156, 0.1733530, 1.58…
#> $ ARI1       <dbl> 1.0000000, 1.0000000, 1.0000000, 1.0000000, 0.9358563, 1.00…
#> $ ARI2       <dbl> 1.0000000, 1.0000000, 1.0000000, 1.0000000, 0.9358563, 1.00…
#> $ ARI3       <dbl> 1.0000000, 0.6832404, 1.0000000, 0.9647853, 0.6114647, 1.00…
#> $ ARI4       <dbl> 1.0000000, 1.0000000, 1.0000000, 0.9647853, 0.5216425, 1.00…
#> $ ARI5       <dbl> 1.0000000, 1.0000000, 0.7788250, 1.0000000, 0.3604172, 1.00…
#> $ ARI6       <dbl> 1.0000000, 1.0000000, 1.0000000, 1.0000000, 0.9358563, 1.00…
#> $ cl_diff1   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
#> $ cl_diff2   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
#> $ cl_diff3   <dbl> 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 2, 0, 1,…
#> $ cl_diff4   <dbl> 0, 0, 0, 1, -1, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, …
#> $ cl_diff5   <dbl> 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,…
#> $ cl_diff6   <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
#> $ exICL1     <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
#> $ exICL2     <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
#> $ exICL3     <dbl> -1387.622, -1524.083, -1530.212, -1122.322, -1513.102, -121…
#> $ exICL4     <dbl> -1387.622, -1509.239, -1530.212, -1122.322, -1508.015, -121…
#> $ exICL5     <dbl> -1387.622, -1509.239, -1544.984, -1112.549, -1530.568, -121…
#> $ exICL6     <dbl> -1387.622, -1509.239, -1530.212, -1112.549, -1494.416, -121…
```

``` r
fig1 = tibble(heterogeneite = res$heterogen1, ARI_VEM = res$ARI1, ARI_ICL = res$ARI3) %>%
  pivot_longer(c(ARI_VEM,ARI_ICL),names_to = 'color', values_to = 'ARI')
ggplot(fig1)+aes(heterogeneite,ARI,color=color)+geom_point()+geom_smooth()
#> `geom_smooth()` using method = 'gam' and formula = 'y ~ s(x, bs = "cs")'
```

<img src="man/figures/README-fig1-1.png" width="100%" />

## Example 2 : gen_mix_sbm

A collection of graphs is generated according to a mixture of SBM,
parameters are contained in thetaMix. If not given thetaMix is generated
randomly at each loop.

``` r
thetaMix = list(
  list(pi=c(0.45,0.55),gamma=matrix(c(.6,.3,.2,.5),2,2),prop=0.5),
  list(pi=c(0.3,0.3,0.4),gamma=matrix(c(.1,.3,.5,.01,.5,.1,.1,.5,.6), 3,3),prop=0.5)
)
res = SBMStudy::gen_mix_sbm(N=30,K=0,M=10,C=0,R=2,thetaMix = thetaMix)
glimpse(res)
#> Rows: 2
#> Columns: 7
#> $ heterogen1  <dbl> 0, 0
#> $ clust_ARIs1 <dbl> 1, 1
#> $ clust_ARIs2 <dbl> 0, 1
#> $ clust_diff1 <dbl> 0, 0
#> $ clust_diff2 <dbl> -1, 0
#> $ nodes_ARIs1 <dbl> 0.200000, 0.650973
#> $ nodes_ARIs2 <dbl> 0.8570959, 0.8848316
```

Criteria :

- heterogen (for each cluster, measure the heterogeneity of interactions
  between classes, return the minimum)

- nodes_ARIs (Average nodes ARI over all graphs)

- clust_ARIs (Ajusted Rank Index between infered network memberships and
  truth)

Models :

1.  package colSBM

2.  package graphclust

## example 3 : gen_vfold

This is a cross-validation method applied to graph. Nodes are divided in
a train group and a test group. inference of nodes membership is done on
the train set. Nodes in the test set are affected to a class thanks to
their connections with the train set. Finally edges probabilities
between nodes in the test set are derived from the infered parameters on
the train set.

``` r
res = SBMStudy::gen_vfold(50,10)
glimpse(res)
#> Rows: 10
#> Columns: 3
#> $ heterogen1 <dbl> 0.38505342, 1.51900805, 1.08422542, 0.69694953, 1.81027991,…
#> $ vfolds1    <dbl> 0.2734110, 0.1548747, 0.1870352, 0.3367920, 0.2424309, 0.47…
#> $ vfolds2    <dbl> 0.2966009, 0.1630528, 0.2089842, 0.2862581, 0.2504147, 0.43…
```

Criteria :

- heterogen (measure the heterogeneity of interactions between classes)

- vfolds (Prediction error over egdes, nodes are split between train and
  test set)

Models :

1.  package sbm

2.  exICL with a random initialization

## Bibliography

- exICL algorithm : P. Louche 2012

- vfold algorithm : L. Zhang 2019

- package sbm : J. Chiquet, S. Donnet, P. Barbillon

- package colSBM : S-C. Chabert-Liddell, P. Barbillon, S. Donnet

- package graphclust : T. Rebafka
