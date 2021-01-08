[![GitHub license](https://img.shields.io/github/license/ucsf-ferguson-lab/syndRomics)](https://github.com/ucsf-ferguson-lab/syndRomics/blob/master/LICENSE)
[![GitHub issues](https://img.shields.io/github/issues/ucsf-ferguson-lab/syndRomics)](https://github.com/ucsf-ferguson-lab/syndRomics/issues)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://github.com/ucsf-ferguson-lab/syndRomics/graphs/commit-activity)
[![GitHub release](https://img.shields.io/github/release/ucsf-ferguson-lab/syndRomics)](https://github.com/ucsf-ferguson-lab/syndRomics/releases/)

# syndRomics
R syndRomics package

    **The package is in development, please open an issue if you find errors or have ideas about improving the package**

# Introduction

The package implements functions for helping in the process of disease patterns analysis by means of principal components. These include component visualization, interpretation and stability analysis. The goal of the analysis is to understand complex disease states or patterns as common factors (syndromes) that can be quantified from measured variables through the use of computational and analytic approaches (Ferguson et al., 2011, Torres-Espin et al., 2020). In particular, principal component analysis (PCA) and related multivariate statistics has been used as primary toolkit for extracting disease patterns. For details on the analysis, please see the manuscript:

    Abel Torres-Espin, Austin Chou, J. Russell Huie, Nikos Kyritsis, Pavan Upadhyayula, and Adam R. Ferguson. Reproducible analysis of disease space via principal components: a brief tutorial and R package (syndRomics). Under review

Here you will find the basics of the package. There is more, so we encourage you to explore further.

# Installation

syndRomics package can be installed from GitHub using the following code:

```r
install.packages('devtools')
devtools::install_github(repo = "ucsf-ferguson-lab/syndRomics@*release")
```

This will install the latest release. To install previous versions, you can use:

    devtools::install_github(repo = "ucsf-ferguson-lab/syndRomics@tagName")

For example:

```r
devtools::install_github(repo = "ucsf-ferguson-lab/syndRomics@0.0.1.9000")
```

# Use Example

For an example of using the package we will use a public dataset accessible through the Open Data Commons for Spinal Cord Injury [ODC-SCI:26](http://dx.doi.org/10.7295/W9T72FMZ). The dataset consist on studies of subjects with cervical spinal cord injury and a battery of functional test to assess neurological function after injury. For the example we will be using the data at 6 weeks after injury.

The same data have been used in the past by Adam Ferguson team to perform SCI syndromics analysis:

    Ferguson AR, Irvine K-A, Gensel JC, Nielson JL, Lin A, Ly J, et al. Derivation of Multivariate Syndromic Outcome Metrics for Consistent Testing across Multiple Models of Cervical Spinal Cord Injury in Rats. PLOS ONE. 2013 Mar 27;8(3):e59712.

We will use the following libraries:

```r
library(syndRomics)
library(ggplot2)
library(dplyr)
library(stringr)
```

Loading the data and extracting the variables at 6 weeks after injury. To simplify
the analysis for this example, a list-wise deletion has been performed to deal with missing values. For an example of combining syndRomics with multiple imputation strategies see the paper.


```r
odc.df<-read.csv('odc-sci_26.csv', na.strings = '')
odc.df.filter<-odc.df[-1,str_detect(colnames(odc.df),"_wk6")]
odc.df.filter<-odc.df.filter[complete.cases(odc.df.filter),] #delete missing
colnames(odc.df.filter)<-str_remove(colnames(odc.df.filter), "_wk6")
colnames(odc.df.filter)<-str_remove(colnames(odc.df.filter), "AbsDev")#shorten var names
odc.df.filter[1:10,1:6]
#>    wtChng   RFSL  RFPA StepDistRF LFSL  LFPA
#> 35     13 170.02 22.44   5.087542  134 60.14
#> 36     32 145.07 24.45   2.099426  105 97.71
#> 37     29 166.41 55.35    1.79259  166  76.6
#> 38     55 155.09 75.99   0.336176  151 85.75
#> 39     28 101.25 11.74  10.103212  108 42.23
#> 40     21 165.83  9.47   9.364293  127 84.48
#> 41     18 141.25 16.06   6.347048  112 87.41
#> 42      5 168.76 53.99   2.148853  171 58.63
#> 43     32 140.96 12.25   5.114431  103  84.2
#> 44      8 100.08 11.77  14.979624   98 30.27
```

# PCA

All functions in the package with the *pca* argument would accept the results of running either *prcomp()* for linear PCA, or *Gifi::princals()* for non-linear PCA.

Here we illustrate the package using *prcomp()*. We first perform PCA on the data and calculate the standardized loadings.


```r
pca_data<-as.matrix(apply(odc.df.filter, 2, as.numeric))
pca<-prcomp(pca_data, center = T, scale. = T)
```

We can extract the standardized loadings with a helper function:


```r
original_loadings<-stand_loadings(pca, pca_data)
original_loadings[1:10,1:5]
#>                   PC1         PC2         PC3         PC4         PC5
#> wtChng      0.3694971 -0.34980235 -0.45740732  0.04051693  0.05442085
#> RFSL        0.5725462 -0.45817668 -0.21806491 -0.31595720 -0.10395008
#> RFPA       -0.8910488 -0.05184157 -0.22922424 -0.15418749  0.08447155
#> StepDistRF  0.6514693  0.26603212  0.22732227 -0.14192281  0.51307380
#> LFSL        0.2934467 -0.65236778 -0.22635885 -0.40471445 -0.30454885
#> LFPA       -0.5421113  0.25936100 -0.48804336 -0.44179228  0.07193321
#> StepDistLF  0.7663339  0.32175144 -0.09272541  0.05901254  0.38585925
#> RHSL        0.8549303 -0.29425913 -0.24705967 -0.13588432  0.08150684
#> RHPA       -0.7422290  0.39382052 -0.30838385 -0.17046472  0.21732442
#> StepDistRH  0.2446815  0.19815387 -0.41653460  0.63917743 -0.04635492
```


# How many components to keep?

The first question to respond after a PCA solution has been found is usually about knowing how many PCs are relevant. The goal is to determine the minimal set of components that can be used to describe the disease space. 

The syndRomics package implements a nonparametric permutation test for the variance accounted for (VAF), previously suggested as a hypothesis test of whether a PC is generated at random or not. We can generate the null distribution for such test by means of permutation, which breaks the structure of the data. Note that there are several other methods broadly used to decide the "importance" of components. The number of permutations (P) affects the results and it should be chosen to be big enough.


```r
set.seed(500)
per<-permut_pc_test(pca, pca_data, P=1000, ndim=5, statistic = "VAF")
#> Permuting 1000 times for VAF using permD method
#> 0s[==========================================================] 100% ~remaining  0s
#> 
#> Calculating VAF...
#> 
#> DONE
per$results
#>       original       mean     ci_low    ci_high      pvalue adj.p.value
#> PC1 0.33406747 0.09645861 0.08905394 0.10544937 0.000999001 0.001665002
#> PC2 0.18161756 0.08775306 0.08194524 0.09503754 0.000999001 0.001665002
#> PC3 0.09788590 0.08124190 0.07626416 0.08716409 0.000999001 0.001665002
#> PC4 0.07557533 0.07600810 0.07152236 0.08131586 0.554445554 0.693056943
#> PC5 0.06292209 0.07109182 0.06687637 0.07550236 1.000000000 1.000000000
```

We can see that the three first components are statistically significant by the permutation test.

We plot the original VAF vs. the permuted VAF

```r
plot(per, plot_resample = T)
```

<img width="30%" height= "30%" src="./vignettes/figure/permuted VAF-1.png">

# Component Interpretation (identity)

**permut_pc_test()**

The next step in the analysis is to interpret the meaning of the components. That is, providing the PCs with a "identity". This is done by looking at metrics that measure the impact of PCs into variables. The syndRomic package works with standardized loadings and communalities. Standardized loadings are interpreted as the correlation coefficient between a PC and a variable, indicative of the contribution of variables in magnitude and direction of each PC. Communalities are the sum of squared loadings for each variable for the PCs retained during component selection, and they represent how much of the variance of each variable can be explained by the total number of kept components. These can be interpreted as the impact of a variable in the chosen PCA solution formed by the retained PCs. 

To decide PC identity, usually only the higher loadings are used. The cut-off to determine if a loading is "important" or high enough in the interpretation of a PC can be decided in several ways. Using the knowledge of the researcher, several authors suggested rules of thumb for such threshold. Others have suggested a more analytic approach. The syndRomic package implements a non-parametric permutation test for loadings as previously suggested (Buja and Eyuboglu, 1992; Linting et al., 2011).


```r
s_per<-permut_pc_test(pca, pca_data, P=1000, ndim=3, statistic = 's.loadings', perm.method = 'permV')
#> Permuting 1000 x 18 times for s.loadings using permV method
#> [==========================================================] 100% ~remaining  0s
#> 
#> Calculating loadings...
#> 
#> DONE
s_per$results
#> # A tibble: 54 x 8
#> # Groups:   component [3]
#>    Variables         component original     mean ci_low ci_high   pvalue adj.p.value
#>    <chr>             <chr>        <dbl>    <dbl>  <dbl>   <dbl>    <dbl>       <dbl>
#>  1 BBB_FergTrans     PC1       -0.493   -0.00535 -0.214   0.198 0.000999     0.00128
#>  2 BBB_FergTrans     PC2       -0.644   -0.0117  -0.270   0.236 0.000999     0.00225
#>  3 BBB_FergTrans     PC3       -0.00985  0.00473 -0.405   0.385 0.954        0.954  
#>  4 ForelimbOpenField PC1       -0.511    0.00278 -0.203   0.219 0.000999     0.00128
#>  5 ForelimbOpenField PC2       -0.315   -0.00387 -0.254   0.253 0.0120       0.0196 
#>  6 ForelimbOpenField PC3       -0.0302  -0.00939 -0.375   0.389 0.882        0.954  
#>  7 Groom             PC1       -0.389   -0.00543 -0.206   0.193 0.000999     0.00128
#>  8 Groom             PC2       -0.643   -0.00593 -0.280   0.248 0.000999     0.00225
#>  9 Groom             PC3        0.0248  -0.0128  -0.411   0.371 0.924        0.954  
#> 10 LFPA              PC1       -0.542   -0.00671 -0.226   0.203 0.000999     0.00128
#> # ... with 44 more rows
```

We can see that the cut-off of significance (for alpha =0.05) using the adjusted p value is around |0.25| for PC1, around |0.3| for PC2 and around |0.45| for PC3 based on the permutation test.

We can get a barmap plot of the loadings and the permuted distribution (error bars represents mean +/- 95%CI of permuted samples)


```r
plot(s_per, plot_resample = T)
```

<img width="50%" height= "50%" src="./vignettes/figure/s plot 1-1.png">


The same analysis can be performed for communalities. See *?permut_pc_test* for details

# Component stability

**pc_stability()**

Component stability refers to the sensitivity of the PCs to variation, and it is important to study to understand the robustness of the analysis. The package implements functions to help in the process of studying the component stability to data selection variability using resampling methods such as bootstrapping. A robust PC would have small sensitivity to data variations, and thus the goal of the stability analysis is to determine such sensitivity. The package uses nonparametric bootstrapping as previously described for component stability (Babamoradi et al., 2013; Linting et al., 2007; Timmerman et al., 2007; Zientek & Thompson, 2007).

## Stability measures

Component stability can be measured at the component level or at the level of loadings or communalities. At the component level, the package includes: Congruence coefficient, Pearson's r, root mean square error (RMSE) and Cattell's s index (see *?component_similarity for details*). At the level of loadings and communalities, the package allows for computing nonparametric confident intervals. A narrow confident region means that the loading or communality is stable to resampling variations. The *pc_stability()* function returns the accelerated and bias-corrected 95% confident interval for the loadings by default as it has been shown to be robust in different simulation conditions (Efron, B. 1987). For the similarity metrics, the confident interval by the percentile method is returned.

See *?pc_stability* for more options.


```r
booted_PCA<-pc_stability(pca, pca_data,ndim = 3,B =1000, similarity_metric = "all")
#> Bootstrapping 1000 times
#> [==========================================================] 100% ~remaining  0s
#> Calculating confident intervals by bca method
#> [==========================================================] 100% ~remaining  0s
#> Calculating PC similarities
#> 0s[==========================================================] 100% ~remaining  0s
#> 
#> 
#> Final B iterations: 1000
booted_PCA$results
#> # A tibble: 54 x 6
#>    Variables  component Original_loading   mean ci_low ci_high
#>    <chr>      <chr>                <dbl>  <dbl>  <dbl>   <dbl>
#>  1 wtChng     PC1                  0.369  0.366  0.249   0.477
#>  2 RFSL       PC1                  0.573  0.584  0.412   0.708
#>  3 RFPA       PC1                 -0.891 -0.891 -0.918  -0.851
#>  4 StepDistRF PC1                  0.651  0.652  0.576   0.717
#>  5 LFSL       PC1                  0.293  0.293  0.179   0.403
#>  6 LFPA       PC1                 -0.542 -0.543 -0.622  -0.453
#>  7 StepDistLF PC1                  0.766  0.766  0.707   0.817
#>  8 RHSL       PC1                  0.855  0.853  0.809   0.892
#>  9 RHPA       PC1                 -0.742 -0.743 -0.789  -0.684
#> 10 StepDistRH PC1                  0.245  0.250  0.107   0.370
#> # ... with 44 more rows
```

**Component similarity**

```r
booted_PCA$PC_similarity
#> $similarity_mean
#> # A tibble: 3 x 6
#>   PC    cc_index r_correlation   rmse s_index    s_HP
#>   <chr>    <dbl>         <dbl>  <dbl>   <dbl>   <dbl>
#> 1 1        0.996         0.997 0.0503   0.997 0.00309
#> 2 2        0.995         0.995 0.0474   0.981 0.0650 
#> 3 3        0.992         0.987 0.0418   0.952 0.260  
#> 
#> $similarity_ci_low
#> # A tibble: 3 x 6
#>   PC    cc_index r_correlation   rmse s_index   s_HP
#>   <chr>    <dbl>         <dbl>  <dbl>   <dbl>  <dbl>
#> 1 1        0.993         0.993 0.0306   0.971 0     
#> 2 2        0.989         0.989 0.0292   0.941 0.0278
#> 3 3        0.982         0.971 0.0258   0.889 0.222 
#> 
#> $similarity_ci_high
#> # A tibble: 3 x 6
#>   PC    cc_index r_correlation   rmse s_index   s_HP
#>   <chr>    <dbl>         <dbl>  <dbl>   <dbl>  <dbl>
#> 1 1        0.999         0.999 0.0744       1 0.0278
#> 2 2        0.998         0.998 0.0700       1 0.0833
#> 3 3        0.997         0.995 0.0628       1 0.306
```

**Plotting loadings and CI**

```r
plot(booted_PCA, plot_resample = T)
```

<img width="50%" height= "50%" src="./vignettes/figure/boot 1-1.png">


**Plotting communalities and CI**

```r
plot(booted_PCA, communalities = T, plot_resample = T)
```

<img width="50%" height= "50%" src="./vignettes/figure/boot 2-1.png">

# Let's get visual

Visualizations might aid in the process of component interpretation. The syndRomics package implements 4 different visualizations: syndromic plots, heatmap of loadings and barmap of loadings or communalities, and VAF plots.

Any of these visualization can be obtained from either the results of the *prcomp()* or *princals()* functions, or a data.frame containing the loadings or VAF. That allows for generating the plots from loadings obtained by any other function/package/software. In the case of pass a data.frame to the plotting functions, the structure should be as follow: 

    First column named "Variables" and a column for each PCs. One row per variable. The values are the loadings.

In case a data.frame of loadings is passed to the plotting functions, the VAF argument must be given as a character vector of the same lenght of the PCs to plot (ndim) and format: c("60.1%","20.2%") for example. If the results of *prcomp()* or *princals()* are passed, VAF is internally calculated.

## Syndromic plots

The syndromic plot were designed and first published by Ferguson et al, 2013. The plot consist on a middle convex triangle (the intersection of a Venn diagram with 3 sets) displaying the variance accounted for (VAF) for a PC and arrows pointing to the center of the triangle representing each variable. The width of each arrow and the color are proportional to the standardized loading they represent. The variables to plot are determined by a threshold of “importance” for the loadings (cutoff). Syndromic plots are especially useful to report PC identity in publication.

Note that cutoff can be of different values for each extracted plot.


```r
splot<-syndromic_plot(pca, pca_data, ndim = 3,cutoff = 0.45, text_size = 5)
splot$PC1
```

<img width="50%" height= "50%" src="./vignettes/figure/s plot 3-1.png">

The *syndromic_plot()* function returns a ggplot2 object that can be further tuned using the ggplot2 package. You will notice that the arrows do not look quite right in Windows machines. This is a problem with the plotting device. Saving the plot to pdf makes it nicer and it can be easily changes in case need further tuning for publication.


```r
ggsave(plot = splot$PC1, filename = 'figure/PC1.pdf', width = 9, height = 9)
```

<embed src="./vignettes/figure/PC1.pdf" title="plot of chunk unnamed-chunk-9" alt="plot of chunk unnamed-chunk-9" width="500" height="500" type="application/pdf" />

We can customize some aspects of the plot:

**Colors**

The colors argument take three colors in the order: negative loading, zero loading, positive loading


```r
syndromic_plot(pca, pca_data, ndim=1,cutoff = 0.45, text_size = 5, colors = c("orange", "white","purple"))
#> $PC1
```
<img width="50%" height= "50%" src="./vignettes/figure/s plot 4-1.png">

**Arrow width**

The size of the arrows can be controlled with the *arrow_size_multi* argument. This value is a multiplayer such that all arrows will be proportional to loadings by this factor. Default is 10.


```r
syndromic_plot(pca, pca_data, ndim=1,cutoff = 0.45, text_size = 5, arrow_size_multi = 2)
```
<img width="50%" height= "50%" src="./vignettes/figure/s plot 5-1.png">

```r
syndromic_plot(pca, pca_data, ndim=1,cutoff = 0.45, text_size = 5, arrow_size_multi = 15)
```
<img width="50%" height= "50%" src="./vignettes/figure/s plot 6-1.png">


**Reorder variables**

Variables can be ordered either by loadings or by a specified order. The options for ordering by loadings are (starting at 12 o’clock and moving counterclockwise): 

*   'abs decreasing': plot by decreasing absolute value
*   'abs increasing': plot by increasing absolute value 
*   'decreasing'
*   'increasing’


```r
syndromic_plot(pca, pca_data, ndim=1,cutoff = 0.45, text_size = 5, var_order = "increasing")
```
<img width="50%" height= "50%" src="./vignettes/figure/s plot 7-1.png">


```r
syndromic_plot(pca, pca_data, ndim=1,cutoff = 0.45, text_size = 5, var_order = "abs increasing")
#> $PC1
```
<img width="50%" height= "50%" src="./vignettes/figure/s plot 8-1.png">

Variables can also be reordered as you wish if we specify a vector of numbers in the order that we want (starting at 12 o’clock and moving counterclockwise). For example, alphabetically.


```r
var_names<-colnames(pca_data)
syndromic_plot(pca, pca_data, ndim=1,cutoff = 0.45, text_size = 5, var_order =order(var_names))
#> $PC1
```
<img width="50%" height= "50%" src="./vignettes/figure/s plot 9-1.png">

## Heatmap and Barmap plots

Heatmap and barmap plots are other visualizations of the loadings (and communalities in case of barmaps) offering other options beyond the syndromic plot, also returning ggplot2 objects. The major difference between these two plots and the syndromic plot is that both heatmap and barmap plots display all variables (although a subset can be specified using the *vars* argument). This is particularly useful when there are too many variables to plot that might make syndromic plots too crowded, or to compare loadings between PCs more easily. These plots are implemented in the *heatmap_loading()*, the *barmap_loading()*, and the *barmap_commun()* functions. 


```r
heatmap_loading(pca, pca_data, ndim=1:3, cutoff = 0.45, star_values = T, text_values = F)
```

<img width="50%" height= "50%" src="./vignettes/figure/s plot 10-1.png">

Note that for barmap plots we specify a vector of dimensions to plot (e.g. 1:3)

```r
barmap_loading(pca, pca_data, ndim=1:3, cutoff = 0.45)
```
<img width="50%" height= "50%" src="./vignettes/figure/s plot 11-1.png">

Although we can also plot a single PC

```r
barmap_loading(pca, pca_data, ndim=2, cutoff = 0.45)
```

<img width="50%" height= "50%" src="./vignettes/figure/s plot 12-1.png">

Colors and variable order can be modified as we have seen for the syndromic plots.


```r
barmap_loading(pca, pca_data, ndim=1:3, cutoff = 0.45, colors= c("orange2", "white","purple"))
```

<img width="50%" height= "50%" src="./vignettes/figure/s plot 13-1.png">


```r
heatmap_loading(pca, pca_data, ndim=1:3, cutoff = 0.45, star_values = T, text_values = F,
                colors= c("orange2", "white","purple"), vars = var_names)
```

<img width="50%" height= "50%" src="./vignettes/figure/s plot 14-1.png">


```r
barmap_commun(pca, pca_data, ndim=1:3, colors= c("white","purple"))
```

<img width="50%" height= "50%" src="./vignettes/figure/s plot 15-1.png">

**Values.** Whether to plot the value of the loadings or a star indicating the |loadings| >= cutoff can be controlled by the *star_value* and the *text_value* arguments.


```r
barmap_loading(pca, pca_data, ndim=1:3, cutoff = 0.45, star_values = F,text_values = T)
```

<img width="50%" height= "50%" src="./vignettes/figure/s plot 16-1.png">


```r
heatmap_loading(pca, pca_data, ndim=1:3, cutoff = 0.45, star_values = T, text_values = T,
                colors= c("orange2", "white","purple"), vars = var_names)
```

<img width="50%" height= "50%" src="./vignettes/figure/s plot 17-1.png">

## VAF plots

There are two styles of VAF plots:

**line**


```r
VAF_plot(pca, pca_data, ndim=1:10)
```

<img width="30%" height= "30%" src="./vignettes/figure/s plot 18-1.png">

**VAF reduced**


```r
VAF_plot(pca, pca_data, ndim=1:5, style = "reduced")
```

<img width="30%" height= "30%" src="./vignettes/figure/s plot 19-1.png">

## The plot method

There are two ways to plot the results of the *pc_stability()* and *permut_pc_test()*. The results for these can be pass to the VAF plot or the barmap plot functions or using the generic *plot()* function.

For example, permutation test of loadings can be obtained by either way:


```r
barmap_loading(pca, pca_data, ndim=1:3, resample_ci = s_per$results)
```

<img width="50%" height= "50%" src="./vignettes/figure/s plot 20-1.png">


```r
plot(s_per, plot_resample = T, ndim = 1:3)
```

<img width="50%" height= "50%" src="./vignettes/figure/s plot 21-1.png">

This is because *pc_stability()* and *permut_pc_test()* return an object of class "syndromics" that is recognized by *plot()* and redirected to the relevant plot. The same arguments used in the respective plot function can be passed to *plot()*


```r
plot(s_per, plot_resample = T, ndim = 1:3, colors=c("orange2", "white","purple"))
```

<img width="50%" height= "50%" src="./vignettes/figure/s plot 22.0-1.png">


## Tuning up the plots

The plots generated by the **syndRomics** package are ggplot2 objects that can be further modified using the ggplot2 package. See the package vignette for some ideas!

# More information

You can read more about the specific functions using the *help()* or *?*. A more detailed explanation of the implementation can be find in the paper.

    Abel Torres-Espin, Austin Chou, J. Russell Huie, Nikos Kyritsis, Pavan Upadhyayula, and Adam R. Ferguson. Reproducible analysis of disease space via principal components: a brief tutorial and R package (syndRomics). Under review

# References

*   Ferguson AR, Stück ED, Nielson JL. Syndromics: A Bioinformatics Approach for Neurotrauma Research. Transl Stroke Res. 2011 Dec;2(4):438–54.

*   Buja A, Eyuboglu N. Remarks on Parallel Analysis. Multivar Behav Res. 1992 Oct 1;27(4):509–40

*   Linting M, van Os BJ, Meulman JJ. Statistical Significance of the Contribution of Variables to the PCA solution: An Alternative Permutation Strategy. Psychometrika. 2011 Jul 1;76(3):440–60

*   Linting M, Meulman JJ, Groenen PJF, van der Kooij AJ. Stability of nonlinear principal components analysis: An empirical study using the balanced bootstrap. Psychol Methods. 2007;12(3):359–79. 

*   Babamoradi H, van den Berg F, Rinnan Å. Bootstrap based confidence limits in principal component analysis — A case study. Chemom Intell Lab Syst. 2013 Jan 15;120:97–105. 

*   Timmerman ME, Kiers HAL, Smilde AK. Estimating confidence intervals for principal component loadings: a comparison between the bootstrap and asymptotic results. Br J Math Stat Psychol. 2007 Nov;60(Pt 2):295–314. 

*   Zientek LR, Thompson B. Applying the bootstrap to the multivariate case: Bootstrap component/factor analysis. Behav Res Methods. 2007 May;39(2):318–25.

*   Efron, B. (1987). Better Bootstrap Confidence Intervals. Journal of the American Statistical Association, 82(397), 171–185. https://doi.org/10.1080/01621459.1987.10478410
  
*   Ferguson AR, Irvine K-A, Gensel JC, Nielson JL, Lin A, Ly J, et al. Derivation of Multivariate Syndromic Outcome Metrics for Consistent Testing across Multiple Models of Cervical Spinal Cord Injury in Rats. PLOS ONE. 2013 Mar 27;8(3):e59712.
