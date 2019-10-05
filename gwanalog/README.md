# README (./percSVexp/gwanalog/README.md)

#### Author:           William H. Asquith
#### Point of contact: William H. Asquith (wasquith@usgs.gov)
#### Date:             September, 2019

***

# INTRODUCTION

This directory is designed to document source code used in the study of the importance of support vectors of Support Vector Machines (SVMs) in an draft academic paper (Asquith, 2019).  The script `ProbSVMexp.R` requires the external _R_ packages **kernlab** and **mgcv**. These can be found on the Comprehensive R Archive Network (CRAN) (https://cran.r-project.org). The script is technically demanding and to be seen by users as a recipe of sorts for some research ideas on part of the author.


The author used an interactive session to compute the results and then port them to the **LaTeX** sources of the paper. Specifically, the author's manner of operation to support the Asquith (2019) paper is to change the seed to different values (2,3,4) and rerun the script. This manner of operation means than with each pass various output facts or results are then retained for the narrative in the paper. This is not an ideal situation but further design of interface is beyond the scope of this research software and scope of documentation. The order of operation of the script closely follows the greater narrative and mathematical presentation in the aforementioned paper.

The script reports many helpful messages to the console. These report material that is included in tables or narrative of the paper. The tables are identified by their labels in the LaTeX source of the paper. For example, `tab:wholemodel` refers to what is table 1 in the manuscript. The script is also important because it creates portable document format (PDF) figures along a relative path (up and down) at `../draftfigures`. Therein, the figures `figNN*.pdf` where `NN` is 01, 02, ..., 08 for the eight figure numbers.

The script `ProbSVMexp_withRVM.R` represents abandoned experimental efforts involving research on Relavence Vector Machines (RVMs). The script is only kept here for some type of legacy reference should someone be interested. It should run out of the box, but be warned, it will overwrite the figures generated aby the `ProbSVMexp.R` script. The `ProbSVMexp_withRVM.R` does not produce any results directly referenced in the paper. No further discussion is made about the script.

There is a file `figs01ab_refmeas.pdf` that represents a combination of the `fig01_rawsurface.pdf`, `fig02_rawimage.pdf`, and `fig02_rawlegend.pdf` files that are produced by the script. The `fig01ab_refmeas.pdf` file was created in a vector editing software package. This file is emplaced in this directory as a reference version hold the "measure" (the margins) and axis titles and labels as well as the legend on the far right. In the event that the author needs to repeat the construction of the figure in the research paper, then a reference is available.

# ORDER OF OPERATION

The character string `# STEP no.` is emplaced as comments in the script `ProbSVM.R` to help guide the summaries listed below.

1. STEP 1 --- The number of simulations and utility functions for Nash--Sutcliffe efficiency (`NSE`) and generation of a surface (`xymod`) are created.

2. STEP 2 --- A 200 by 200 grid is defined into the matrix `M` with the surface simulated by the `xymod()` calls.

3. STEP 3 --- Creation of the observation network of 800 sites. Here is it important that the seed is set to 62 so for each operation of the script, the same observation sites are always present. The locations and observations of this network are in the variables `X` and `Y` (location) and `Z` (observation).

4. STEP 4 --- Creation of figures `fig01_rawsurface.pdf` and `fig02_rawimage.pdf`. (Recall that these are written to the `../draftfigures` directory. The first figure shows a contour map with the true surface along with the observation sites where realization of observation are made. The shaded region is a special subset to be used in later figures. The second figure is a terrain map representation. A raw legend of the terrain map is `fig02_rawlegend.pdf`. These three figures are branded as `raw` to indicate that a passage through a vector editing software was used to create the final figures 1 and 2 in the paper.

5. STEP 5 --- Creation of figure `fig03_horzmarginA.pdf` uses the data shaded in figures 1 and 2 where the northing coordinate is flatted to the horizontal axis. A generalized additive model (GAM) (`gam <- mgcv::gam(Zp~s(Xp, bs="tp"))`) is fit and draw (blue line, `col=4`). A SVM is fit (`svm <- kernlab::ksvm(Zp~Xp, epsilon=0.3)`) to the data and drawn (red line, `col=2`). It is important to note the use of the plotting character (`pch`) in the `points()` calls. The support vectors are plotted as filled grey dots and the nonsupport vectors are the open circles.

6. STEP 6 --- Creation of figure `fig04_horzmarginB.pdf` is basically the same as in STEP 5 with the exception that a data gap has been artifically created by the `skip` variable. The GAM is again the blue line and the SVM is the red line.

7. STEP 7 --- Creation of figure `fig05_horzmarginC.pdf` is basically the same as in STEP 6, but what has been added is the demonstration of the stochastic nature of the SVM. It is thought that the authors of **kernlab** use a random starting row in the input data as part of the SVM fitting procedure. This produces, potentially, a different combination of support vectors and nonsupport vectors. The red lines in the figure are the individual SVM passes and the green line prepresents the SVM in figure `fig04_horzmarginB.pdf` that was drawn therein in red.

8. STEP 8 --- Creation of figure `fig06_gam2d.pdf` represents the results of a `plot.gam(GAM2d)` call of the GAM model to the surface using a 2-dimensional smooth on the `X` and `Y` coordinates (`GAM2d <- mgcv::gam(Z~s(X,Y, bs="tp"), data=data.frame(X=X, Y=Y))`). There is a second component. The SVM is fit to model the surface (`SVM <- kernlab::ksvm(Z~X+Y, epsilon=0.1)`) and the support vectors isolated and drawn as concentric red circles on the figure. As the figure creation completes a `message()` like the following results that reports entries (random number seed, Nash-Sutcliffe [NSE], Root-mean-square error [RMSE], and sample size [n]) for a table in the Asquith (2019) manuscript with **LaTeX** label `tab:wholemodel`:
```{r}
  SVM Whole (seed, NSE, RMSE, n): 1, 0.986, 3.033, 369 {tab:wholemodel}
```
Highly related to this figure, the next figure (`fig07_svmgam.pdf`) is created that shows the relation between the GAM predictions to the SVM predictions. The purpose of this figure is to demonstrate a general consistency between the two models. This step completes the the saving of the `Z`, `X`, `Y`, and `M` variables to the `ProbSVMexp.RData` (binary) file. This saving is intended to help the user experiment and study the script's operations without having to run the whole script each time. However, this is intended for interactive use by the user experimenting and not as some type of information delivery mechanism.

9. STEP 9 --- Performs computations to evaluate the performance of a GAM if only the support vectors were the "sample": `mgcv::gam(Z[six]~s(X,Y, bs="tp")`. Various results of NSE and RMSE are then reported for the GAM using the sample (support vectors in `six`) as well as those for the nonsupport vectors as `-six`. The **LaTeX** label `tab:submodel` within this step represent results heading to that table in the Asquith (2019) manuscript. The key code components to pay attention to are the `Z[six]` and `Z[-six]` subsettings of the observations for the support and nonsupport vectors, respectively.

10. STEP 10 --- Performs computations to study the performance of random samples to feed to GAM modeling of the surface of sample sizes equal the number of support vectors (`svm.n`) identified earlier. There is no use of an SVM in the step. The **LaTeX** label `tab:simrmsense` within this step represent results heading to that table in the Asquith (2019) manuscript.

## Concerning the Use of Seed Setting on the Pseudo-Random Number Generator

The script is set up with a seeds on the pseudo-random number generated as `set.seed()` for purposes of trying foster reproducibility. The author uses a `set.seed(62)` associated with the "observation network" so that for each run of the script the same 800 locations are used. The `set.seed(seed)` calls reset back when need to a starting state for the seed. Users will note the strategic use of `set.seed(seed)` in several locations in the script and not just a one time setting.

Late testing (October 2019) seems to indicate that the seed setting is not holding some output results entirely constant compared to the period of original code development in January 2019. Concerning the author's computational infrastructure, R versions at the minimum have changed or some other nuance is present. To add in showing reproducibility of sorts, several screen shots are shown. The first two both use a seed of 1 and the next three use seeds of 2, 3, and 4, respectively. The files are located at `./percSVexp/www` with file names `ProbSVMexp_Seed1runA.png`, `ProbSVMexp_Seed1runB.png`, `ProbSVMexp_Seed2run.png`, `ProbSVMexp_Seed3run.png`, and `ProbSVMexp_Seed4run.png`.

### Seed = 1 (Run A)

<img src='/www/ProbSVMexp_Seed1runA.png'  width='900' align="middle" />

### Seed = 1 (Run B)

<img src='/www/ProbSVMexp_Seed1runB.png'  width='900' align="middle" />

### Seed = 2

<img src='/www/ProbSVMexp_Seed2run.png'  width='900' align="middle" />

### Seed = 3

<img src='/www/ProbSVMexp_Seed3run.png'  width='900' align="middle" />

### Seed = 4

<img src='/www/ProbSVMexp_Seed4run.png' width='900' align="middle" />




# SESSION INFO (Fri Sep 27 09:04:02 2019)
```{r}
  R version 3.6.1 (2019-07-05)
  Platform: x86_64-apple-darwin15.6.0 (64-bit)
  Running under: macOS High Sierra 10.13.6

  Matrix products: default
  BLAS:     /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/
                                       vecLib.framework/Versions/A/libBLAS.dylib
  LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib

  locale:
  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

  attached base packages:
  [1] stats     graphics  grDevices utils     datasets  methods   base     

  other attached packages:
  [1] mgcv_1.8-29         nlme_3.1-140        feather_0.3.5       GISTools_0.7-4     
  [5] rgeos_0.5-1         MASS_7.3-51.4       RColorBrewer_1.1-2  maptools_0.9-5     
  [9] rgdal_1.4-4         sp_1.3-1            kernlab_0.9-27      dataRetrieval_2.7.5

 loaded via a namespace (and not attached):
  [1] Rcpp_1.0.2      rstudioapi_0.10 xml2_1.2.2      splines_3.6.1   hms_0.5.1      
  [6] lattice_0.20-38 R6_2.4.0        lmomco_2.3.2    rlang_0.4.0     httr_1.4.1     
 [11] tools_3.6.1     Lmoments_1.3-1  grid_3.6.1      goftest_1.1-1   tibble_2.1.3   
 [16] crayon_1.3.4    Matrix_1.2-17   readr_1.3.1     vctrs_0.2.0     curl_4.2       
 [21] zeallot_0.1.0   compiler_3.6.1  pillar_1.4.2    backports_1.1.4 foreign_0.8-71 
 [26] pkgconfig_2.0.3
```

# REFERENCES

Asquith, W.H., 2019, Assessing Site Importance using Support Vectors for Hydrometeorologic Network Analyses: in reconciliation review with Journal of Hydrology. [William H. Asquith, 0000-0002-7400-1861; USGS Information Product Data System (IPDS) no. IP-104552 (internal agency tracking)]

