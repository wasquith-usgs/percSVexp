# README (./percSVexp/README.md)
<img src='/www/bluebonnets.jpg' height='160' align="right" />

#### Author:           William H. Asquith
#### Point of contact: William H. Asquith (wasquith@usgs.gov)
#### Date:             January 2020

***

# DESCRIPTION

## Introduction

The **percSVexp** repository contains all computer code and data used in the development of a research manuscript that has been submitted to an academic journal. The code herein is considered experimental and research in nature. The suggested citation of the paper at the time of this documentation follows:

* Asquith, W.H., 2020, The use of support vectors from support vector machines for hydrometeorologic monitoring network analysis: Journal of Hydrology
(in press, January 14, 2020)
[William H. Asquith, 0000-0002-7400-1861; USGS Information Product Data System (IPDS) no. IP-104552 (internal agency tracking)]

The home directory of the repository contains two primary subordinate directories: `./exampapp/` and `./gwanalog`. Each of these contains _R_ source code used to conduct the study of **Support Vector Machines** (SVMs) within the main body of the paper as well as the _R_ source code to conduct the example application.

The scripts referenced herein use several add-on packages to the _R_ language. It is the user's responsibility to install these (_e.g._ `install.packages("kernlab")` for the SVM algorithms). Some general notes follow. The aforementioned paper is written in the **LaTeX** format. That format supports use of macros that can be set to report repeated material to mitigate for verificaiton mistakes. The **LaTeX** variables can be identified inside the `.R` scripts herein by comments such as

```{r}
# In LaTeX source \numgageslessfivepercentSVs is the variable less_than_5
```

Connections such as shown above to the manuscript sources are potentially in understanding how the `.R` scripts contribute the manuscript. The **LaTeX** source for the manuscript is not presently provided in this archive. It is anticipated that at some point, the documentation herein or the code will have reference to table and or page numbers in the final paper.


## Subdirectories of the Source Code Archive

1. Directory `./figures/` --- Contains reference copies of the dated figures that upon successful operation of the scripts in this archive are written (hard wired) into the `./draftfigures/` directory. The figures contained herein are planned to be those used in the aforementioned paper and minor editing in vector editing software is anticipated. To clarify, it is the author's intent to place in this directory the final PDF files submitted to the journal.

2. Directory `./draftfigures/` --- Contains results of the scripts contained herein. This directory is hard-wired in the code through notation such as `../draftfigures` to lift the path up from either `./gwanalog/ProbSVMexp.R` or `./examapp/examapp.R` and down into `./draftfigures/`. A user running the scripts can then compare the output to the reference results in `./figures` with a caveat that some differences could be seen because of seeds on the random number generator in `ProbSVMexp.R` or general nuances of SVM use (`examapp.R`).

3. Directory `./gwanalog/` --- Contains the code `ProbSVMexp.R` (probability support vector machine experiment) that is intended to run stand alone within the `./gwanalog/` as the working directory. The script creates eight portable document format (PDF) files that become figures in the aforementioned paper. The figures are written into the `./draftfigures` directory of this archive.

4. Directory `./examapp/` --- Contains the code `examapp.R` along input data in a subdirectory `./exampapp/data` that are described with their own README files. The directory name is to mean "Example Application" as a section in the aforementioned paper.  The script does not need the contents of `./gwanalog/`.
