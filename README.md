# R Installation - Windows / Debian / Ubuntu Version

**Last tested : R 3.5.0, 2018/05/02 (May 02, 2018)**

R packages for installation, the Windows / Debian / Ubuntu version.

Works well for Windows. Works well for Linux (Debian/Ubuntu-like).

This document helps you install over 1,000 packages.

Validation on Windows operating systems:

| Operating System | Success | R Version |
| --- | --- | --- |
| Windows 10 (1803) | :interrobang: Unknown... | :trident: None yet! |
| Windows 10 (1709) | :heavy_check_mark: Pass! | :100: R 3.5.0, R 3.4.4 |
| Windows 10 (1703) | :heavy_check_mark: Pass! | :100: R 3.5.0, R 3.4.4 |
| Windows 10 (1607) | :heavy_check_mark: Pass! | :100: R 3.5.0, R 3.4.4 |
| Windows 10 (1511) | :heavy_check_mark: Pass! | :100: R 3.5.0, R 3.4.4 |
| Windows 10 (1507) | :heavy_check_mark: Pass! | :100: R 3.5.0, R 3.4.4 |
| Windows Server 2016 | :heavy_check_mark: Pass! | :100: R 3.5.0, R 3.4.4 |
| Windows Server 2012 R2 | :heavy_check_mark: Pass! | :100: R 3.5.0, R 3.4.4 |
| Windows 8.1 | :heavy_check_mark: Pass! | :100: R 3.5.0, R 3.4.4 |
| Windows Server 2012 | :heavy_check_mark: Pass! | :100: R 3.5.0, R 3.4.4 |
| Windows 7 | :heavy_check_mark: Pass! | :100: R 3.5.0, R 3.4.4 |
| Windows Server 2008 R2 | :heavy_check_mark: Pass! | :100: R 3.5.0, R 3.4.4 |
| Windows Vista | :heavy_exclamation_mark: Not fully passing! | :100: R 3.5.0, R 3.4.4 |
| Windows XP | :heavy_exclamation_mark: Not fully passing! | :100: R 3.5.0, R 3.4.4 |

Validation on Ubuntu operating systems:

| Operating System | Success | R Version |
| --- | --- | --- |
| Ubuntu 18.04 | :heavy_exclamation_mark: Fails! | :100: R 3.5.0, R 3.4.4 |
| Ubuntu 17.10 | :heavy_check_mark: Pass! | :100: R 3.5.0, R 3.4.4 |
| Ubuntu 17.04 | :heavy_check_mark: Pass! | :100: R 3.5.0, R 3.4.4 |
| Ubuntu 16.10 | :heavy_check_mark: Pass! | :100: R 3.5.0, R 3.4.4 |
| Ubuntu 16.04 | :heavy_check_mark: Pass! | :100: R 3.5.0, R 3.4.4 |
| Ubuntu 15.10 | :interrobang: Unknown... | :trident: None yet! |
| Ubuntu 14.10 | :interrobang: Unknown... | :trident: None yet! |
| Ubuntu 14.04 | :interrobang: Unknown... | :trident: None yet! |

Validation on SUSE operating systems:

| Operating System | Success | R Version |
| --- | --- | --- |
| SUSE 12 | :heavy_multiplication_x: On the road to test... | :100: R 3.5.0 |

## Windows Pre-Requisites

Windows requires the following for this task:

* Intel CPU, from Nehalem (iX-xxx / Xeon v1, early2011) to Coffee Lake / Scalable Processors series (iX-8xxx / Silver / Gold / Platinum, late2017/early2018)
* R (R >= 3.4.0, 64-bit **only**) : https://cran.r-project.org/bin/windows/base/
* RStudio : https://www.rstudio.com/products/rstudio/download2/
* MinGW (Rtools, 64-bit **only**) : http://cran.us.r-project.org/bin/windows/Rtools/
* cmake (3.8, 64-bit) : https://cmake.org/files/v3.8/ (required in PATH)
* Git Bash : https://gitforwindows.org/ (required in PATH)
* Visual Studio 2017 Community with the appropriate SDK (use Windows 10 SDK if you are under Windows 10, Windows 8 SDK if you are under Windows 8 even though Windows 10 SDK is partially retrocompatible) : https://www.visualstudio.com/downloads/ - not required if you don't want accelerated xgboost/LightGBM, nor GPU xgboost.

If using with NVIDIA GPU:

* CUDA 9.0 : https://developer.nvidia.com/cuda-90-download-archive
* CuDNN 7.0 : https://developer.nvidia.com/cudnn (requires registration, use a throwaway email account)
* C++ Build Tools 2015 Update 3 (choose the appropriate SDK to your Windows Operating System version at the custom installation part, otherwise you screw yourself) : https://www.visualstudio.com/vs/older-downloads/

Non-GPU will activate the following:

* Intel optimized Python (sometimes 1000x faster on very large servers than regular Python, tested using Dual Xeon E5-2697v3)
* Intel optimized Tensorflow 1.6, VS 2017 + AVX wheel (usually 2x faster than regular Tensorflow, please download Visual Studio 2017 Community beforehand, tested using Dual Xeon E5-2697v3)

GPU will activate the following:

* GPU enabled Tensorflow 1.6
* GPU enabled xgboost

## Getting Maximum Performance in R

To get the maximum performance in R (+0 to 10%, depending on the library), please change the following in the file `/etc/x64/Makeconf`:

```
CFLAGS = -O2 -Wall $(DEBUGFLAG) -std=gnu99 -mtune=generic
CXXFLAGS = -O2 -Wall $(DEBUGFLAG) -mtune=generic
CXX98FLAGS = -O2 -Wall $(DEBUGFLAG) -mtune=generic
CXX11FLAGS = -O2 -Wall $(DEBUGFLAG) -mtune=generic
FCFLAGS = -O2 $(DEBUGFLAG) -mtune=generic
FFLAGS = -O2 $(DEBUGFLAG) -mtune=generic
```

to the following:

```
CFLAGS = -O3 -Wall $(DEBUGFLAG) -std=gnu99 -mtune=native
CXXFLAGS = -O3 -Wall $(DEBUGFLAG) -mtune=native
CXX98FLAGS = -O3 -Wall $(DEBUGFLAG) -mtune=native
CXX11FLAGS = -O3 -Wall $(DEBUGFLAG) -mtune=native
FCFLAGS = -O3 $(DEBUGFLAG) -mtune=native
FFLAGS = -O3 $(DEBUGFLAG) -mtune=native
```

If you are getting unexpected mathematical errors while working in R, reinstall R using the default flags. The known highest performance safest flags for R are the following, but please revert to the initial flags instead of these ones:

```
-O2 -msse2 -mfpmath=sse
```

## Linux Pre-Requisites

On Linux (Debian/Ubuntu-like), I don't support GPU installation but you can do it by yourself easily.

Update Linux and install OpenJDK (or any Java SDK you wish).

```r
sudo apt update
sudo apt upgrade
sudo apt-get install openjdk-9-jre-headless
```

Add Ubuntu R repository for the latest R: 

```r
sudo add-apt-repository "deb http://cran.rstudio.com/bin/linux/ubuntu $(lsb_release -sc)/"
```

If using Ubuntu, it is stuborn and requires to bypass the security checks:

```r
gpg --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys E084DAB9
gpg -a --export E084DAB9 | sudo apt-key add -
```

We can now install R and RStudio Server:

```r
sudo apt-get update
sudo apt-get install cmake libcurl4-gnutls-dev libcurl4-openssl-dev libxml2-dev git-core libssl-dev libssh2-1-dev gdebi-core openjdk-9-jre-headless libwebp-dev libprotobuf-dev libjq-dev libpoppler-cpp-dev libcairo2-dev librsvg2-dev libv8-3.14-dev libgtk2.0-dev default-jre default-jdk libgmp3-dev libgsl-dev jags libudunits2-dev protobuf-compiler mesa-common-dev libglu1-mesa-dev coinor-libsymphony-dev libtiff5-dev tcl-dev tk-dev libmpfr-dev ggobi libgdal-dev libglpk-dev libgeos-dev netcdf-bin libfftw3-dev libopenmpi-dev bwidget mpi-default-bin
sudo apt-get install r-base
sudo apt-get install r-cran-rmpi
sudo R CMD javareconf
wget https://download2.rstudio.org/rstudio-server-1.1.442-amd64.deb
sudo gdebi rstudio-server-1.1.442-amd64.deb
```

If you get an error on `libjq-dev`, this is normal: it is a Ubuntu 17+ exclusive package.

## Install Python (Intel Distribution)

For GPU, supposes CUDA 9, CuDNN 7, and Python 3.5 (Anaconda 5.1).

Linux can download and install Anaconda using the following:

```r
curl -O https://repo.continuum.io/archive/Anaconda3-5.1.0-Linux-x86_64.sh
bash Anaconda3-5.1.0-Linux-x86_64.sh
```

Run using Anaconda shell:

```py
conda update conda
conda config --add channels intel
conda create -n idp intelpython3_full python=3
activate idp
pip install --upgrade pip
```

If `activate idp` fails, run `source activate idp`.

Depending on Linux, CPU or GPU, choose.

Linux:

```py
pip install --ignore-installed --upgrade https://anaconda.org/intel/tensorflow/1.6.0/download/tensorflow-1.6.0-cp36-cp36m-linux_x86_64.whl
```

CPU:

```py
pip install --ignore-installed --upgrade https://github.com/fo40225/tensorflow-windows-wheel/raw/master/1.6.0/py36/CPU/avx2/tensorflow-1.6.0-cp36-cp36m-win_amd64.whl
```

GPU:

```py
pip install --ignore-installed --upgrade tensorflow-gpu==1.6.0
```

Keras installation:

```py
pip install h5py requests pyyaml Pillow
conda install scipy -c intel --no-update-deps
pip install keras==2.1.5
```

(Linux specific) Jupyter Notebook IP casting can be defined using the following in Linux:

```py
jupyter notebook --generate-config
vi ./.jupyter/jupyter_notebook_config.py
```

(Linux specific) Then, modify c.NotebookApp.ip = '0.0.0.0' to allow accessing Jupyter Notebook from anywhere (not recommended in Internet).

Esc + : + wq! + Enter can help a lot to quit vi.

## Install lot of packages

If using Linux: test Java compilation so you do not get that absurd `fatal error: jni.h: No such file or directory`:

```r
R CMD javareconf
```

DO NOT DO the following using RStudio. Please use a raw R session (RGui for Windows).

Install two specific BioConductor packages:

```r
source("http://bioconductor.org/biocLite.R")
biocLite(c("graph", "RBGL"))
```

Accept all prompts to force compilation.

For everything below, accept all required prompts. Use 0-Cloud when prompted if you have no idea what you are doing to download packages.

Some packages might fail. Ignore if it is not critical to your use-case.

```r
packages <- c("abind", "acepack", "actuar", "ada", "adabag", "ade4", "ade4TkGUI", 
"adegraphics", "adehabitatLT", "adehabitatMA", "ADGofTest", "AER", 
"AGD", "agricolae", "AICcmodavg", "akima", "alabama", "AlgDesign", 
"alphahull", "alr3", "alr4", "amap", "Amelia", "anchors", "animation", 
"aod", "aods3", "ape", "aplpack", "argparse", "arm", "arules", 
"arulesViz", "ascii", "assertthat", "AUC", "BaBooN", "backports", 
"barcode", "bartMachine", "bartMachineJARs", "base64", "base64enc", 
"BatchJobs", "BayesFactor", "bayesplot", "BayesX", "BayesXsrc", 
"BB", "BBmisc", "bbmle", "BCA", "bcp", "BDgraph", "bdsmatrix", 
"betareg", "BH", "BHH2", "BiasedUrn", "bibtex", "biclust", "biganalytics", 
"biglm", "bigmemory", "bigmemory.sri", "bigRR", "bigtabulate", 
"binda", "bindr", "bindrcpp", "binGroup", "bisoreg", "bit", "bit64", 
"bitops", "blme", "blob", "BMA", "boot", "bootstrap", "Boruta", 
"BradleyTerry2", "breakpoint", "brew", "brglm", "brnn", "broom", 
"BsMD", "bst", "btergm", "C50", "ca", "Cairo", "cairoDevice", 
"CALIBERrfimpute", "calibrator", "candisc", "caper", "car", "CARBayes", 
"CARBayesdata", "carData", "care", "caret", "caretEnsemble", 
"catdata", "caTools", "cba", "ccaPP", "cclust", "CDM", "CDVine", 
"cellranger", "cem", "censReg", "CEoptim", "changepoint", "checkmate", 
"checkpoint", "chemometrics", "chron", "circlize", "CircStats", 
"citr", "Ckmeans.1d.dp", "class", "classInt", "clue", "cluster", 
"clusterSim", "clustvarsel", "clv", "clValid", "cmaes", "cmprsk", 
"coda", "codetools", "coin", "colorplaner", "colorspace", "colourpicker", 
"combinat", "commonmark", "CompQuadForm", "compute.es", "conf.design", 
"config", "contfrac", "contrast", "copula", "CORElearn", "corpcor", 
"corrgram", "corrplot", "covr", "cowplot", "CoxBoost", "coxme", 
"cplm", "crayon", "crosstalk", "crossval", "crp.CSFP", "crrstep", 
"crs", "cshapes", "cubature", "Cubist", "curl", "cvAUC", "CVST", 
"cvTools", "d3heatmap", "d3Network", "DAAG", "dagitty", "data.table", 
"data.tree", "DatABEL", "dataframes2xls", "date", "dbConnect", 
"DBI", "dbscan", "ddalpha", "debugme", "Deducer", "DeducerExtras", 
"deepnet", "degreenet", "deldir", "dendextend", "dendroextras", 
"DendSer", "denstrip", "DEoptim", "DEoptimR", "depthTools", "Deriv", 
"desc", "descr", "DescTools", "deSolve", "Devore7", "devtools", 
"dfoptim", "diagram", "DiagrammeR", "DiagrammeRsvg", "DiceDesign", 
"DiceKriging", "DiceOptim", "dichromat", "digest", "dimRed", 
"diptest", "directlabels", "discretization", "DiscriMiner", "distr", 
"distrEx", "DistributionUtils", "diveMove", "dlm", "DMwR", "doBy", 
"DoE.base", "DoE.wrapper", "doMPI", "doParallel", "doRedis", 
"DoseFinding", "dotCall64", "downloader", "dplR", "dplyr", "drat", 
"DRR", "DT", "dtplyr", "dtw", "dygraphs", "dynamicTreeCut", "dynlm", 
"e1071", "eaf", "earth", "Ecdat", "Ecfun", "ecodist", "effects", 
"eha", "elasticnet", "ElemStatLearn", "ellipse", "elliptic", 
"elmNN", "emdbook", "emoa", "emulator", "energy", "ENmisc", "entropy", 
"EntropyExplorer", "Epi", "EpiModel", "epitools", "erer", "ergm", 
"ergm.count", "ergm.userterms", "eRm", "estimability", "etm", 
"evaluate", "evd", "expint", "ExplainPrediction", "expm", "extrafont", 
"extrafontdb", "extraTrees", "factoextra", "FactoMineR", "Fahrmeir", 
"fail", "faraway", "fAssets", "fastcluster", "fastdigest", "fastICA", 
"fastmatch", "fastR", "fBasics", "fCopulae", "fda", "fdrtool", 
"FeaLect", "feather", "FeatureHashing", "fExoticOptions", "fExtremes", 
"ff", "ffbase", "FFTrees", "fftw", "fGarch", "fields", "filehash", 
"fImport", "findpython", "fit.models", "fitdistrplus", "flare", 
"flashClust", "flexclust", "flexmix", "flexsurv", "FME", "FMStable", 
"fMultivar", "FNN", "fNonlinear", "fontcm", "fOptions", "forcats", 
"foreach", "forecast", "foreign", "formatR", "formattable", "Formula", 
"fortunes", "forward", "fpc", "fPortfolio", "fracdiff", "FRB", 
"frbs", "fRegression", "FrF2", "FrF2.catlg128", "FSelector", 
"fst", "fTrading", "fts", "futile.logger", "futile.options", 
"future", "GA", "gam", "gamair", "GAMBoost", "gamboostLSS", "gamlss", 
"gamlss.data", "gamlss.dist", "gamm4", "gapminder", "gbm", "gclus", 
"gdata", "gdtools", "gee", "geeM", "geepack", "GenABEL", "GenABEL.data", 
"GeneralizedHyperbolic", "genetics", "GenSA", "geoR", "geoRglm", 
"geosphere", "GERGM", "getopt", "GGally", "ggcorrplot", "ggdendro", 
"ggeffects", "ggExtra", "ggformula", "ggfortify", "ggiraph", 
"ggm", "ggplot2", "ggplot2movies", "ggpubr", "ggrepel", "ggsci", 
"ggsignif", "ggThemeAssist", "ggthemes", "ggvis", "giphyr", "git2r", 
"gitgadget", "glasso", "glmmML", "glmmTMB", "glmnet", "glmulti", 
"GlobalOptions", "globals", "glue", "gmailr", "Gmedian", "gmm", 
"gmodels", "gmp", "gnm", "gof", "goftest", "googleVis", "gower", 
"gpairs", "GPArotation", "GPfit", "gplots", "gRbase", "GREA", 
"gridBase", "gridExtra", "grouped", "gsl", "gss", "gstat", "gsubfn", 
"gtable", "gtools", "Guerry", "gWidgets", "gWidgetsRGtk2", "gWidgetstcltk", 
"h2o", "haplo.stats", "haven", "hdi", "heatmaply", "heplots", 
"hergm", "hexbin", "hglm", "hglm.data", "HH", "HiClimR", "highlight", 
"highr", "hmeasure", "Hmisc", "hms", "hrbrthemes", "HSAUR", "HSAUR2", 
"HSAUR3", "htmlTable", "htmltools", "htmlwidgets", "httpuv", 
"httr", "huge", "hunspell", "hwriter", "hypergeo", "ibdreg", 
"ic.infer", "ICS", "ICSNP", "igraph", "igraphdata", "import", 
"imputeTS", "ineq", "influenceR", "Information", "infotheo", 
"inline", "inlinedocs", "intergraph", "intervals", "intsvy", 
"iplots", "ipred", "irace", "irlba", "irr", "isa2", "Iso", "ISOcodes", 
"isotone", "ISwR", "iterators", "itertools", "JavaGD", "JGR", 
"jomo", "jpeg", "jsonlite", "kappalab", "kdecopula", "Kendall", 
"keras", "kernlab", "KernSmooth", "KFAS", "kinship2", "kknn", 
"klaR", "kmi", "knitcitations", "knitr", "kohonen", "koRpus", 
"ks", "labeling", "labelled", "laeken", "LaF", "laGP", "Lahman", 
"lambda.r", "largeVis", "lars", "lasso2", "latentnet", "lattice", 
"latticeExtra", "lava", "lava.tobit", "lavaan", "lavaan.survey", 
"lawstat", "lazyeval", "LCA", "lcopula", "leaflet", "leaps", 
"LearnBayes", "lfda", "lfe", "lhs", "LiblineaR", "likert", "linprog", 
"lintr", "lisrelToR", "listenv", "littleboxes", "lme4", "lmerTest", 
"lmodel2", "lmtest", "loa", "locfit", "logcondens", "LogicReg", 
"logistf", "logspline", "lokern", "longmemo", "loo", "lpSolve", 
"lpSolveAPI", "lqa", "lqmm", "lsmeans", "lubridate", "MAc", "MAd", 
"magrittr", "mail", "manipulate", "mapdata", "mapproj", "maps", 
"maptools", "maptree", "markdown", "MASS", "Matching", "MatchIt", 
"mathgraph", "matlab", "Matrix", "matrixcalc", "MatrixModels", 
"matrixStats", "maxLik", "maxlike", "MBA", "MBESS", "mboost", 
"mc2d", "mcclust", "mcgibbsit", "mclust", "mcmc", "MCMCglmm", 
"MCMCpack", "mco", "mda", "MDSGUI", "mediation", "memisc", "memoise", 
"MEMSS", "merTools", "MetABEL", "metafor", "Metrics", "mets", 
"mgcv", "mi", "mice", "miceadds", "microbenchmark", "microplot", 
"mime", "minerva", "miniUI", "minpack.lm", "minqa", "mirt", "mirtCAT", 
"misc3d", "miscTools", "missForest", "missMDA", "mitml", "mitools", 
"mix", "mlbench", "MLmetrics", "mlmRev", "mlogit", "mlr", "mlrMBO", 
"mnlogit", "mnormt", "modeest", "ModelMetrics", "modelr", "modeltools", 
"mondate", "monreg", "moonBook", "mosaic", "mosaicCalc", "mosaicCore", 
"mosaicData", "movMF", "MplusAutomation", "mpmi", "MPV", "mratios", 
"mRMRe", "msm", "mstate", "MSwM", "muhaz", "multcomp", "multcompView", 
"multicool", "multiwayvcov", "MuMIn", "munsell", "mvinfluence", 
"mvnormtest", "mvoutlier", "mvtnorm", "NbClust", "ncdf4", "ncvreg", 
"ndtv", "network", "networkD3", "networkDynamic", "networkDynamicData", 
"networksis", "neuralnet", "NeuralNetTools", "NHANES", "nlme", 
"nloptr", "NLP", "NMF", "nnet", "nnls", "nodeHarvest", "nor1mix", 
"norm", "nortest", "np", "numbers", "numDeriv", "nws", "nycflights13", 
"obliqueRF", "odfWeave", "officer", "OpenMx", "openssl", "openxlsx", 
"optextras", "optimx", "optmatch", "orcutt", "ordinal", "ore", 
"orloca", "orloca.es", "orthopolynom", "outliers", "oz", "packrat", 
"pageviews", "pamr", "pan", "pander", "parallelMap", "ParamHelpers", 
"partitions", "party", "partykit", "pastecs", "pbapply", "pbivnorm", 
"pbkrtest", "pbmcapply", "PBSmapping", "PBSmodelling", "pcalg", 
"pcaPP", "pec", "penalized", "PerformanceAnalytics", "permute", 
"pgirmess", "pixmap", "pkgconfig", "pkgKitten", "pkgmaker", "PKI", 
"PKPDmodels", "playwith", "plm", "plogr", "plot3D", "plotly", 
"plotmo", "plotrix", "pls", "plyr", "PMCMR", "pmml", "pmmlTransformations", 
"png", "poistweedie", "poLCA", "polspline", "polyclip", "polycor", 
"polynom", "prabclus", "pracma", "praise", "PredictABEL", "prediction", 
"prefmod", "prettyunits", "prim", "pROC", "processx", "prodlim", 
"profdpm", "profileModel", "propagate", "proto", 
"proxy", "pryr", "pscl", "pso", "pspline", "psych", "psychotools", 
"psychotree", "purrr", "pvclust", "pwr", "qap", "qcc", "qgraph", 
"QRAGadget", "qrng", "quadprog", "quantmod", "quantreg", "questionr", 
"qvcalc", "R.cache", "R.devices", "R.matlab", "R.methodsS3", 
"R.oo", "R.rsp", "R.utils", "R2BayesX", "R2Cuba", "R2HTML", "R2jags", 
"R2OpenBUGS", "R2PPT", "R2wd", "R2WinBUGS", "R6", "radiant", 
"radiant.basics", "radiant.data", "radiant.design", "radiant.model", 
"radiant.multivariate", "RandomFields", "RandomFieldsUtils", 
"randomForest", "randomForestSRC", "randtests", "randtoolbox", 
"ranger", "RankAggreg", "RANN", "rappdirs", "RArcInfo", "rARPACK", 
"RaschSampler", "raster", "rasterVis", "rattle", "rbenchmark", 
"rbounds", "rbvs", "Rcgmin", "Rcmdr", "RcmdrMisc", "RcmdrPlugin.BCA", 
"RcmdrPlugin.coin", "RcmdrPlugin.depthTools", "RcmdrPlugin.DoE", 
"RcmdrPlugin.doex", "RcmdrPlugin.epack", "RcmdrPlugin.Export", 
"RcmdrPlugin.FactoMineR", "RcmdrPlugin.HH", "RcmdrPlugin.IPSUR", 
"RcmdrPlugin.KMggplot2", "RcmdrPlugin.mosaic", "RcmdrPlugin.orloca", 
"RcmdrPlugin.pointG", "RcmdrPlugin.qual", "RcmdrPlugin.SLC", 
"RcmdrPlugin.sos", "RcmdrPlugin.steepness", "RcmdrPlugin.survival", 
"RcmdrPlugin.TeachingDemos", "RcmdrPlugin.UCA", "RColorBrewer", 
"Rcpp", "RcppArmadillo", "RcppCNPy", "RcppDE", "RcppEigen", "RcppParallel", 
"RcppProgress", "RcppRoll", "Rcsdp", "RCurl", "readr", "readstata13", 
"readxl", "recipes", "recommenderlab", "ref", "RefManageR", "registry", 
"relaimpo", "relations", "relax", "relevent", "reliaR", "relimp", 
"rem", "rematch", "reportr", "repr", "reshape", "reshape2", "reticulate", 
"rex", "rFerns", "rgdal", "rgenoud", "rgeos", "rgexf", "rggobi", 
"rgl", "Rglpk", "rglwidget", "RgoogleMaps", "RGtk2", "RGtk2Extras", 
"RH2", "rio", "riskRegression", "RItools", "rjags", "rJava", 
"RJDBC", "rjson", "RJSONIO", "rknn", "rlang", "rlecuyer", "rmarkdown", 
"rmeta", "Rmpfr", "Rmpi", "rms", "RMySQL", "rneos", "rngtools", 
"rngWELL", "robCompositions", "robust", "robustbase", "rockchalk", 
"ROCR", "RODBC", "Rook", "rootSolve", "rotationForest", "roxygen2", 
"rpanel", "rpart", "rpart.plot", "rpf", "rpivotTable", "RPostgreSQL", 
"rprojroot", "rrcov", "rredis", "RRF", "rrlda", "RSclient", "rsconnect", 
"Rserve", "RSiena", "RSKC", "rsm", "RSNNS", "Rsolnp", "RSpectra", 
"RSQLite", "rstan", "rstanarm", "rstantools", "rsvg", "Rsymphony", 
"rtiff", "Rtsne", "Rttf2pt1", "rugarch", "RUnit", "Runuran", 
"rversions", "rvest", "rvg", "Rvmmin", "RWeka", "RWekajars", 
"Ryacas", "sampleSelection", "sampling", "sandwich", "scagnostics", 
"scales", "scalreg", "scatterplot3d", "sda", "SEL", "selectr", 
"sem", "semiArtificial", "semPlot", "semTools", "sendmailR", 
"sendplot", "SensoMineR", "seriation", "setRNG", "sets", "sfsmisc", 
"sgeostat", "shape", "shapefiles", "shapes", "shiny", "shinyAce", 
"shinyjs", "shinystan", "shinythemes", "signal", "SimComp", "SimDesign", 
"simecol", "simex", "simsem", "sirt", "SIS", "sjlabelled", "sjmisc", 
"sjPlot", "sjstats", "SkewHyperbolic", "skmeans", "slackr", "slam", 
"SLC", "Sleuth2", "sm", "smbinning", "smoof", "sn", "sna", "snakecase", 
"snow", "SnowballC", "snowfall", "snowFT", "som", "soma", "sos", 
"sourcetools", "sp", "spacetime", "spam", "sparcl", "SparseGrid", 
"sparseLDA", "SparseM", "sparsio", "spatial", "spatstat", "spatstat.data", 
"spatstat.utils", "spc", "spd", "spdep", "speedglm", "sphet", 
"splancs", "splm", "spls", "sqldf", "sROC", "stabledist", "stabs", 
"StanHeaders", "startupmsg", "StatMatch", "statmod", "statnet", 
"statnet.common", "steepness", "stepPlr", "stinepack", "stringdist", 
"stringi", "stringr", "strucchange", "subselect", "subsemble", 
"sudoku", "SuperLearner", "superpc", "SuppDists", "survey", "survival", 
"svd", "svglite", "svGUI", "svUnit", "svyPVpack", "SwarmSVM", 
"SweaveListingUtils", "systemfit", "tables", "tabplot", "tabplotd3", 
"TAM", "tcltk2", "tclust", "TeachingDemos", "tensor", "tensorA", 
"tensorflow", "tergm", "testit", "testthat", "texreg", "tfestimators", 
"tfruns", "tgp", "TH.data", "threejs", "tibble", "tidyr", "tidyselect", 
"tikzDevice", "timeDate", "timereg", "timeSeries", "tis", "tkrplot", 
"tm", "tmap", "TMB", "tmvtnorm", "tnam", "TransferEntropy", "tree", 
"trimcluster", "tripack", "truncdist", "truncnorm", "truncreg", 
"trust", "TSA", "tseries", "tseriesEntropy", "tsna", "TSP", "TTR", 
"tufte", "tuneR", "tweedie", "ucminf", "uniReg", "unmarked", 
"urca", "uuid", "V8", "VarianceGamma", "vars", "vcd", "vcdExtra", 
"Vdgraph", "vegan", "verification", "VGAM", "VGAMdata", "VIM", 
"VIMGUI", "VineCopula", "vioplot", "viridis", "viridisLite", 
"visNetwork", "vtreat", "wavelets", "waveslim", "wbstats", "webp", 
"webshot", "WGCNA", "WhatIf", "whisker", "whoami", "withr", "woe", 
"wordcloud", "WrightMap", "WriteXLS", "wskm", "wsrf", "xergm", 
"xergm.common", "xkcd", "XLConnect", "XLConnectJars", "XML", 
"xml2", "xtable", "xts", "YaleToolkit", "yaml", "yarrr", "zeallot", 
"Zelig", "zip", "zipcode", "zoo", "ztable", "getPass", "lineprof", 
"mapmate", "miniCRAN", "NMOF", "odbc", "recosystem", "redpen", 
"rgeoapi", "rgp", "rgpui", "RSAP", "scrypt", "smooth", "stR", 
"Boom", "BoomSpikeSlab", "bsts", "CausalImpact", "cli", "ClusterR", 
"emmeans", "FD", "fromo", "gdalUtils", "geojson", "geojsonio", 
"geojsonlint", "geometry", "ggridges", "installr", "inum", "jqr", 
"jsonvalidate", "libcoin", "magic", "manipulateWidget", "mapview", 
"moments", "NADA", "OceanView", "OpenImageR", "osmar", "pillar", 
"plot3Drgl", "protolite", "ReporteRs", "ReporteRsjars", "rmapshaper", 
"satellite", "sf", "spData", "SQUAREM", "tiff", "tmaptools", 
"translations", "udunits2", "units", "uroot", "utf8", "xfun", 
"zCompositions")
```

Choose the method below suiting you the best to install the packages.

Failsafe loop:

```r
ret <- lapply(packages, function(x) {
  if (!(x %in% rownames(installed.packages()))) {
    install.packages(x)
  } else {
    cat("Skipping ", x, "\n", sep = "")
  }
})
dput(packages[which(!(packages %in% rownames(installed.packages())))]) # Prints failed packages
```

Install all packages from binaries, Windows:

```r
install.packages(packages, dependencies = TRUE, Ncpus = 8) # Change Ncpus to your number of threads
dput(packages[which(!(packages %in% rownames(installed.packages())))]) # Prints failed packages
```

Install all packages from source, Windows or Linux, please specify number of cores because it could take several hours and/or might end up in an endless loop (it will be an endless loop in Windows in most cases):

```r
install.packages(packages, dependencies = TRUE, type = "source", Ncpus = 8) # Change Ncpus to your number of threads
dput(packages[which(!(packages %in% rownames(installed.packages())))]) # Prints failed packages
```

This might endlessly loop when attempting to install packages one by one:

```r
ret <- lapply(packages, function(x) {
  if (!(x %in% rownames(installed.packages()))) {
    install.packages(x, dependencies = TRUE, type = "source", Ncpus = 8) # Change Ncpus to your number of threads
  } else {
    cat("Skipping ", x, "\n", sep = "")
  }
})
dput(packages[which(!(packages %in% rownames(installed.packages())))]) # Prints failed packages
```

## Attempt Again to Install Packages

This section applies only if you really need to install the failed packages. To attempt to reinstall again the failed packages, you can choose one of the following:

Failsafe loop:

```r
ret <- lapply(packages[which(!(packages %in% rownames(installed.packages())))], function(x) {
  if (!(x %in% rownames(installed.packages()))) {
    install.packages(x)
  } else {
    cat("Skipping ", x, "\n", sep = "")
  }
})
dput(packages[which(!(packages %in% rownames(installed.packages())))]) # Prints failed packages
```

Install all packages from binaries, Windows:

```r
install.packages(packages[which(!(packages %in% rownames(installed.packages())))], dependencies = TRUE, Ncpus = 8) # Change Ncpus to your number of threads
dput(packages[which(!(packages %in% rownames(installed.packages())))]) # Prints failed packages
```

Install all packages from source, Windows or Linux, please specify number of cores because it could take several hours and/or might end up in an endless loop:

```r
install.packages(packages[which(!(packages %in% rownames(installed.packages())))], dependencies = TRUE, type = "source", Ncpus = 8) # Change Ncpus to your number of threads
dput(packages[which(!(packages %in% rownames(installed.packages())))]) # Prints failed packages
```

Suicidal failsafe loop which attempts to install packages one by one:

```r
ret <- lapply(packages[which(!(packages %in% rownames(installed.packages())))], function(x) {
  if (!(x %in% rownames(installed.packages()))) {
    install.packages(x, dependencies = TRUE, type = "source", Ncpus = 8) # Change Ncpus to your number of threads
  } else {
    cat("Skipping ", x, "\n", sep = "")
  }
})
dput(packages[which(!(packages %in% rownames(installed.packages())))]) # Prints failed packages
```

## Extra packages

```r
devtools::install_github("r-lib/progress@a2678e8") # Progress bars with bug fix
devtools::install_github("Laurae2/woe")
devtools::install_github("Laurae2/xgbdl")
devtools::install_github("Laurae2/lgbdl")
devtools::install_github("twitter/AnomalyDetection")
devtools::install_github("rstudio/tensorflow@a73c8d6") # reinstall again
devtools::install_github("rstudio/keras@bc775ac") # reinstall again
install.packages("reticulate") # reinstall again
```

Get even more packages below (success rate not guaranteed for Linux systems, as RInno is Windows only for instance):

```r
devtools::install_github("cmpolis/datacomb", subdir = "pkg", ref = "1.1.2")
devtools::install_github("ficonsulting/RInno", build_vignettes = TRUE)
```

## Windows: Make sure to select the right parameters

Follow Linux instructions if using MinGW only.

Run in RStudio, not Rgui (the xgboost step can fail in Rgui but not in RStudio).

If CPU, you get xgboost enhanced GLM and AVX instructions (change the compiler, commit, and AVX if needed):

```r
xgbdl::xgb.dl(compiler = "Visual Studio 15 2017 Win64", commit = "017acf5", use_avx = TRUE, use_gpu = FALSE)
```

If GPU, you get GPU enabled xgboost, enhanced GLM, and AVX instructions (change the compiler, commit, and AVX if needed):

```r
xgbdl::xgb.dl(compiler = "Visual Studio 14 2015 Win64", commit = "017acf5", use_avx = TRUE, use_gpu = TRUE)
```

Then run for a standard LightGBM installation along with some of my packages to make life easier (change the compiler and commit if needed):

```r
lgbdl::lgb.dl(commit = "4967709", compiler = "vs")
```

```r
devtools::install_github("Laurae2/Laurae")
devtools::install_github("Laurae2/LauraeParallel")
devtools::install_github("Laurae2/LauraeDS")
devtools::install_github("Laurae2/LauraeCE")
install.packages("https://cran.r-project.org/src/contrib/Archive/tabplot/tabplot_1.1.tar.gz", repos=NULL, type="source") # Further versions are too bad / not reliable / generated unreadable plots
```

## Linux: Make sure to select the right parameters

Run in RStudio, not Rgui (the xgboost step can fail in Rgui but not in RStudio).

If CPU, you get xgboost enhanced GLM and AVX instructions:

```r
xgbdl::xgb.dl(compiler = "gcc", commit = "017acf5", use_avx = TRUE, use_gpu = TRUE)
```

If GPU, you get GPU enabled xgboost, enhanced GLM, and AVX instructions:

```r
xgbdl::xgb.dl(compiler = "gcc", commit = "017acf5", use_avx = TRUE, use_gpu = TRUE)
```

Then run for a standard LightGBM installation along with some of my packages to make life easier:

```r
lgbdl::lgb.dl(commit = "4967709", compiler = "gcc")
devtools::install_github("Laurae2/Laurae")
devtools::install_github("Laurae2/LauraeParallel")
devtools::install_github("Laurae2/LauraeDS")
devtools::install_github("Laurae2/LauraeCE")
install.packages("https://cran.r-project.org/src/contrib/Archive/tabplot/tabplot_1.1.tar.gz", repos=NULL, type="source") # Further versions are too bad / not reliable / generated unreadable plots
```

## Confirm Tensorflow works

```r
library(reticulate)
use_condaenv("idp")
reticulate::py_module_available("tensorflow")
library(tensorflow)
sess <- tf$Session()
hello <- tf$constant("OKAY")
sess$run(hello)
```

If not, the error lies here to run in Python:

```py
python
import tensorflow as tf
print(tf.__version__)
okay = tf.constant("OKAY")
sess = tf.Session()
print(sess.run(okay))
```
