# R Installation - Windows / Debian / Ubuntu Version

**Last tested : R 3.5.0, 2018/05/25 (May 21, 2018)**

R packages for installation, the Windows / Debian / Ubuntu version.

Works well for Windows. Works well for Linux (Debian/Ubuntu-like).

This document helps you install over 1,000 packages.

Validation on Windows Subsystem for Linux:

| Operating System | Success | R Version |
| --- | --- | --- |
| Ubuntu 16.04 | :heavy_check_mark: Pass! | :100: R 3.5.0, R 3.4.4 |

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
| SUSE 12 | :heavy_check_mark: Pass! (Not public) | :100: R 3.5.0 |

## Windows Subsystem for Linux (WSL)

<details><summary>REVEAL Windows Subsystem for Linux (WSL) steps</summary>
<p>

Pre-requisites: activate Windows Subsystem for Linux in Additional Features.

### Step 1: install Ubuntu in Windows Subsystem for Linux (WSL)

**Open PowerShell as an Administrator** (right-click the Windows icon on the taskbar > Windows PowerShell Admin):

```ps
Invoke-WebRequest -Uri https://aka.ms/wsl-ubuntu-1604 -OutFile Ubuntu.zip -UseBasicParsing
Expand-Archive Ubuntu.zip Ubuntu
Ubuntu/ubuntu.exe
```

### Step 2: update Ubuntu core packages and install JDK

To perform in Bash shell.

```sh
sudo apt-get update
sudo apt-get upgrade
sudo apt-get install openjdk-9-jdk-headless
sudo apt-get install htop sysstat gedit
export DISPLAY=localhost:0.0
```

### Step 3: do all the Anaconda stuff required

To perform in Bash shell.

```sh
curl -O https://repo.continuum.io/archive/Anaconda3-5.1.0-Linux-x86_64.sh
bash Anaconda3-5.1.0-Linux-x86_64.sh
conda update conda
conda create -n r-tensorflow anaconda
source activate r-tensorflow
pip install --upgrade pip
pip install --ignore-installed --upgrade tensorflow==1.6.0
pip install h5py requests pyyaml Pillow
pip install keras==2.1.5
source deactivate
```

### Step 4: Update bashrc

We are editing `~/.bashrc` file using `gedit`. You might want to use `vi` directly actually.

**Note that we are using `/mnt/e/WSL/R-lib/R-3.5.0` as our R library directory for easy management. Change it to what you want (make sure to Ctrl+F all the required places to change it)**

To perform in Bash shell.

```sh
gedit ~/.bashrc
export DISPLAY="localhost:0.0"
export JAVA_HOME=/usr/lib/jvm/java-9-openjdk-amd64
export PATH=$JAVA_HOME/bin:$PATH
export PATH="$PATH:/usr/local/lib/R"
export R_LIBS_USER="/mnt/e/WSL/R-lib/R-3.5.0"
export R_HOME="/usr/local/lib/R"
export PATH="$PATH:/path/to/anaconda3/bin"
```

Then reset `~/.bashrc`:

```sh
source ~/.bashrc
```

### Step 5: Prepare R pre-requisites

To perform in Bash shell.

```sh
sudo apt-get install libcurl4-openssl-dev libreadline-dev libbz2-dev cmake libxml2-dev git-core libssl-dev libssh2-1-dev gdebi-core libwebp-dev libprotobuf-dev libpoppler-cpp-dev libcairo2-dev librsvg2-dev libv8-3.14-dev libgtk2.0-dev default-jre default-jdk libgmp3-dev libgsl-dev jags libudunits2-dev protobuf-compiler mesa-common-dev libglu1-mesa-dev coinor-libsymphony-dev libtiff5-dev tcl-dev tk-dev libmpfr-dev ggobi libgdal-dev libglpk-dev libgeos-dev netcdf-bin libfftw3-dev libopenmpi-dev bwidget mpi-default-bin libx11-dev ratfor libproj-dev libmagick++-dev coinor-libsymphony-dev coinor-libcgl-dev
```

`jq` specific for geo packages:

```sh
sudo add-apt-repository -y ppa:opencpu/jq
sudo apt-get update
sudo apt-get install libjq-dev
```

`gdal` specific for geo packages:

```sh
sudo add-apt-repository "deb http://ppa.launchpad.net/ubuntugis/ppa/ubuntu $(lsb_release -sc) main"
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 314DF160
sudo apt-get update
sudo apt-get install libgdal-dev
```

R documentation preparation:

```sh
sudo apt-get install texlive-latex-base texlive-fonts-extra
```

### Step 6: Compile R

To perform in Bash shell. `--enable-R-shlib` is mandatory for RStudio Server.

```sh
mkdir R
cd R
wget https://cran.r-project.org/src/base/R-3/R-3.5.0.tar.gz
tar zxvf R-3.5.0.tar.gz
cd R-3.5.0
./configure --enable-R-shlib
```

If everything looks OK, compile R. In this example, we are using 64 cores (using `-j 64`), feel free to modify it to your own usage:

```sh
make -j 64
```

We can now install R if `make` worked:

```sh
sudo make install
```

### Step 7: Check R installation (optional)

We can check the R installation. **This will fail on strict tests.** You can modify the `R_HOME` to point to your R home.

To perform in Bash shell.

```sh
Sys.setenv(LC_COLLATE = "C", LC_TIME = "C", LANGUAGE = "en", R_HOME = "/mnt/e/WSL/R/R-3.5.0")
tools::testInstalledBasic("both")
tools::testInstalledPackages(scope = "base")
tools::testInstalledPackages(scope = "recommended")
```

### Step 8: Install RStudio Server

To perform in Bash shell.

```sh
wget https://download2.rstudio.org/rstudio-server-1.1.453-amd64.deb
sudo gdebi rstudio-server-1.1.453-amd64.deb
```

Inside `/etc/rstudio/rsession.conf`, add the following:

```sh
r-libs-user=/mnt/e/WSL/R-lib/R-3.5.0
session-timeout-minutes=0
```

### Step 9: Firewall protection

If your R Server's IP is publicly available on Internet, it is better to secure the RStudio Server port.

Follow the following steps:

* Windows Defender Firewall with Advanced Security >
* Inbound Rules >
* New Rule (TCP 8787) >
* Properties >
* Scope >
* Remote IP Address >
* 127.0.0.1

You can check by pasting `127.0.0.1:8787` in your Internet browser (use the public IP if front-facing Internet).

### Step 10: Install R packages

You can now install R packages. In a R console from Bash:

And now the long list...:

```r
install.packages("devtools", dependencies = TRUE, Ncpus = parallel::detectCores())
install.packages("tcltk2", dependencies = TRUE, Ncpus = parallel::detectCores())
source("http://bioconductor.org/biocLite.R")
biocLite(c("graph", "RBGL"))

packages <- c("abind", "acepack", "actuar", "ada", "adabag", "ade4", "ade4TkGUI")
install.packages(packages[which(!(packages %in% rownames(installed.packages())))], dependencies = TRUE, Ncpus = parallel::detectCores())

packages <- c("adegraphics", "adehabitatLT", "adehabitatMA", "ADGofTest", "AER", 
              "AGD", "agricolae", "AICcmodavg", "akima", "alabama", "AlgDesign", 
              "alphahull", "alr3", "alr4", "amap", "Amelia", "anchors", "animation", 
              "aod", "aods3", "ape", "aplpack", "argparse", "arm", "arules")
install.packages(packages[which(!(packages %in% rownames(installed.packages())))], dependencies = TRUE, Ncpus = parallel::detectCores())

packages <- c("arulesViz", "ascii", "assertthat", "AUC", "BaBooN", "backports", 
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
              "cvTools", "d3heatmap", "d3Network", "DAAG", "dagitty", "data.table")
install.packages(packages[which(!(packages %in% rownames(installed.packages())))], dependencies = TRUE, Ncpus = parallel::detectCores())

packages <- c("data.tree", "DatABEL", "dataframes2xls", "date", "dbConnect", 
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
              "httr", "huge", "hunspell", "hwriter", "hypergeo", "ibdreg")
install.packages(packages[which(!(packages %in% rownames(installed.packages())))], dependencies = TRUE, Ncpus = parallel::detectCores())

packages <- c("ic.infer", "ICS", "ICSNP", "igraph", "igraphdata", "import", 
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
              "orloca", "orloca.es", "orthopolynom", "outliers", "oz", "packrat")
install.packages(packages[which(!(packages %in% rownames(installed.packages())))], dependencies = TRUE, Ncpus = parallel::detectCores())

packages <- c("pageviews", "pamr", "pan", "pander", "parallelMap", "ParamHelpers", 
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
              "SweaveListingUtils", "systemfit", "tables", "tabplot", "tabplotd3")
install.packages(packages[which(!(packages %in% rownames(installed.packages())))], dependencies = TRUE, Ncpus = parallel::detectCores())

packages <- c("TAM", "tclust", "TeachingDemos", "tensor", "tensorA", 
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
install.packages(packages[which(!(packages %in% rownames(installed.packages())))], dependencies = TRUE, Ncpus = parallel::detectCores())

devtools::install_github("r-lib/progress@a2678e8")
devtools::install_github("Laurae2/woe")
devtools::install_github("Laurae2/xgbdl")
devtools::install_github("Laurae2/lgbdl")
devtools::install_github("twitter/AnomalyDetection")
devtools::install_github("rstudio/tensorflow@a73c8d6")
devtools::install_github("rstudio/keras@bc775ac")
install.packages("reticulate")

devtools::install_github("cmpolis/datacomb", subdir = "pkg", ref = "1.1.2")

# MODIFY THIS FILE TO ADD -O3 to CXXFLAGS: /usr/lib/R/etc/Makeconf
xgbdl::xgb.dl(compiler = "gcc", commit = "8f6aadd", use_avx = TRUE, use_gpu = TRUE)
lgbdl::lgb.dl(commit = "3f54429", compiler = "gcc")
# MODIFY THIS FILE TO PUT BACK -O2 to CXXFLAGS: /usr/lib/R/etc/Makeconf

devtools::install_github("Laurae2/Laurae")
devtools::install_github("Laurae2/LauraeParallel")
devtools::install_github("Laurae2/LauraeDS")
devtools::install_github("Laurae2/LauraeCE")
install.packages("https://cran.r-project.org/src/contrib/Archive/tabplot/tabplot_1.1.tar.gz", repos = NULL, type = "source")
```

</p>
</details>

## Windows Client / Windows Server

<details><summary>REVEAL Windows Client / Server steps</summary>
<p>

### Step 1: Check hardware and software pre-requisites

You will need the following:

* R (R >= 3.4.0, 64-bit **only**) : https://cran.r-project.org/bin/windows/base/
* RStudio : https://www.rstudio.com/products/rstudio/download2/
* MinGW (Rtools, 64-bit **only**) : http://cran.us.r-project.org/bin/windows/Rtools/
* cmake (3.8, 64-bit) : https://cmake.org/files/v3.8/ (required in PATH)
* Git Bash : https://gitforwindows.org/ (required in PATH)
* Java SE Development Kit (JDK) : http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html (download from Oracle)
* Visual Studio 2017 Community with the appropriate SDK (use Windows 10 SDK if you are under Windows 10, Windows 8 SDK if you are under Windows 8 even though Windows 10 SDK is partially retrocompatible) : https://www.visualstudio.com/downloads/ - not required if you don't want accelerated xgboost/LightGBM, nor GPU xgboost.

With GPU support, requires the following extras:

* CUDA 9.0 : https://developer.nvidia.com/cuda-90-download-archive
* CuDNN 7.0 : https://developer.nvidia.com/cudnn
* VC++ 2015 : https://www.visualstudio.com/vs/older-downloads/

The PATH environment variables has the following (make sure to modify the R version accordingly):

```sh
c:\Rtools\bin
c:\Rtools\mingw_64\bin
C:\Program Files\R\R-3.5.0\bin\x64
C:\ProgramData\Oracle\Java\javapath
C:\Program Files\CMake\bin
C:\Program Files\Git\cmd
```

If using GPU, make sure to have the following:

```sh
C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v9.0\bin
C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v9.0\include
C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v9.0\libnvvp
```

### Step 2: Python packages

Windows can download and install Anaconda here: https://repo.anaconda.com/archive/Anaconda3-5.1.0-Windows-x86_64.exe

Run using Anaconda shell:

```py
conda update conda
conda create -n r-tensorflow anaconda
source activate r-tensorflow
pip install --upgrade pip
pip install --ignore-installed --upgrade tensorflow==1.6.0
pip install h5py requests pyyaml Pillow
pip install keras==2.1.5
```

If `activate r-tensorflow` fails, run `source activate r-tensorflow`.

Jupyter Notebook IP casting can be defined using the following in Linux:

```py
jupyter notebook --generate-config
vi ./.jupyter/jupyter_notebook_config.py
```

Then, modify c.NotebookApp.ip = '0.0.0.0' to allow accessing Jupyter Notebook from anywhere (not recommended in Internet).

Esc + : + wq! + Enter can help a lot to quit vi

## Step 3: Install R packages

From an R console (non RStudio), run the following:

```r
install.packages("devtools", dependencies = TRUE, Ncpus = parallel::detectCores())
install.packages("tcltk2", dependencies = TRUE, Ncpus = parallel::detectCores())

source("http://bioconductor.org/biocLite.R")
biocLite(c("graph", "RBGL"))

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
"TAM", "tclust", "TeachingDemos", "tensor", "tensorA", 
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
install.packages(packages[which(!(packages %in% rownames(installed.packages())))], dependencies = TRUE, Ncpus = parallel::detectCores())
```

More R packages:

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

Get even more packages below:

```r
devtools::install_github("cmpolis/datacomb", subdir = "pkg", ref = "1.1.2")
devtools::install_github("ficonsulting/RInno", build_vignettes = TRUE)
```

For xgboost with CPU:

```r
xgbdl::xgb.dl(compiler = "Visual Studio 15 2017 Win64", commit = "8f6aadd", use_avx = TRUE, use_gpu = FALSE)
```

For xgboost with GPU:

```r
xgbdl::xgb.dl(compiler = "Visual Studio 15 2017 Win64", commit = "8f6aadd", use_avx = TRUE, use_gpu = TRUE)
```

Installining now LightGBM...:

```r
lgbdl::lgb.dl(commit = "3f54429", compiler = "vs")
```

Install even more packages...:

```r
devtools::install_github("Laurae2/Laurae")
devtools::install_github("Laurae2/LauraeParallel")
devtools::install_github("Laurae2/LauraeDS")
devtools::install_github("Laurae2/LauraeCE")
install.packages("https://cran.r-project.org/src/contrib/Archive/tabplot/tabplot_1.1.tar.gz", repos=NULL, type="source") # Further versions are too bad / not reliable / generated unreadable plots
```

Confirm Tensorflow works:

```r
library(reticulate)
use_condaenv("r-tensorflow")
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

</p>
</details>

## Ubuntu 17.10

<details><summary>REVEAL Ubuntu 17.10 steps</summary>
<p>

### Step 1: Check hardware and software pre-requisites

To perform in Bash shell.

```sh
sudo apt update
sudo apt upgrade
sudo apt-get install openjdk-9-jre-headless
```

### Step 2: Python packages

Linux can download and install Anaconda using the following:

```r
curl -O https://repo.continuum.io/archive/Anaconda3-5.1.0-Linux-x86_64.sh
bash Anaconda3-5.1.0-Linux-x86_64.sh
```

Run using Anaconda shell:

```py
conda update conda
conda create -n r-tensorflow anaconda
source activate r-tensorflow
pip install --upgrade pip
pip install --ignore-installed --upgrade tensorflow==1.6.0
pip install h5py requests pyyaml Pillow
pip install keras==2.1.5
```

If `activate r-tensorflow` fails, run `source activate r-tensorflow`.

Jupyter Notebook IP casting can be defined using the following in Linux:

```py
jupyter notebook --generate-config
vi ./.jupyter/jupyter_notebook_config.py
```

Then, modify c.NotebookApp.ip = '0.0.0.0' to allow accessing Jupyter Notebook from anywhere (not recommended in Internet).

Esc + : + wq! + Enter can help a lot to quit vi

### Step 3: Use precompiled R

To perform in Bash shell.

```sh
sudo add-apt-repository "deb http://cran.rstudio.com/bin/linux/ubuntu $(lsb_release -sc)/"
gpg --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys E084DAB9
gpg -a --export E084DAB9 | sudo apt-key add -
sudo apt-get update
```

```sh
sudo apt-get install r-base
sudo apt-get install r-cran-rmpi
sudo apt-get install cmake libcurl4-gnutls-dev libcurl4-openssl-dev libxml2-dev git-core libssl-dev libssh2-1-dev gdebi-core openjdk-9-jre-headless libwebp-dev libprotobuf-dev libjq-dev libpoppler-cpp-dev libcairo2-dev librsvg2-dev libv8-3.14-dev libgtk2.0-dev default-jre default-jdk libgmp3-dev libgsl-dev jags libudunits2-dev protobuf-compiler mesa-common-dev libglu1-mesa-dev coinor-libsymphony-dev libtiff5-dev tcl-dev tk-dev libmpfr-dev ggobi libgdal-dev libglpk-dev libgeos-dev netcdf-bin libfftw3-dev libopenmpi-dev bwidget mpi-default-bin
```

And add RStudio Server:

```sh
wget https://download2.rstudio.org/rstudio-server-1.1.442-amd64.deb
sudo gdebi rstudio-server-1.1.442-amd64.deb
```

### Step 4: add R packages

If you get an error on `libjq-dev`, this is normal: it is a Ubuntu 17+ exclusive package.

From an R console (non RStudio), run the following:

```r
install.packages("devtools", dependencies = TRUE, Ncpus = parallel::detectCores())
install.packages("tcltk2", dependencies = TRUE, Ncpus = parallel::detectCores())

source("http://bioconductor.org/biocLite.R")
biocLite(c("graph", "RBGL"))

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
"TAM", "tclust", "TeachingDemos", "tensor", "tensorA", 
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
install.packages(packages[which(!(packages %in% rownames(installed.packages())))], dependencies = TRUE, Ncpus = parallel::detectCores())
```

More R packages:

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

Get even more packages below:

```r
devtools::install_github("cmpolis/datacomb", subdir = "pkg", ref = "1.1.2")
```

For xgboost with CPU:

```r
xgbdl::xgb.dl(compiler = "gcc", commit = "8f6aadd", use_avx = TRUE, use_gpu = FALSE)
```

For xgboost with GPU:

```r
xgbdl::xgb.dl(compiler = "gcc", commit = "8f6aadd", use_avx = TRUE, use_gpu = TRUE)
```

Installining now LightGBM...:

```r
lgbdl::lgb.dl(commit = "3f54429", compiler = "gcc")
```

Install even more packages...:

```r
devtools::install_github("Laurae2/Laurae")
devtools::install_github("Laurae2/LauraeParallel")
devtools::install_github("Laurae2/LauraeDS")
devtools::install_github("Laurae2/LauraeCE")
install.packages("https://cran.r-project.org/src/contrib/Archive/tabplot/tabplot_1.1.tar.gz", repos=NULL, type="source") # Further versions are too bad / not reliable / generated unreadable plots
```

Confirm Tensorflow works:

```r
library(reticulate)
use_condaenv("r-tensorflow")
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

</p>
</details>

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
