# R Installation - Windows / Debian / Ubuntu Version

**Last tested : R 3.5.0, 2018/05/27 (May 27, 2018)**

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
| Windows 10 (1803) | :heavy_check_mark: Pass! | :100: R 3.5.0, R 3.4.4 |
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

<details><summary>:information_desk_person: CLICK THE ARROW TO REVEAL Windows Subsystem for Linux (WSL) steps</summary>
<p>

Pre-requisites: activate Windows Subsystem for Linux in Additional Features. Make sure to also have VcXsrv installed.

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

To perform in Bash shell. Add Anaconda to PATH, you will need to change it (if you have no idea what are you doing).

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
```

Put the following in the file, change the PATH of R-lib and Anaconda accordingly!:

```sh
export DISPLAY="localhost:0.0"
export JAVA_HOME=/usr/lib/jvm/java-9-openjdk-amd64
export PATH=$JAVA_HOME/bin:$PATH
export PATH="$PATH:/usr/local/lib/R"
export R_LIBS_USER="/mnt/e/WSL/R-lib/R-3.5.0"
export R_HOME="/usr/local/lib/R"
export PATH="$PATH:/home/Lolo/anaconda3/bin"
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

Ignore the errors "rsession: no process found" and "Couldn't find an alternative telinit implementation to spawn.".

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
devtools::install_github("cmpolis/datacomb", subdir = "pkg", ref = "1.1.2")
```

Modify `/usr/lib/R/etc/Makeconf` for maximum performance for xgboost / LightGBM (+0-10% speed): change -O2 to -O3 to CXXFLAGS

Installining xgboost and LightGBM:

```r
xgbdl::xgb.dl(compiler = "gcc", commit = "8f6aadd", use_avx = TRUE, use_gpu = FALSE)
lgbdl::lgb.dl(commit = "3f54429", compiler = "gcc")
```

Modify `/usr/lib/R/etc/Makeconf` for maximum performance for xgboost / LightGBM (+0-10% speed): change back -O3 to -O2 to CXXFLAGS

Install even more packages...:

```r
devtools::install_github("Laurae2/Laurae")
devtools::install_github("Laurae2/LauraeParallel")
devtools::install_github("Laurae2/LauraeDS")
devtools::install_github("Laurae2/LauraeCE")
install.packages("https://cran.r-project.org/src/contrib/Archive/tabplot/tabplot_1.1.tar.gz", repos=NULL, type="source") # Further versions are too bad / not reliable / generated unreadable plots
```

### Step 11: I want to restart from scratch!!!

Run the following as Administrator to fully wipe the Ubuntu installation:

```sh
wslconfig.exe /u Ubuntu
```

</p>
</details>

## Windows Subsystem for Linux (WSL) with Intel Compilers

<details><summary>:information_desk_person: CLICK THE ARROW TO REVEAL Windows Subsystem for Linux (WSL) with Intel Compilers steps</summary>
<p>

Pre-requisites: activate Windows Subsystem for Linux in Additional Features. Make sure to also have VcXsrv installed.

***This will not make R run significantly faster!!!***

**HUGE WARNING: Microsoft intentionally limited WSL stack size limit to 8192 and there is no direct known workaround other than getting into the root user.**

**We are lucky enough to have R running through the root account thanks to the way RStudio Server works.**

**It is also extremely sensitive to the order the R packages are installed, and how many times you tried (once or twice).**

### Step 1: install Ubuntu in Windows Subsystem for Linux (WSL)

**Open PowerShell as an Administrator** (right-click the Windows icon on the taskbar > Windows PowerShell Admin):

```ps
Invoke-WebRequest -Uri https://aka.ms/wsl-ubuntu-1604 -OutFile Ubuntu.zip -UseBasicParsing
Expand-Archive Ubuntu.zip Ubuntu
Ubuntu/ubuntu.exe
```

### Step 2: setup Ubuntu with all libraries required

To perform in Bash shell.

```sh
sudo apt-get update
sudo apt-get upgrade
sudo apt-get install openjdk-9-jdk-headless
sudo apt-get install htop sysstat gedit
sudo apt-get install libcurl4-openssl-dev libreadline-dev libbz2-dev cmake libxml2-dev git-core libssl-dev libssh2-1-dev gdebi-core libwebp-dev libprotobuf-dev libpoppler-cpp-dev libcairo2-dev librsvg2-dev libv8-3.14-dev libgtk2.0-dev default-jre default-jdk libgmp3-dev libgsl-dev jags libudunits2-dev protobuf-compiler mesa-common-dev libglu1-mesa-dev coinor-libsymphony-dev libtiff5-dev tcl-dev tk-dev libmpfr-dev ggobi libgdal-dev libglpk-dev libgeos-dev netcdf-bin libfftw3-dev libopenmpi-dev bwidget mpi-default-bin libx11-dev ratfor libproj-dev libmagick++-dev coinor-libsymphony-dev coinor-libcgl-dev
sudo add-apt-repository -y ppa:opencpu/jq
sudo apt-get update
sudo apt-get install libjq-dev
sudo add-apt-repository "deb http://ppa.launchpad.net/ubuntugis/ppa/ubuntu $(lsb_release -sc) main"
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 314DF160
sudo apt-get update
sudo apt-get install libgdal-dev
sudo apt-get install texlive-latex-base texlive-fonts-extra
export DISPLAY=localhost:0.0
```

### Step 3: Install Intel Compilers

For this example, we are using Intel Parallel Studio XE 2018 Update 2 Cluster Edition.

After you downloaded Intel Parallel Studio, decompress the tgz file then execute the installer:

```sh
tar -xvzf parallel_studio_xe_2018_update2_cluster_edition.tgz
cd parallel_studio_xe_2018_update2_cluster_edition
sudo ./install_GUI.sh
```

Do the following:

* License agreement: check
* Intel Software Improvement Program: decline (I do NOT)
* Prerequisites: skip errors
* License activation: enter your Intel Parallel Studio license
* Options: use only this system (not the cluster), then Customize. Choose whatever you want. VTune will fail, ignore.
* Installation: wait for the massive 10GB+ installation to be done.

Intel Compilers by default, should be in: /opt/intel/compilers_and_libraries_2018.2.199/linux/bin/intel64
/opt/intel/compilers_and_libraries_2018.2.199/linux/bin/compilevars.sh

### Step 4: do all the Anaconda stuff required

To perform in Bash shell. Add Anaconda to PATH, you will need to change it (if you have no idea what are you doing).

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

### Step 5: Update bashrc

We are editing `~/.bashrc` file using `gedit`. You might want to use `vi` directly actually.

**Note that we are using `/mnt/e/WSL/R-lib/R-3.5.0` as our R library directory for easy management. Change it to what you want (make sure to Ctrl+F all the required places to change it)**

To perform in Bash shell.

```sh
gedit ~/.bashrc
```

Put the following in the file, change the PATH of R-lib and Anaconda accordingly!:

```sh
export DISPLAY="localhost:0.0"
export JAVA_HOME=/usr/lib/jvm/java-9-openjdk-amd64
export PATH=$JAVA_HOME/bin:$PATH
export PATH="$PATH:/usr/local/lib/R"
export R_LIBS_USER="/mnt/e/WSL/R-lib/R-3.5.0"
export R_HOME="/usr/local/lib/R"
export PATH="$PATH:/home/Lolo/anaconda3/bin"
source /opt/intel/compilers_and_libraries_2018.2.199/linux/bin/compilervars.sh intel64
export CC="icc"
export CXX="icpc"
export F77="ifort"
export FC="ifort"
export AR="xiar"
export LD="xild"
export CFLAGS="-O3 -ipo -qopenmp -xHost -fPIC"
export CXXFLAGS="-O3 -ipo -qopenmp -xHost -fPIC"
export FFLAGS="-O3 -ipo -qopenmp -xHost -fPIC"
export FCFLAGS="-O3 -ipo -qopenmp -xHost -fPIC"
export MAIN_LDFLAGS="-qopenmp"
export MKL="-L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl"
```

Then reset `~/.bashrc`:

```sh
source ~/.bashrc
```

### Step 6: Workaround stack limit

You will not be able to compile R nor run R without the root user in WSL, because Intel Compilers require a stack limit larger than 8192 for R.

Microsoft added a hard limitation for non root users on the stack limit. It means `ulimit -s unlimited` will always fail (or any value greater than 8192).

Edit the limits.conf file:

```sh
sudo gedit /etc/security/limits.conf
```
Add the following to the file:

```
* hard stack 16384
* soft stack 16384
```

### Step 7: Compile R

To tune Intel MKL (BLAS), use the following link: https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor

Use the following settings as starter:

* Intel MKL
* Linux
* None
* Intel Fortran
* 64-bit
* Dynamic
* 32-bit Integer
* Sequential

To perform in Bash shell. `--enable-R-shlib` is mandatory for RStudio Server.

```sh
mkdir R
cd R
wget https://cran.r-project.org/src/base/R-3/R-3.5.0.tar.gz
tar zxvf R-3.5.0.tar.gz
cd R-3.5.0
```

Now, export some variables. Tune them to your preferences. `-wd308` at CXXFLAGS is mandatory to compile some packages such as `forecast`.

```sh
source /opt/intel/compilers_and_libraries_2018.2.199/linux/bin/compilervars.sh intel64
export CC="icc"
export CXX="icpc"
export F77="ifort"
export FC="ifort"
export AR="xiar"
export LD="xild"
export CFLAGS="-O3 -ipo -qopenmp -xHost -fPIC"
export CXXFLAGS="-O3 -ipo -qopenmp -xHost -fPIC -wd308"
export FFLAGS="-O3 -ipo -qopenmp -xHost -fPIC"
export FCFLAGS="-O3 -ipo -qopenmp -xHost -fPIC"
export LDFLAGS="-qopenmp"
export MKL="-L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl"
```

We can now configure R:

```sh
./configure --enable-R-shlib CC=" icc" CXX=" icpc" FC=" ifort" F77=" ifort" FPICFLAGS=" -fPIC" AR=xiar LD=xild --with-blas="$MKL" --with-lapack
```

This will print you this, it is good:

```sh
  Source directory:          .
  Installation directory:    /usr/local

  C compiler:                 icc  -O3 -ipo -qopenmp -xHost -fPIC
  Fortran 77 compiler:        ifort  -O3 -ipo -qopenmp -xHost -fPIC

  Default C++ compiler:       icpc   -O3 -ipo -qopenmp -xHost -fPIC -wd308
  C++98 compiler:             icpc -std=gnu++98 -O3 -ipo -qopenmp -xHost -fPIC
  C++11 compiler:             icpc -std=gnu++11 -O3 -ipo -qopenmp -xHost -fPIC
  C++14 compiler:
  C++17 compiler:
  Fortran 90/95 compiler:     ifort -O3 -ipo -qopenmp -xHost -fPIC
  Obj-C compiler:

  Interfaces supported:      X11, tcltk
  External libraries:        readline, BLAS(MKL), LAPACK(in blas), curl
  Additional capabilities:   PNG, JPEG, TIFF, NLS, cairo, ICU
  Options enabled:           shared R library, R profiling

  Capabilities skipped:
  Options not enabled:       shared BLAS, memory profiling

  Recommended packages:      yes

configure: WARNING: you cannot build info or HTML versions of the R manuals
```

If everything looks OK, compile R as root. In this example, we are using 64 cores (using `-j 64`), feel free to modify it to your own usage:

```sh
sudo -i
ulimit -s 16384
make -j 64
```

We can now install R if `make` worked:

```sh
make install
```

### Step 7: Check R installation (optional)

Tests will fail using Intel Compilers.

You should get something similar to this:

```sh
make check all

...

running code in 'eval-etc.R' ...Makefile.common:102: recipe for target 'eval-etc.Rout' failed
Attributes: < Component "y": Mean relative difference: 2.029844e-16 >
```

### Step 8: Install RStudio Server

To perform in Bash shell.

```sh
wget https://download2.rstudio.org/rstudio-server-1.1.453-amd64.deb
sudo gdebi rstudio-server-1.1.453-amd64.deb
```

Ignore the errors "rsession: no process found" and "Couldn't find an alternative telinit implementation to spawn.".

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

And now the long list... if it hangs for more than 45 minutes (without anything that seem updating in your console), press Ctrl+C several times and resume from where you were interrupted (you better run block by block the code):

```r
.libPaths("/mnt/e/WSL/R-lib/R-3.5.0")
install.packages("devtools", dependencies = TRUE)
install.packages("tcltk2", dependencies = TRUE)
source("http://bioconductor.org/biocLite.R")
biocLite(c("graph", "RBGL"))

devtools::install_github("davidcsterratt/geometry@d5c98e7", subdir = "pkg")
packages <- c("abind", "acepack", "actuar", "ada", "adabag", "ade4", "ade4TkGUI")
install.packages(packages[which(!(packages %in% rownames(installed.packages())))])

# You will get an error on MCMCpack and Zelig, so we compile MCMCpack by hand
install.packages("https://cran.r-project.org/src/contrib/Archive/MCMCpack/MCMCpack_1.2-4.tar.gz")
packages <- c("adegraphics", "adehabitatLT", "adehabitatMA", "ADGofTest", "AER", 
              "AGD", "agricolae", "AICcmodavg", "akima", "alabama", "AlgDesign", 
              "alphahull", "alr3", "alr4", "amap", "Amelia", "anchors", "animation", 
              "aod", "aods3", "ape", "aplpack", "argparse", "arm", "arules")
install.packages(packages[which(!(packages %in% rownames(installed.packages())))], Ncpus = parallel::detectCores())

packages <- c("arulesViz", "ascii", "assertthat", "AUC", "BaBooN", "backports", 
              "barcode", "base64", "base64enc", "BatchJobs", "BayesFactor", "bayesplot", "BayesX", 
              "BB", "BBmisc", "bbmle", "BCA", "bcp", "BDgraph", "bdsmatrix", 
              "betareg", "BH", "BHH2", "BiasedUrn", "bibtex", "biclust", "biganalytics", 
              "biglm", "bigmemory", "bigmemory.sri", "bigtabulate", 
              "binda", "bindr", "bindrcpp", "binGroup", "bisoreg", "bit", "bit64", 
              "bitops", "blme", "blob", "BMA", "boot", "bootstrap", "Boruta", 
              "BradleyTerry2", "breakpoint", "brew", "brglm", "brnn", "broom", 
              "BsMD", "bst", "C50", "ca", "Cairo", "cairoDevice", 
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
install.packages(packages[which(!(packages %in% rownames(installed.packages())))], Ncpus = parallel::detectCores())

packages <- c("data.tree", "dataframes2xls", "date", "dbConnect", 
              "DBI", "dbscan", "ddalpha", "debugme", 
              "deepnet", "degreenet", "deldir", "dendextend", "dendroextras", 
              "DendSer", "denstrip", "DEoptim", "DEoptimR", "depthTools", "Deriv", 
              "desc", "descr", "DescTools", "deSolve", "Devore7", "devtools", 
              "dfoptim", "diagram", "DiagrammeR", "DiagrammeRsvg", "DiceDesign", 
              "DiceKriging", "DiceOptim", "dichromat", "digest", "dimRed", 
              "diptest", "directlabels", "discretization", "DiscriMiner", "distr", 
              "distrEx", "DistributionUtils", "diveMove", "dlm", "DMwR", "doBy", 
              "DoE.base", "DoE.wrapper", "doParallel", "doRedis", 
              "DoseFinding", "dotCall64", "downloader", "dplR", "dplyr", "drat", 
              "DRR", "DT", "dtplyr", "dtw", "dygraphs", "dynamicTreeCut", "dynlm", 
              "e1071", "eaf", "earth", "Ecdat", "Ecfun", "ecodist", "effects", 
              "eha", "elasticnet", "ElemStatLearn", "ellipse", "elliptic", 
              "elmNN", "emdbook", "emoa", "emulator", "energy", "ENmisc", "entropy", 
              "EntropyExplorer", "Epi", "EpiModel", "epitools", "erer", "ergm", 
              "ergm.count", "ergm.userterms", "eRm", "estimability", "etm", 
              "evaluate", "evd", "expint", "expm", "extrafont", 
              "extrafontdb", "factoextra", "FactoMineR", "Fahrmeir", 
              "fail", "faraway", "fAssets", "fastcluster", "fastdigest", "fastICA", 
              "fastmatch", "fastR", "fBasics", "fCopulae", "fda", "fdrtool", 
              "FeaLect", "feather", "FeatureHashing", "fExoticOptions", "fExtremes", 
              "ff", "ffbase", "FFTrees", "fftw", "fGarch", "fields", "filehash", 
              "fImport", "findpython", "fit.models", "fitdistrplus", "flare", 
              "flashClust", "flexclust", "flexmix", "flexsurv", "FME", "FMStable", 
              "fMultivar", "FNN", "fNonlinear", "fontcm", "fOptions", "forcats", 
              "foreach", "forecast", "foreign", "formatR", "formattable", "Formula", 
              "fortunes", "fpc", "fPortfolio", "fracdiff", "FRB", 
              "frbs", "fRegression", "FrF2", "FrF2.catlg128", 
              "fst", "fTrading", "fts", "futile.logger", "futile.options", 
              "future", "GA", "gam", "gamair", "GAMBoost", "gamboostLSS", "gamlss", 
              "gamlss.data", "gamlss.dist", "gamm4", "gapminder", "gbm", "gclus", 
              "gdata", "gdtools", "gee", "geeM", "geepack", 
              "GeneralizedHyperbolic", "genetics", "GenSA", "geoR", "geoRglm", 
              "geosphere", "GERGM", "getopt", "GGally", "ggcorrplot", "ggdendro", 
              "ggeffects", "ggExtra", "ggformula", "ggfortify", "ggiraph", 
              "ggm", "ggplot2", "ggplot2movies", "ggpubr", "ggrepel", "ggsci", 
              "ggsignif", "ggThemeAssist", "ggthemes", "ggvis", "giphyr", "git2r", 
              "gitgadget", "glmmML", "glmmTMB", "glmnet", 
              "GlobalOptions", "globals", "glue", "gmailr", "Gmedian", "gmm", 
              "gmodels", "gmp", "gnm", "gof", "goftest", "googleVis", "gower", 
              "gpairs", "GPArotation", "GPfit", "gplots", "gRbase", 
              "gridBase", "gridExtra", "grouped", "gsl", "gss", "gstat", "gsubfn", 
              "gtable", "gtools", "Guerry", "gWidgets", "gWidgetsRGtk2", "gWidgetstcltk", 
              "h2o", "haplo.stats", "haven", "hdi", "heatmaply", "heplots", 
              "hergm", "hexbin", "hglm", "hglm.data", "HH", "HiClimR", "highlight", 
              "highr", "hmeasure", "Hmisc", "hms", "hrbrthemes", "HSAUR", "HSAUR2", 
              "HSAUR3", "htmlTable", "htmltools", "htmlwidgets", "httpuv", 
              "httr", "huge", "hunspell", "hwriter", "hypergeo", "ibdreg")
install.packages(packages[which(!(packages %in% rownames(installed.packages())))], Ncpus = parallel::detectCores())

packages <- c("ic.infer", "ICS", "ICSNP", "igraph", "igraphdata", "import", 
              "imputeTS", "ineq", "influenceR", "Information", "infotheo", 
              "inline", "inlinedocs", "intergraph", "intervals", "intsvy", 
              "ipred", "irace", "irlba", "irr", "isa2", "Iso", "ISOcodes", 
              "isotone", "ISwR", "iterators", "itertools", "JavaGD", 
              "jomo", "jpeg", "jsonlite", "kappalab", "kdecopula", "Kendall", 
              "keras", "kernlab", "KernSmooth", "kinship2", "kknn", 
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
              "matrixStats", "maxLik", "maxlike", "MBA", "mboost", 
              "mc2d", "mcclust", "mcgibbsit", "mclust", "mcmc", "MCMCglmm", 
              "MCMCpack", "mco", "mda", "MDSGUI", "mediation", "memisc", "memoise", 
              "MEMSS", "merTools", "MetABEL", "metafor", "Metrics", "mets", 
              "mgcv", "mi", "mice", "miceadds", "microbenchmark", 
              "mime", "minerva", "miniUI", "minpack.lm", "minqa", "mirt", "mirtCAT", 
              "misc3d", "miscTools", "missForest", "missMDA", "mitml", "mitools", 
              "mix", "mlbench", "MLmetrics", "mlmRev", "mlogit", "mlr", "mlrMBO", 
              "mnlogit", "mnormt", "modeest", "ModelMetrics", "modelr", "modeltools", 
              "mondate", "monreg", "moonBook", "mosaic", "mosaicCalc", "mosaicCore", 
              "mosaicData", "movMF", "MplusAutomation", "MPV", "mratios", 
              "mRMRe", "msm", "mstate", "MSwM", "muhaz", "multcomp", "multcompView", 
              "multicool", "multiwayvcov", "MuMIn", "munsell", "mvinfluence", 
              "mvnormtest", "mvoutlier", "mvtnorm", "NbClust", "ncdf4", "ncvreg", 
              "ndtv", "network", "networkD3", "networkDynamic", "networkDynamicData", 
              "networksis", "neuralnet", "NeuralNetTools", "NHANES", "nlme", 
              "nloptr", "NLP", "NMF", "nnet", "nnls", "nodeHarvest", "nor1mix", 
              "norm", "nortest", "np", "numbers", "numDeriv", "nws", "nycflights13", 
              "obliqueRF", "odfWeave", "officer", "openssl", "openxlsx", 
              "optextras", "optimx", "optmatch", "orcutt", "ordinal", "ore", 
              "orloca", "orloca.es", "orthopolynom", "outliers", "oz", "packrat")
install.packages(packages[which(!(packages %in% rownames(installed.packages())))], Ncpus = parallel::detectCores())

packages <- c("pageviews", "pamr", "pan", "pander", "parallelMap", "ParamHelpers", 
              "partitions", "party", "partykit", "pastecs", "pbapply", "pbivnorm", 
              "pbkrtest", "pbmcapply", "PBSmapping", "PBSmodelling", "pcalg", 
              "pcaPP", "pec", "penalized", "PerformanceAnalytics", "permute", 
              "pgirmess", "pixmap", "pkgconfig", "pkgKitten", "pkgmaker", "PKI", 
              "PKPDmodels", "playwith", "plm", "plogr", "plot3D", "plotly", 
              "plotmo", "plotrix", "pls", "plyr", "PMCMR", "pmml", "pmmlTransformations", 
              "png", "poistweedie", "poLCA", "polspline", "polyclip", "polycor", 
              "polynom", "prabclus", "pracma", "praise", "PredictABEL", "prediction", 
              "prettyunits", "prim", "pROC", "processx", "prodlim", 
              "profdpm", "profileModel", "propagate", "proto", 
              "proxy", "pryr", "pscl", "pso", "pspline", "psych", "psychotools", 
              "psychotree", "purrr", "pvclust", "pwr", "qap", "qcc", 
              "QRAGadget", "qrng", "quadprog", "quantmod", "quantreg", "questionr", 
              "qvcalc", "R.cache", "R.devices", "R.matlab", "R.methodsS3", 
              "R.oo", "R.rsp", "R.utils", "R2Cuba", "R2HTML", "R2jags", 
              "R2OpenBUGS", "R2WinBUGS", "R6", "RandomFields", "RandomFieldsUtils", 
              "randomForest", "randtests", "randtoolbox", 
              "ranger", "RankAggreg", "RANN", "rappdirs", "RArcInfo", "rARPACK", 
              "RaschSampler", "raster", "rasterVis", "rattle", "rbenchmark", 
              "rbounds", "rbvs", "Rcgmin", "RColorBrewer", 
              "Rcpp", "RcppArmadillo", "RcppCNPy", "RcppDE", "RcppEigen", "RcppParallel", 
              "RcppProgress", "RcppRoll", "Rcsdp", "RCurl", "readr", "readstata13", 
              "readxl", "recipes", "recommenderlab", "ref", "RefManageR", "registry", 
              "relaimpo", "relations", "relevent", "reliaR", "relimp", 
              "rem", "rematch", "reportr", "repr", "reshape", "reshape2", "reticulate", 
              "rex", "rFerns", "rgdal", "rgenoud", "rgeos", "rgexf", "rggobi", 
              "rgl", "Rglpk", "rglwidget", "RgoogleMaps", "RGtk2", "RGtk2Extras", 
              "RH2", "rio", "riskRegression", "RItools", "rjags", "rJava", 
              "RJDBC", "rjson", "RJSONIO", "rknn", "rlang", "rlecuyer", "rmarkdown", 
              "rmeta", "Rmpfr", "rms", "RMySQL", "rneos", "rngtools", 
              "rngWELL", "robCompositions", "robust", "robustbase", "rockchalk", 
              "ROCR", "RODBC", "Rook", "rootSolve", "rotationForest", "roxygen2", 
              "rpanel", "rpart", "rpart.plot", "rpf", "rpivotTable", "RPostgreSQL", 
              "rprojroot", "rrcov", "rredis", "RRF", "RSclient", "rsconnect", 
              "Rserve", "RSKC", "rsm", "Rsolnp", "RSpectra", 
              "RSQLite", "rstan", "rstanarm", "rstantools", "rsvg", "Rsymphony", 
              "rtiff", "Rtsne", "Rttf2pt1", "rugarch", "RUnit", "Runuran", 
              "rversions", "rvest", "rvg", "Rvmmin", 
              "Ryacas", "sampleSelection", "sampling", "sandwich", 
              "scales", "scalreg", "scatterplot3d", "sda", "SEL", "selectr", 
              "sem", "semTools", "sendmailR", 
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
install.packages(packages[which(!(packages %in% rownames(installed.packages())))], Ncpus = parallel::detectCores())

packages <- c("TAM", "tclust", "TeachingDemos", "tensor", "tensorA", 
              "tensorflow", "tergm", "testit", "testthat", "texreg", "tfestimators", 
              "tfruns", "tgp", "TH.data", "threejs", "tibble", "tidyr", "tidyselect", 
              "tikzDevice", "timeDate", "timereg", "timeSeries", "tis", "tkrplot", 
              "tm", "tmap", "TMB", "tmvtnorm", "tnam", "TransferEntropy", "tree", 
              "trimcluster", "tripack", "truncdist", "truncnorm", "truncreg", 
              "trust", "TSA", "tseries", "tsna", "TSP", "TTR", 
              "tufte", "tuneR", "tweedie", "ucminf", "uniReg", "unmarked", 
              "urca", "uuid", "V8", "VarianceGamma", "vars", "vcd", "vcdExtra", 
              "Vdgraph", "vegan", "verification", "VGAM", "VGAMdata", "VIM", 
              "VIMGUI", "VineCopula", "vioplot", "viridis", "viridisLite", 
              "visNetwork", "vtreat", "wavelets", "waveslim", "wbstats", "webp", 
              "webshot", "WhatIf", "whisker", "whoami", "withr", "woe", 
              "wordcloud", "WrightMap", "WriteXLS", "wskm", "wsrf", 
              "xergm.common", "xkcd", "XML", 
              "xml2", "xtable", "xts", "YaleToolkit", "yaml", "yarrr", "zeallot", 
              "Zelig", "zip", "zipcode", "zoo", "ztable", "getPass", 
              "miniCRAN", "NMOF", "odbc", "recosystem", "rgeoapi", "smooth", "stR", 
              "Boom", "BoomSpikeSlab", "bsts", "CausalImpact", "cli", "ClusterR", 
              "emmeans", "FD", "fromo", "gdalUtils", "geojson", "geojsonio", 
              "geojsonlint", "geometry", "ggridges", "inum", "jqr", 
              "jsonvalidate", "libcoin", "magic", "manipulateWidget", "mapview", 
              "moments", "NADA", "OceanView", "OpenImageR", "osmar", "pillar", 
              "plot3Drgl", "protolite", "rmapshaper", 
              "satellite", "sf", "spData", "SQUAREM", "tiff", "tmaptools", 
              "udunits2", "units", "uroot", "utf8", "xfun", 
              "zCompositions")
install.packages(packages[which(!(packages %in% rownames(installed.packages())))], Ncpus = parallel::detectCores())
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
devtools::install_github("cmpolis/datacomb", subdir = "pkg", ref = "1.1.2")
```

Installing xgboost is specific to be done in bash:

```
git clone --recursive https://github.com/dmlc/xgboost
cd xgboost
git checkout 8f6aadd
cd R-package
```

In src/install.libs.R, add `-DUSE_AVX=ON` at line 11:

```
gedit src/Makevars.in
```

Then compile xgboost:

```
R
install.packages('.', repos = NULL, type = "source")
```

For LightGBM, this also requires a specific installation in bash:

```
git clone --recursive https://github.com/Microsoft/LightGBM
cd LightGBM
git checkout 3f54429
cd R-package
```

In src/Makevars.in, replace `cmake_cmd` content line 50 by `"cmake -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc "`:

```
gedit src/install.libs.R
```

Then compile LightGBM:

```
R CMD INSTALL --build . --no-multiarch
```

Install even more packages...:

```r
devtools::install_github("Laurae2/Laurae")
devtools::install_github("Laurae2/LauraeParallel")
devtools::install_github("Laurae2/LauraeDS")
devtools::install_github("Laurae2/LauraeCE")
install.packages("https://cran.r-project.org/src/contrib/Archive/tabplot/tabplot_1.1.tar.gz", repos=NULL, type="source") # Further versions are too bad / not reliable / generated unreadable plots
```

### Step 11: I want to restart from scratch!!!

Run the following as Administrator to fully wipe the Ubuntu installation:

```sh
wslconfig.exe /u Ubuntu
```

</p>
</details>

## Windows Client / Windows Server

<details><summary>:information_desk_person: CLICK THE ARROW TO REVEAL Windows Client / Server steps</summary>
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
devtools::install_github("cmpolis/datacomb", subdir = "pkg", ref = "1.1.2")
# devtools::install_github("ficonsulting/RInno", build_vignettes = TRUE) # If you want to build R executable standalones
```

Installining xgboost and LightGBM:

```r
xgbdl::xgb.dl(compiler = "Visual Studio 15 2017 Win64", commit = "8f6aadd", use_avx = TRUE, use_gpu = FALSE)
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

<details><summary>:information_desk_person: CLICK THE ARROW TO REVEAL Ubuntu 17.10 steps</summary>
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
devtools::install_github("cmpolis/datacomb", subdir = "pkg", ref = "1.1.2")
```

Installining xgboost and LightGBM:

```r
xgbdl::xgb.dl(compiler = "gcc", commit = "8f6aadd", use_avx = TRUE, use_gpu = FALSE)
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

### Compilation Flags

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

### Compiler Choice

You may want to update the compiler to a more recent version to benefit from newer optimizations, especially on newer CPUs. You may also choose another compiler for performance, such as Intel Compiler.

In addition, you can use a different BLAS to improve the performance of R. Intel MKL is probably the best you can get, but it is not "free".

Please have a read if you intend to use Intel Compilers. Intel Compilers vs gcc: https://sites.google.com/view/lauraepp/intel-compiler (tl;dr: do not waste all those hours trying it, the benefit is very small unless you are running a very CPU-limited machine such as a laptop)

I provide below the list of compiler flags with Intel Compilers to use:

```sh
export CC="icc"
export CXX="icpc"
export F77="ifort"
export FC="ifort"
export AR="xiar"
export LD="xild"
export CFLAGS="-O3 -ipo -qopenmp -xHost -fPIC"
export CXXFLAGS="-O3 -ipo -qopenmp -xHost -fPIC -wd308"
export FFLAGS="-O3 -ipo -qopenmp -xHost -fPIC"
export FCFLAGS="-O3 -ipo -qopenmp -xHost -fPIC"
export LDFLAGS="-qopenmp"
export MKL="-L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl"
```
