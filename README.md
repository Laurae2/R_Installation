# R Installation - Windows / Debian / Ubuntu Version

**Last tested : R 3.5.2, 2019/01/27 (27 Jan, 2019)**

R packages for installation, the Windows / Debian / Ubuntu version.

Works well for Windows. Works well for Linux (Debian/Ubuntu-like).

This document helps you install over 1,000 packages.

Validation on Windows Subsystem for Linux:

| Operating System | Success | R Version |
| --- | --- | --- |
| Ubuntu 18.04 | :heavy_check_mark: Pass! | :100: R 3.5.x, R 3.4.4 |
| Ubuntu 16.04 | :heavy_check_mark: Pass! | :100: R 3.5.x, R 3.4.4 |

Validation on Windows operating systems:

| Operating System | Success | R Version |
| --- | --- | --- |
| Windows 10 (1809) | :heavy_check_mark: Pass! | :100: R 3.5.x, R 3.4.4 |
| Windows 10 (1803) | :heavy_check_mark: Pass! | :100: R 3.5.x, R 3.4.4 |
| Windows 10 (1709) | :heavy_check_mark: Pass! | :100: R 3.5.x, R 3.4.4 |
| Windows 10 (1703) | :heavy_check_mark: Pass! | :100: R 3.5.x, R 3.4.4 |
| Windows 10 (1607) | :heavy_check_mark: Pass! | :100: R 3.5.x, R 3.4.4 |
| Windows 10 (1511) | :heavy_check_mark: Pass! | :100: R 3.5.x, R 3.4.4 |
| Windows 10 (1507) | :heavy_check_mark: Pass! | :100: R 3.5.x, R 3.4.4 |
| Windows Server 2016 | :heavy_check_mark: Pass! | :100: R 3.5.x, R 3.4.4 |
| Windows Server 2012 R2 | :heavy_check_mark: Pass! | :100: R 3.5.x, R 3.4.4 |
| Windows 8.1 | :heavy_check_mark: Pass! | :100: R 3.5.x, R 3.4.4 |
| Windows Server 2012 | :heavy_check_mark: Pass! | :100: R 3.5.x, R 3.4.4 |
| Windows 7 | :heavy_check_mark: Pass! | :100: R 3.5.x, R 3.4.4 |
| Windows Server 2008 R2 | :heavy_check_mark: Pass! | :100: R 3.5.x, R 3.4.4 |
| Windows Vista | :heavy_exclamation_mark: Not fully passing! | :100: R 3.5.x, R 3.4.4 |
| Windows XP | :heavy_exclamation_mark: Not fully passing! | :100: R 3.5.x, R 3.4.4 |

Validation on Ubuntu operating systems:

| Operating System | Success | R Version |
| --- | --- | --- |
| pop!_OS 18.10 | :heavy_check_mark: Pass! | :100: R 3.5.x |
| Ubuntu 18.10 | :heavy_check_mark: Pass! | :100: R 3.5.x |
| pop!_OS 18.04 | :heavy_check_mark: Pass! (Not public) | :100: R 3.5.x, R 3.4.4 |
| Ubuntu 18.04 | :heavy_check_mark: Pass! (Not public) | :100: R 3.5.x, R 3.4.4 |
| pop!_OS 17.10 | :heavy_check_mark: Pass! (Not public) | :100: R 3.5.x, R 3.4.4 |
| Ubuntu 17.10 | :heavy_check_mark: Pass! | :100: R 3.5.x, R 3.4.4 |
| Ubuntu 17.04 | :heavy_check_mark: Pass! | :100: R 3.5.x, R 3.4.4 |
| Ubuntu 16.10 | :heavy_check_mark: Pass! | :100: R 3.5.x, R 3.4.4 |
| Ubuntu 16.04 | :heavy_check_mark: Pass! | :100: R 3.5.x, R 3.4.4 |
| Ubuntu 15.10 | :interrobang: Unknown... | :trident: None yet! |
| Ubuntu 14.10 | :interrobang: Unknown... | :trident: None yet! |
| Ubuntu 14.04 | :interrobang: Unknown... | :trident: None yet! |

N.B: Ubuntu 18.10 and pop!_OS 18.10 requires a `gfortran` update. We provide a different way of installing R and Python for 18.10: Intel Distribution for Python, and R with Intel MKL. We also provide GPU steps for Ubuntu 18.10 and pop!_OS 18.10.

Validation on SUSE operating systems:

| Operating System | Success | R Version |
| --- | --- | --- |
| SUSE 12 | :heavy_check_mark: Pass! (Not public) | :100: R 3.5.x |

The custom parts to change for your own use case are found by searching `# Change this`.

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

To perform in Bash shell. Add Anaconda to PATH so you can do the Python-specific stuff early, we will remove it afterwards.

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

Put the following in the file, change the PATH of R-lib and get rid of Anaconda!:

```sh
export DISPLAY="localhost:0.0"
export JAVA_HOME=/usr/lib/jvm/java-9-openjdk-amd64
export PATH=$JAVA_HOME/bin:$PATH
export PATH="$PATH:/usr/local/lib/R"
export R_HOME="/usr/local/lib/R"
export R_LIBS_USER="/mnt/e/WSL/R-lib/R-3.5.0" # Change this
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
r-libs-user=/mnt/e/WSL/R-lib/R-3.5.0 # Change this
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

install.packages("BiocManager")
BiocManager::install()
BiocManager::install(c("graph", "RBGL"))

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
devtools::install_github("wrathematics/float", upgrade_dependencies = FALSE)
```

Modify `/usr/lib/R/etc/Makeconf` for maximum performance for xgboost / LightGBM (+0-10% speed): change -O2 to -O3 to CXXFLAGS

Installining xgboost and LightGBM:

```r
xgbdl::xgb.dl(compiler = "gcc", commit = "4fac987", use_avx = TRUE, use_gpu = FALSE)
lgbdl::lgb.dl(commit = "f9a1465", compiler = "gcc")
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

### Step 12: I want RStudio Server and Desktop Preview!!!

Run the following for RStudio Desktop Preview:

```sh
sudo apt-get install libclang-3.8-dev libclang-common-3.8-dev libclang-dev libclang1-3.8 libllvm3.8 libobjc-5-dev libobjc4
wget https://s3.amazonaws.com/rstudio-ide-build/desktop/trusty/amd64/rstudio-1.2.747-amd64.deb
sudo gdebi rstudio-1.2.747-amd64.deb
```

Run the following for RStudio Server Preview:

```sh
sudo apt-get install libclang-3.8-dev libclang-common-3.8-dev libclang-dev libclang1-3.8 libllvm3.8 libobjc-5-dev libobjc4
wget https://s3.amazonaws.com/rstudio-ide-build/server/trusty/amd64/rstudio-server-1.2.747-amd64.deb
sudo gdebi rstudio-server-1.2.747-amd64.deb
```

Do not forget to set R preferences: inside `/etc/rstudio/rsession.conf`, add the following:

```sh
r-libs-user=/mnt/e/WSL/R-lib/R-3.5.0 # Change this
session-timeout-minutes=0
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

To perform in Bash shell. Add Anaconda to PATH so you can do the Python-specific stuff early, we will remove it afterwards.

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

Put the following in the file, change the PATH of R-lib and get rid of Anaconda!:

```sh
export DISPLAY="localhost:0.0"
export JAVA_HOME=/usr/lib/jvm/java-9-openjdk-amd64
export PATH=$JAVA_HOME/bin:$PATH
export PATH="$PATH:/usr/local/lib/R"
export R_LIBS_USER="/mnt/e/WSL/R-lib/R-3.5.0" # Change this
export R_HOME="/usr/local/lib/R"
source /opt/intel/compilers_and_libraries_2018.2.199/linux/bin/compilervars.sh intel64 # Change this
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
r-libs-user=/mnt/e/WSL/R-lib/R-3.5.0 # Change this
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
.libPaths("/mnt/e/WSL/R-lib/R-3.5.0") # Change this
install.packages("devtools", dependencies = TRUE)
install.packages("tcltk2", dependencies = TRUE)

install.packages("BiocManager")
BiocManager::install()
BiocManager::install(c("graph", "RBGL"))

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
devtools::install_github("wrathematics/float", upgrade_dependencies = FALSE)
```

Installing xgboost is specific to be done in bash:

```
git clone --recursive https://github.com/dmlc/xgboost
cd xgboost
git checkout 4fac987
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
git checkout f9a1465
cd R-package
```

In src/Makevars.in, replace `cmake_cmd` content line 50 by `"cmake -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc "`:

```
gedit src/install.libs.R
```

Then compile LightGBM:

```
Rscript build_r.R
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

### Can I get rid of RStudio Server / Desktop?

Yes:

```sh
sudo apt-get remove rstudio-server
sudo apt-get remove rstudio
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

Windows can download and install Anaconda here: https://repo.anaconda.com/archive/Anaconda3-5.1.0-Windows-x86_64.exe (do not add to PATH)

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

install.packages("BiocManager")
BiocManager::install()
BiocManager::install(c("graph", "RBGL"))

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
devtools::install_github("wrathematics/float", upgrade_dependencies = FALSE)
# devtools::install_github("ficonsulting/RInno", build_vignettes = TRUE) # If you want to build R executable standalones
```

Installining xgboost and LightGBM:

```r
xgbdl::xgb.dl(compiler = "Visual Studio 15 2017 Win64", commit = "4fac987", use_avx = TRUE, use_gpu = FALSE)
lgbdl::lgb.dl(commit = "f9a1465", compiler = "vs")
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

## Ubuntu 18.10 / Pop_OS! 18.10

<details><summary>:information_desk_person: CLICK THE ARROW TO REVEAL Ubuntu 18.10 / Pop_OS! 18.10 steps</summary>
<p>

This setup is for the high performance mode of R and Python. This is not recommended as a daily development driver, but as a high performance remote productive setup.

We use the R 3.5.2 instead of R 3.5.1 build.

### Step 1: Setup SSH

```sh
sudo passwd root
passwd
sudo apt-get update
sudo apt-get install openssh-server fail2ban
sudo gedit /etc/ssh/sshd_config => change port to anything you want
sudo systemctl restart ssh
reboot
```

### Step 2: Unstall x2go and xfce4

```sh
sudo add-apt-repository ppa:x2go/stable
sudo apt-get update
sudo apt-get install x2goserver x2goserver-xsession screen
sudo apt-get install xfce4 xfce4-goodies
```

When x2go crashes (server-sided), you can use `sudo service x2goserver restart` to restart x2go.

### Step 3: Disable Network and Printer Discovery

Disabling printer and network discovery can be very useful.

```sh
sudo systemctl disable avahi-daemon
sudo systemctl stop avahi-daemon
sudo systemctl disable cups-browsed
sudo systemctl stop cups-browsed
```

You can do `ss -lntu` to know what ports are still in use

### Step 4: Get some monitoring

htop and iotop can help for troubleshooting issues...:

```sh
sudo apt-get install htop iotop
```

### Step 5: Performance mode of the processors

Setup performance mode to the processors:

```sh
sudo apt-get install cpufrequtils
echo 'GOVERNOR="performance"' | sudo tee /etc/default/cpufrequtils
sudo systemctl disable ondemand
```

### Step 6: Customize kernel parameters to get better performance

Install and launch `grub-customizer`:

```sh
sudo add-apt-repository ppa:danielrichter2007/grub-customizer
sudo apt-get update
sudo apt-get install grub-customizer
sudo grub-customizer
```

Disable most Meltdown and Spectre mitigations by adding the following to kernel parameters:

```sh
pti=off spectre_v2=off spec_store_bypass_disable=off l1tf=off noibrs noibpb nopti no_stf_barrier 
```

Edit `/etc/sysctl.conf` and add the following at the bottom to disable ASLR:

```sh
kernel.randomize_va_space=0
```

Reboot:

```sh
sudo reboot
```

Test your computer against vulnerabilities, two should be negative (red):

```sh
cd Desktop
wget https://raw.githubusercontent.com/speed47/spectre-meltdown-checker/master/spectre-meltdown-checker.sh
sudo sh spectre-meltdown-checker.sh
```

### Step 7: Setup Python (Intel Distribution for Python)

To setup Intel Distribution for Python, we are using Miniconda, which is a lightweight Anaconda (do yes, enter, and no for the Miniconda installation when requested):

```sh
cd Downloads
mkdir R
cd R
sudo apt-get install curl
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Enable symlinks to access `conda` from anywhere (:

```sh
sudo ln -s /home/laurae/miniconda3/bin/conda /usr/bin/conda  # Change this user
sudo ln -s /home/laurae/miniconda3/bin/activate /usr/bin/activate  # Change this user
sudo ln -s /home/laurae/miniconda3/bin/deactivate /usr/bin/deactivate  # Change this user
```

Then we can install Intel Distribution for Python from conda:

```
conda config --add channels intel
conda create -n r-tensorflow intelpython3_full
```

To install keras, do the following:

```sh
source activate r-tensorflow
conda install keras
source deactivate
```

Unhappy with the environment? Uninstall it:

```sh
conda remove --name r-tensorflow --all
```

Unhappy with Intel Python? Remove the channel:

```sh
conda config --remove channels intel
```

Want to destroy Conda fully? Open a console and run the following:

```sh
rm -rf miniconda
rm -rf ./condarc
```

### Step 7 Alternative: Setup Python (Conda)

To setup Conda for Python, we are using Miniconda, which is a lightweight Anaconda (do yes, enter, and no for the Miniconda installation when requested):

```sh
cd Downloads
mkdir R
cd R
sudo apt-get install curl
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Enable symlinks to access `conda` from anywhere (:

```sh
sudo ln -s /home/laurae/miniconda3/bin/conda /usr/bin/conda  # Change this user
sudo ln -s /home/laurae/miniconda3/bin/activate /usr/bin/activate  # Change this user
sudo ln -s /home/laurae/miniconda3/bin/deactivate /usr/bin/deactivate  # Change this user
```

Then we can install Tensorflow and all Anaconda packages from conda:

```
conda create -n r-tensorflow anaconda tensorflow keras spyder scikit-learn=0.20
```

To install keras, do the following:

```sh
source activate r-tensorflow
conda install keras
source deactivate
```

Unhappy with the environment? Uninstall it:

```sh
conda remove --name r-tensorflow --all
```

Want to destroy Conda fully? Open a console and run the following:

```sh
rm -rf miniconda
rm -rf ./condarc
```

### Step 8: Install R pre-requisites

Lot of dependencies...:

```sh
sudo add-apt-repository -y ppa:opencpu/jq
sudo add-apt-repository "deb http://ppa.launchpad.net/ubuntugis/ppa/ubuntu $(lsb_release -sc) main"
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 314DF160
sudo apt-get update
sudo apt-get install cmake curl libxml2-dev git gfortran libssl-dev libssh2-1-dev gdebi-core libwebp-dev libprotobuf-dev libjq-dev libpoppler-cpp-dev libcairo2-dev librsvg2-dev libv8-3.14-dev libgtk2.0-dev default-jre default-jdk libgmp3-dev libgsl-dev jags libudunits2-dev protobuf-compiler mesa-common-dev libglu1-mesa-dev coinor-libsymphony-dev libtiff5-dev tcl-dev tk-dev libmpfr-dev ggobi libgdal-dev libglpk-dev libgeos-dev libgeos++-dev netcdf-bin libfftw3-dev libopenmpi-dev bwidget mpi-default-bin libcurl4-openssl-dev libreadline-dev libbz2-dev cmake libxml2-dev git-core libssl-dev libssh2-1-dev gdebi-core libwebp-dev libprotobuf-dev libpoppler-cpp-dev libcairo2-dev librsvg2-dev libv8-3.14-dev libgtk2.0-dev default-jre default-jdk libgmp3-dev libgsl-dev jags libudunits2-dev protobuf-compiler mesa-common-dev libglu1-mesa-dev coinor-libsymphony-dev libtiff5-dev tcl-dev tk-dev libmpfr-dev ggobi libgdal-dev libglpk-dev libgeos-dev netcdf-bin libfftw3-dev libopenmpi-dev bwidget mpi-default-bin libx11-dev ratfor libproj-dev libmagick++-dev coinor-libsymphony-dev coinor-libcgl-dev libjq-dev libgdal-dev texlive-latex-base texlive-fonts-extra
```

### Step 9: Intel MKL

Add Intel MKL to your Linux installation for speed:

```sh
wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
sudo apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
sudo wget https://apt.repos.intel.com/setup/intelproducts.list -O /etc/apt/sources.list.d/intelproducts.list
sudo apt-get update
sudo apt-get install intel-mkl-64bit-2019.0-045
```

### Step 10: Download, Configure, and Compile R

If you followed exactly the steps, it should work out of the box for maximum performance R:

```sh
wget https://cran.r-project.org/src/base/R-3/R-3.5.2.tar.gz
tar zxvf R-3.5.2.tar.gz
cd R-3.5.2
export CFLAGS="-O3 -g -mtune=native -fno-stack-protector"
export CXXFLAGS="-O3 -g -mtune=native -fno-stack-protector"
export FFLAGS="-O3 -g -mtune=native -fno-stack-protector"
export FCFLAGS="-O3 -g -mtune=native -fno-stack-protector"
source /opt/intel/mkl/bin/mklvars.sh intel64 mod lp64
MKL=" -L$MKLROOT/lib/intel64 -Wl,--start-group -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lm -Wl,--end-group"
./configure --enable-R-shlib --with-blas="$MKL" --with-lapack
make -j 36 # Change this core count
sudo make install
```

From now, make sure you always run R with sudo unless you want shared libraries... RStudio Server always does but RStudio Desktop requires to be run with sudo so you were warned!

### Step 11: Setup ~/.bashrc

You can change ~/.bashrc:

```sh
sudo gedit ~/.bashrc
```

Add the following:

```sh
export PATH="$PATH:/usr/local/lib/R"
export R_HOME="/usr/local/lib/R"
export R_LIBS_USER="usr/local/lib/R/library"
source /opt/intel/mkl/bin/mklvars.sh intel64 mod lp64
```

### Step 12: Install a lot of R packages

Run an R console using `sudo R` (or use `R` if you want to separate user libraries).

Install a crazy lot of packages from R afterwards...:

```r
install.packages("devtools", dependencies = TRUE, Ncpus = parallel::detectCores())
install.packages("tcltk2", dependencies = TRUE, Ncpus = parallel::detectCores())

install.packages("BiocManager")
BiocManager::install()
BiocManager::install(c("graph", "RBGL"))

packages <- c("abind", "acepack", "actuar", "ada", "adabag", "ade4", "ade4TkGUI")
if (length(packages[which(!(packages %in% rownames(installed.packages())))]) > 0) {install.packages(packages[which(!(packages %in% rownames(installed.packages())))], dependencies = TRUE, Ncpus = parallel::detectCores())}
if (length(packages[which(!(packages %in% rownames(installed.packages())))]) > 0) {install.packages(packages[which(!(packages %in% rownames(installed.packages())))], dependencies = TRUE, Ncpus = parallel::detectCores())}

packages <- c("adegraphics", "adehabitatLT", "adehabitatMA", "ADGofTest", "AER", 
              "AGD", "agricolae", "AICcmodavg", "akima", "alabama", "AlgDesign", 
              "alphahull", "alr3", "alr4", "amap", "Amelia", "anchors", "animation", 
              "aod", "aods3", "ape", "aplpack", "argparse", "arm", "arules")
if (length(packages[which(!(packages %in% rownames(installed.packages())))]) > 0) {install.packages(packages[which(!(packages %in% rownames(installed.packages())))], dependencies = TRUE, Ncpus = parallel::detectCores())}
if (length(packages[which(!(packages %in% rownames(installed.packages())))]) > 0) {install.packages(packages[which(!(packages %in% rownames(installed.packages())))], dependencies = TRUE, Ncpus = parallel::detectCores())}

packages <- c("arulesViz", "ascii", "assertthat", "AUC", "BaBooN", "backports", 
              "barcode", "bartMachine", "bartMachineJARs", "base64", "base64enc", 
              "BatchJobs", "BayesFactor", "bayesplot", "BayesX", "BayesXsrc", 
              "BB", "BBmisc", "bbmle", "BCA", "bcp", "BDgraph", "bdsmatrix", 
              "betareg", "BH", "BHH2", "BiasedUrn", "bibtex", "biclust", "biganalytics", 
              "biglm", "bigmemory", "bigmemory.sri", "bigtabulate", 
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
              "clustvarsel", "clv", "clValid", "cmaes", "cmprsk", 
              "coda", "codetools", "coin", "colorspace", "colourpicker", 
              "combinat", "commonmark", "CompQuadForm", "compute.es", "conf.design", 
              "config", "contfrac", "contrast", "copula", "CORElearn", "corpcor", 
              "corrgram", "corrplot", "covr", "cowplot", "CoxBoost", "coxme", 
              "cplm", "crayon", "crosstalk", "crossval", "crp.CSFP", "crrstep", 
              "crs", "cshapes", "cubature", "Cubist", "curl", "cvAUC", "CVST", 
              "cvTools", "d3heatmap", "d3Network", "DAAG", "dagitty", "data.table")
if (length(packages[which(!(packages %in% rownames(installed.packages())))]) > 0) {install.packages(packages[which(!(packages %in% rownames(installed.packages())))], dependencies = TRUE, Ncpus = parallel::detectCores())}
if (length(packages[which(!(packages %in% rownames(installed.packages())))]) > 0) {install.packages(packages[which(!(packages %in% rownames(installed.packages())))], dependencies = TRUE, Ncpus = parallel::detectCores())}

packages <- c("data.tree", "dataframes2xls", "date", "dbConnect", 
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
              "emdbook", "emoa", "emulator", "energy", "ENmisc", "entropy", 
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
              "fortunes", "forward", "fpc", "fPortfolio", "fracdiff", 
              "frbs", "fRegression", "FrF2", "FrF2.catlg128", "FSelector", 
              "fst", "fTrading", "fts", "futile.logger", "futile.options", 
              "future", "GA", "gam", "gamair", "GAMBoost", "gamboostLSS", "gamlss", 
              "gamlss.data", "gamlss.dist", "gamm4", "gapminder", "gbm", "gclus", 
              "gdata", "gdtools", "gee", "geeM", "geepack", 
              "GeneralizedHyperbolic", "genetics", "GenSA", "geoR", "geoRglm", 
              "geosphere", "GERGM", "getopt", "GGally", "ggcorrplot", "ggdendro", 
              "ggeffects", "ggExtra", "ggformula", "ggfortify", "ggiraph", 
              "ggm", "ggplot2", "ggplot2movies", "ggpubr", "ggrepel", "ggsci", 
              "ggsignif", "ggThemeAssist", "ggthemes", "ggvis", "giphyr", "git2r", 
              "gitgadget", "glasso", "glmmML", "glmmTMB", "glmnet", "glmulti", 
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
if (length(packages[which(!(packages %in% rownames(installed.packages())))]) > 0) {install.packages(packages[which(!(packages %in% rownames(installed.packages())))], dependencies = TRUE, Ncpus = parallel::detectCores())}
if (length(packages[which(!(packages %in% rownames(installed.packages())))]) > 0) {install.packages(packages[which(!(packages %in% rownames(installed.packages())))], dependencies = TRUE, Ncpus = parallel::detectCores())}

packages <- c("ic.infer", "ICS", "ICSNP", "igraph", "igraphdata", "import", 
              "imputeTS", "ineq", "influenceR", "Information", "infotheo", 
              "inline", "inlinedocs", "intergraph", "intervals", "intsvy", 
              "iplots", "ipred", "irace", "irlba", "irr", "isa2", "Iso", "ISOcodes", 
              "isotone", "ISwR", "iterators", "itertools", "JavaGD", "JGR", 
              "jomo", "jpeg", "jsonlite", "kappalab", "kdecopula", "Kendall", 
              "keras", "kernlab", "KernSmooth", "KFAS", "kinship2", "kknn", 
              "klaR", "kmi", "knitcitations", "knitr", "kohonen", "koRpus", 
              "ks", "labeling", "labelled", "laeken", "LaF", "laGP", "Lahman", 
              "lambda.r", "lars", "lasso2", "latentnet", "lattice", 
              "latticeExtra", "lava", "lava.tobit", "lavaan", "lavaan.survey", 
              "lawstat", "lazyeval", "LCA", "lcopula", "leaflet", "leaps", 
              "LearnBayes", "lfda", "lfe", "lhs", "LiblineaR", "likert", "linprog", 
              "lintr", "lisrelToR", "listenv", "lme4", "lmerTest", 
              "lmodel2", "lmtest", "loa", "locfit", "logcondens", "LogicReg", 
              "logistf", "logspline", "lokern", "longmemo", "loo", "lpSolve", 
              "lpSolveAPI", "lqmm", "lsmeans", "lubridate", "MAc", "MAd", 
              "magrittr", "mail", "manipulate", "mapdata", "mapproj", "maps", 
              "maptools", "maptree", "markdown", "MASS", "Matching", "MatchIt", 
              "matlab", "Matrix", "matrixcalc", "MatrixModels", 
              "matrixStats", "maxLik", "maxlike", "MBA", "MBESS", "mboost", 
              "mc2d", "mcclust", "mcgibbsit", "mclust", "mcmc", "MCMCglmm", 
              "MCMCpack", "mco", "mda", "mediation", "memisc", "memoise", 
              "MEMSS", "merTools", "MetABEL", "metafor", "Metrics", "mets", 
              "mgcv", "mi", "mice", "miceadds", "microbenchmark", "microplot", 
              "mime", "minerva", "miniUI", "minpack.lm", "minqa", "mirt", "mirtCAT", 
              "misc3d", "miscTools", "missForest", "missMDA", "mitml", "mitools", 
              "mix", "mlbench", "MLmetrics", "mlmRev", "mlogit", "mlr", "mlrMBO", 
              "mnlogit", "mnormt", "ModelMetrics", "modelr", "modeltools", 
              "mondate", "monreg", "moonBook", "mosaic", "mosaicCalc", "mosaicCore", 
              "mosaicData", "movMF", "MplusAutomation", "mpmi", "MPV", "mratios", 
              "mRMRe", "msm", "mstate", "MSwM", "muhaz", "multcomp", "multcompView", 
              "multicool", "multiwayvcov", "MuMIn", "munsell", "mvinfluence", 
              "mvnormtest", "mvoutlier", "mvtnorm", "NbClust", "ncdf4", "ncvreg", 
              "ndtv", "network", "networkD3", "networkDynamic", "networkDynamicData", 
              "networksis", "neuralnet", "NeuralNetTools", "NHANES", "nlme", 
              "nloptr", "NLP", "NMF", "nnet", "nnls", "nodeHarvest", "nor1mix", 
              "norm", "nortest", "np", "numbers", "numDeriv", "nws", "nycflights13", 
              "obliqueRF", "officer", "OpenMx", "openssl", "openxlsx", 
              "optextras", "optimx", "optmatch", "orcutt", "ordinal", "ore", 
              "orloca", "orloca.es", "orthopolynom", "outliers", "oz", "packrat")
if (length(packages[which(!(packages %in% rownames(installed.packages())))]) > 0) {install.packages(packages[which(!(packages %in% rownames(installed.packages())))], dependencies = TRUE, Ncpus = parallel::detectCores())}
if (length(packages[which(!(packages %in% rownames(installed.packages())))]) > 0) {install.packages(packages[which(!(packages %in% rownames(installed.packages())))], dependencies = TRUE, Ncpus = parallel::detectCores())}

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
              "R.oo", "R.rsp", "R.utils", "R2BayesX", "R2HTML", "R2jags", 
              "R2OpenBUGS", "R2WinBUGS", "R6", "radiant", 
              "radiant.basics", "radiant.data", "radiant.design", "radiant.model", 
              "radiant.multivariate", "RandomFields", "RandomFieldsUtils", 
              "randomForest", "randomForestSRC", "randtests", "randtoolbox", 
              "ranger", "RankAggreg", "RANN", "rappdirs", "RArcInfo", "rARPACK", 
              "RaschSampler", "raster", "rasterVis", "rattle", "rbenchmark", 
              "rbounds", "rbvs", "Rcgmin", "Rcmdr", "RcmdrMisc", 
              "RcmdrPlugin.coin", "RcmdrPlugin.depthTools", "RcmdrPlugin.DoE", 
              "RcmdrPlugin.Export", 
              "RcmdrPlugin.FactoMineR", "RcmdrPlugin.HH", "RcmdrPlugin.IPSUR", 
              "RcmdrPlugin.KMggplot2", "RcmdrPlugin.mosaic", "RcmdrPlugin.orloca", 
              "RcmdrPlugin.pointG", "RcmdrPlugin.qual", "RcmdrPlugin.SLC", 
              "RcmdrPlugin.sos", "RcmdrPlugin.steepness", "RcmdrPlugin.survival", 
              "RcmdrPlugin.TeachingDemos", "RcmdrPlugin.UCA", "RColorBrewer", 
              "Rcpp", "RcppArmadillo", "RcppCNPy", "RcppDE", "RcppEigen", "RcppParallel", 
              "RcppProgress", "RcppRoll", "Rcsdp", "RCurl", "readr", "readstata13", 
              "readxl", "recipes", "recommenderlab", "RefManageR", "registry", 
              "relaimpo", "relations", "relevent", "reliaR", "relimp", 
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
              "sourcetools", "sp", "spacetime", "spam", "SparseGrid", 
              "sparseLDA", "SparseM", "sparsio", "spatial", "spatstat", "spatstat.data", 
              "spatstat.utils", "spc", "spd", "spdep", "speedglm", "sphet", 
              "splancs", "splm", "spls", "sqldf", "sROC", "stabledist", "stabs", 
              "StanHeaders", "startupmsg", "StatMatch", "statmod", "statnet", 
              "statnet.common", "steepness", "stepPlr", "stinepack", "stringdist", 
              "stringi", "stringr", "strucchange", "subselect", "subsemble", 
              "sudoku", "SuperLearner", "superpc", "SuppDists", "survey", "survival", 
              "svd", "svglite", "svGUI", "svUnit", "svyPVpack", "SwarmSVM", 
              "systemfit", "tables", "tabplot", "tabplotd3")
if (length(packages[which(!(packages %in% rownames(installed.packages())))]) > 0) {install.packages(packages[which(!(packages %in% rownames(installed.packages())))], dependencies = TRUE, Ncpus = parallel::detectCores())}
if (length(packages[which(!(packages %in% rownames(installed.packages())))]) > 0) {install.packages(packages[which(!(packages %in% rownames(installed.packages())))], dependencies = TRUE, Ncpus = parallel::detectCores())}

packages <- c("TAM", "tclust", "TeachingDemos", "tensor", "tensorA", 
              "tensorflow", "tergm", "testit", "testthat", "texreg", "tfestimators", 
              "tfruns", "tgp", "TH.data", "threejs", "tibble", "tidyr", "tidyselect", 
              "tikzDevice", "timeDate", "timereg", "timeSeries", "tis", "tkrplot", 
              "tm", "tmap", "TMB", "tmvtnorm", "tnam", "tree", 
              "trimcluster", "tripack", "truncdist", "truncnorm", "truncreg", 
              "trust", "TSA", "tseries", "tsna", "TSP", "TTR", 
              "tufte", "tuneR", "tweedie", "ucminf", "uniReg", "unmarked", 
              "urca", "uuid", "V8", "VarianceGamma", "vars", "vcd", "vcdExtra", 
              "Vdgraph", "vegan", "verification", "VGAM", "VGAMdata", "VIM", 
              "VIMGUI", "VineCopula", "vioplot", "viridis", "viridisLite", 
              "visNetwork", "vtreat", "wavelets", "waveslim", "wbstats", "webp", 
              "webshot", "WhatIf", "whisker", "whoami", "withr", "woe", 
              "wordcloud", "WrightMap", "WriteXLS", "wskm", "wsrf", "xergm", 
              "xergm.common", "xkcd", "XLConnect", "XLConnectJars", "XML", 
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
if (length(packages[which(!(packages %in% rownames(installed.packages())))]) > 0) {install.packages(packages[which(!(packages %in% rownames(installed.packages())))], dependencies = TRUE, Ncpus = parallel::detectCores())}
if (length(packages[which(!(packages %in% rownames(installed.packages())))]) > 0) {install.packages(packages[which(!(packages %in% rownames(installed.packages())))], dependencies = TRUE, Ncpus = parallel::detectCores())}

devtools::install_github("r-lib/progress@54a14f5", upgrade_dependencies = FALSE)
devtools::install_github("elbamos/largeVis", upgrade_dependencies = FALSE)
devtools::install_github("Laurae2/woe", upgrade_dependencies = FALSE)
devtools::install_github("Laurae2/xgbdl", upgrade_dependencies = FALSE)
devtools::install_github("Laurae2/lgbdl", upgrade_dependencies = FALSE)
devtools::install_github("twitter/AnomalyDetection", upgrade_dependencies = FALSE)
devtools::install_github("cmpolis/datacomb", subdir = "pkg", ref = "1.1.2", upgrade_dependencies = FALSE)
devtools::install_github("wrathematics/float", upgrade_dependencies = FALSE)

xgbdl::xgb.dl(compiler = "gcc", commit = "4fac987", use_avx = FALSE, use_gpu = FALSE)
lgbdl::lgb.dl(commit = "f9a1465", compiler = "gcc", cores = 1)

devtools::install_github("Laurae2/Laurae", upgrade_dependencies = FALSE)
devtools::install_github("Laurae2/LauraeParallel", upgrade_dependencies = FALSE)
devtools::install_github("Laurae2/LauraeDS", upgrade_dependencies = FALSE)
devtools::install_github("Laurae2/LauraeCE", upgrade_dependencies = FALSE)
install.packages("https://cran.r-project.org/src/contrib/Archive/tabplot/tabplot_1.1.tar.gz", repos=NULL, type="source") # Further versions are too bad / not reliable / generated unreadable plots
```

### Step 13: Can't believe I'm running Intel MKL in R

You can run the following in R:

```r
sessionInfo()
```

This should give you MKL:

```r
R version 3.5.2 (2018-12-20) -- "Eggshell Igloo"
Copyright (C) 2018 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu (64-bit)

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

  Natural language support but running in an English locale

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

> sessionInfo()
R version 3.5.2 (2018-12-20)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Pop!_OS 18.10

Matrix products: default
BLAS/LAPACK: /opt/intel/compilers_and_libraries_2019.0.117/linux/mkl/lib/intel64_lin/libmkl_gf_lp64.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=fr_FR.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=fr_FR.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=fr_FR.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=fr_FR.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

loaded via a namespace (and not attached):
[1] compiler_3.5.2
```

You can also run the following in your `Download/R/R-3.5.2` folder:

```sh
ldd ./lib/libR.so
```

You should get the following mentionning MKL:

```sh
laurae@laurae:~/Downloads/R/R-3.5.2$ ldd ./lib/libR.so
	linux-vdso.so.1 (0x00007fff67559000)
	libmkl_gf_lp64.so => /opt/intel/compilers_and_libraries_2019.0.117/linux/mkl/lib/intel64_lin/libmkl_gf_lp64.so (0x00007fb28ffa3000)
	libmkl_sequential.so => /opt/intel/compilers_and_libraries_2019.0.117/linux/mkl/lib/intel64_lin/libmkl_sequential.so (0x00007fb28ea0b000)
	libmkl_core.so => /opt/intel/compilers_and_libraries_2019.0.117/linux/mkl/lib/intel64_lin/libmkl_core.so (0x00007fb28a8d2000)
	libm.so.6 => /lib/x86_64-linux-gnu/libm.so.6 (0x00007fb28a71e000)
	libreadline.so.7 => /lib/x86_64-linux-gnu/libreadline.so.7 (0x00007fb28a4d5000)
	libpcre.so.3 => /lib/x86_64-linux-gnu/libpcre.so.3 (0x00007fb28a45f000)
	liblzma.so.5 => /lib/x86_64-linux-gnu/liblzma.so.5 (0x00007fb28a239000)
	libbz2.so.1.0 => /lib/x86_64-linux-gnu/libbz2.so.1.0 (0x00007fb28a226000)
	libz.so.1 => /lib/x86_64-linux-gnu/libz.so.1 (0x00007fb28a009000)
	librt.so.1 => /lib/x86_64-linux-gnu/librt.so.1 (0x00007fb289fff000)
	libdl.so.2 => /lib/x86_64-linux-gnu/libdl.so.2 (0x00007fb289ff9000)
	libicuuc.so.60 => /usr/lib/x86_64-linux-gnu/libicuuc.so.60 (0x00007fb289c40000)
	libicui18n.so.60 => /usr/lib/x86_64-linux-gnu/libicui18n.so.60 (0x00007fb28979f000)
	libgomp.so.1 => /usr/lib/x86_64-linux-gnu/libgomp.so.1 (0x00007fb28976e000)
	libpthread.so.0 => /lib/x86_64-linux-gnu/libpthread.so.0 (0x00007fb28974d000)
	libc.so.6 => /lib/x86_64-linux-gnu/libc.so.6 (0x00007fb289563000)
	/lib64/ld-linux-x86-64.so.2 (0x00007fb290f91000)
	libtinfo.so.6 => /lib/x86_64-linux-gnu/libtinfo.so.6 (0x00007fb289338000)
	libicudata.so.60 => /usr/lib/x86_64-linux-gnu/libicudata.so.60 (0x00007fb28778d000)
	libstdc++.so.6 => /usr/lib/x86_64-linux-gnu/libstdc++.so.6 (0x00007fb287603000)
	libgcc_s.so.1 => /lib/x86_64-linux-gnu/libgcc_s.so.1 (0x00007fb2875e9000)

```

Still not sure you are running Intel MKL? Placebo effect? Run the following:

```r
X <- matrix(rnorm(1000000), 1000, 1000)
system.time(svd(X))
```

If you are getting under 1 second, you are not using the standard R BLAS (because you need a computer with over 7 GHz to go under 1 second with the standard R BLAS). Here is MKL (OpenBLAS fares similarly):

```r
> X <- matrix(rnorm(1000000), 1000, 1000)
> system.time(svd(X))
   user  system elapsed 
  0.410   0.024   0.434 
```

With standard R BLAS, speed is doomed to favor (near) perfect reproducibility:

```r
> X <- matrix(rnorm(1000000), 1000, 1000)
> system.time(svd(X))
   user  system elapsed 
   2.80    0.02    2.84 
```

### Step 14: Install RStudio

Run the following for RStudio Desktop Preview:

```sh
sudo apt-get install lib32gcc1 lib32stdc++6 libc6-i386 libclang-7-dev libclang-common-7-dev libclang-dev libclang1-7 libgc1c2 libobjc-8-dev libobjc4
wget https://s3.amazonaws.com/rstudio-ide-build/desktop/trusty/amd64/rstudio-1.2.1070-amd64.deb
sudo gdebi rstudio-1.2.1070-amd64.deb
```

Run the following for RStudio Server Preview:

```sh
sudo apt-get install lib32gcc1 lib32stdc++6 libc6-i386 libclang-7-dev libclang-common-7-dev libclang-dev libclang1-7 libgc1c2 libobjc-8-dev libobjc4
wget https://s3.amazonaws.com/rstudio-ide-build/server/trusty/amd64/rstudio-server-1.2.1070-amd64.deb
sudo gdebi rstudio-server-1.2.1070-amd64.deb
```

Do not forget to set R preferences (`sudo nano /etc/rstudio/rsession.conf`): inside `/etc/rstudio/rsession.conf`, add the following:

```sh
r-libs-user=usr/local/lib/R/library
session-timeout-minutes=0
```

And also inside `/etc/rstudio/rserver.conf` to protect RStudio Server (`sudo nano /etc/rstudio/rserver.conf`) if required for front-facing Internet servers:

```sh
www-address=127.0.0.1
```

Then `sudo reboot` your server to make all the changes stick forever. Otherwise, they are not applied immediately.

Now, you can use RStudio Server from your own desktop using a SSH tunnel and port forwarding to 8787.

### Step 15: Unhappy with R or Intel MKL

Uninstall R and keep libraries, by running the following from Downloads/R/R-3.5.2 (where you ran `sudo make install`):

```
sudo make uninstall
```

Get rid definitely of all packages and of the R folder if you want:

```sh
sudo rm -rf usr/local/lib/R
```

### Step 16: Best VNC for remote control

Unhappy with x2go crashing on Windows? Use VNC! It works very well even for very high latency (>300ms) low bandwidth (<1Mbps) networks/routing through Internet. Setup with SSH port forwarding obviously to not give out your passwords in clear text!

TigerVNC, this is the BEST solution:

```sh
sudo apt-get install tigervnc-standalone-server
```

Using `nano ~/.vnc/xstartup` (Ctrl+X + y + Enter to exit once done), insert the following:

```sh
#!/bin/bash
xrdb $HOME/.Xresources
autocutsel -fork
startxfce4 &
```

Allow rights on the newly created file...:

```
sudo chmod +x ~/.vnc/xstartup
```

Edit VNC configuration file using `sudo nano /etc/vnc.conf` and add:

```
$localhost = "yes";
$depth = "16";
```

Now you can run `vncserver` from your user, and connect to VNC using a SSH tunnel with port forwarding to port 5901!

### Step 17: Add GPU (NVIDIA) support (pop_OS! only)

We will be using System76's repositories to install GPU support for NVIDIA.

#### Step 17.1: General Setup for NVIDIA and OpenCL

Run the following:

```sh
sudo apt-add-repository -y ppa:system76-dev/stable
sudo apt-get update
sudo apt-get install -y system76-driver
sudo apt-get install system76-driver-nvidia
sudo apt-get install system76-cuda-10.0
sudo apt-get install system76-cudnn-10.0
sudo apt-get install g++-6
sudo apt-get install ocl-icd-opencl-dev
sudo apt-get install opencl-headers clinfo
sudo apt-get install nvidia-opencl-icd
sudo apt-get install libnvidia-compute-410
sudo apt-get install gcc-8-offload-nvptx
```

It allows to:

* Install NVIDIA drivers + CUDA + CuDNN
* Install the gcc required by NVIDIA (so you do not need to hack `$CUDA_HOME/include/host_config.h`), you can use by calling `gcc-6` instead of `gcc`
* OpenCL headers
* gcc GPU offloading for NVIDIA GPUs (might crash with very strange error... it would be very rare you need it anyways)

Add the following to ~/.bashrc using `sudo nano ~/.bashrc`:

```sh
CUDA_HOME="/usr/lib/cuda"
export PATH="$PATH:/usr/lib/cuda/bin:/usr/lib/x86_64-linux-gnu"
```

If you are using multiple NVIDIA GPUs, download NCCL: https://developer.nvidia.com/nccl - we will assume it is nccl-repo-ubuntu1804-2.3.7-ga-cuda10.0_1-1_amd64.deb. Run the following:

```sh
sudo dpkg -i nccl-repo-ubuntu1804-2.3.7-ga-cuda10.0_1-1_amd64.deb
```

It will tell you the command to run to add the key for the repository:

```sh
sudo apt-key add /var/nccl-repo-2.3.7-ga-cuda10.0/7fa2af80.pub
```

Run the following once done:

```sh
sudo dpkg -i nccl-repo-ubuntu1804-2.3.7-ga-cuda10.0_1-1_amd64.deb
sudo apt-get update
sudo apt-get install libnccl2 libnccl-dev
```

#### Step 17.2: Add Monitoring

Use the following to add monitoring for GPUs:

```sh
mkdir nvtop
cd nvtop
git clone https://github.com/Syllo/nvtop.git
mkdir -p nvtop/build
cd nvtop/build
cmake .. -DNVML_RETRIEVE_HEADER_ONLINE=True -DCMAKE_BUILD_TYPE=Release
make
sudo make install
```

You can run the monitoring using `nvtop`.

#### Step 17.3: Add GPU support for some R packages

In R, run the following to install gpuR, xgboost (without NCCL, check below with NCCL), and LightGBM with GPU support:

```r
devtools::install_github("cdeterman/gpuR@cuda")
xgbdl::xgb.dl(compiler = "gcc", commit = "4fac987", use_avx = FALSE, use_gpu = TRUE, CUDA = list("/usr/lib/cuda", "/usr/bin/gcc-6", "/usr/bin/g++-6"))
lgbdl::lgb.dl(commit = "f9a1465", compiler = "gcc", use_gpu = TRUE)
```

If you are using multiple GPUs, you can install xgboost with NCCL to support multiple GPUs:

```r
xgbdl::xgb.dl(compiler = "gcc", commit = "4fac987", use_avx = FALSE, use_gpu = TRUE, CUDA = list("/usr/lib/cuda", "/usr/bin/gcc-6", "/usr/bin/g++-6"), NCCL = "/usr/lib/x86_64-linux-gnu")
```

#### Step 17.4: Add GPU support for some Python packages

In a command line, run the following to add Tensorflow with GPU support:

```py
conda install tensorflow-gpu
```

#### Step 17.5: Test GPU support for CLI (CUDA + OpenCL)

In a command line:

```sh
nvidia-smi
nvcc --version
clinfo
```

If you are unlucky and it does not work, run the following to get rid of the open source Nouveau drivers in Ubuntu:

```
sudo bash -c "echo blacklist nouveau > /etc/modprobe.d/blacklist-nvidia-nouveau.conf"
sudo bash -c "echo options nouveau modeset=0 >> /etc/modprobe.d/blacklist-nvidia-nouveau.conf"
cat /etc/modprobe.d/blacklist-nvidia-nouveau.conf
sudo update-initramfs -u
sudo reboot
```

#### Step 17.6: Test gpuR

Run from R:

```r
library(gpuR)

detectGPUs()
listContexts()

set.seed(11111)
gpuA <- gpuMatrix(rnorm(262144), nrow = 512, ncol = 512)
gpuB <- gpuA %*% gpuA
as.numeric(gpuB[1, 1]) == (32.51897782759634 + 7.105427e-15)
```

You should get `TRUE` at the end. If not, you should adjust the approximation for the comparison (use `all.equal`?).

#### Step 17.7: Test GPU Tensorflow from R

Run from R:

```r
library(reticulate)
use_condaenv(condaenv = "r-tensorflow", required = TRUE)

library(tensorflow)
sess <- tf$Session()
hello <- tf$constant('Hello, TensorFlow!')
sess$run(hello)
```

#### Step 17.8: Test GPU xgboost from R

Run from R, requires 4GB:

```r
library(xgboost)

set.seed(1)
N <- 500000
p <- 100
pp <- 25
X <- matrix(runif(N * p), ncol = p)
betas <- 2 * runif(pp) - 1
sel <- sort(sample(p, pp))
m <- X[, sel] %*% betas - 1 + rnorm(N)
y <- rbinom(N, 1, plogis(m))

tr <- sample.int(N, N * 0.90)

trainer <- function(n_cpus, n_gpus, n_iterations) {
  
  dtrain <- xgb.DMatrix(X[tr,], label = y[tr])
  dtest <- xgb.DMatrix(X[-tr,], label = y[-tr])
  wl <- list(test = dtest)
  
  if (n_gpus == 0) {
    
    pt <- proc.time()
    model <- xgb.train(list(objective = "reg:logistic", eval_metric = "logloss", subsample = 0.8, nthread = n_cpus, eta = 0.10,
                            max_bin = 64, tree_method = "hist"),
                       dtrain, watchlist = wl, nrounds = n_iterations)
    my_time <- proc.time() - pt
    
  } else {
    
    pt <- proc.time()
    model <- xgb.train(list(objective = "reg:logistic", eval_metric = "logloss", subsample = 0.8, nthread = n_cpus, eta = 0.10,
                            max_bin = 64, tree_method = "gpu_hist", n_gpus = n_gpus),
                       dtrain, watchlist = wl, nrounds = n_iterations)
    my_time <- proc.time() - pt
    
  }
  
  rm(model, dtrain, dtest)
  gc(verbose = FALSE)
  
  return(my_time)
  
}

trainer(1, 1, 50)
trainer(1, 0, 50)
```

#### Step 17.9: Test GPU LightGBM from R

Run from R:

```r
library(lightgbm)
data(agaricus.train, package = "lightgbm")
train <- agaricus.train
train$data[, 1] <- 1:6513
dtrain <- lgb.Dataset(train$data, label = train$label)
data(agaricus.test, package = "lightgbm")
test <- agaricus.test
dtest <- lgb.Dataset.create.valid(dtrain, test$data, label = test$label)
valids <- list(test = dtest)

params <- list(objective = "regression",
               metric = "rmse",
               device = "gpu",
               gpu_platform_id = 0,
               gpu_device_id = 0,
               nthread = 1,
               boost_from_average = FALSE,
               num_tree_per_iteration = 10,
               max_bin = 32)
model <- lgb.train(params,
                   dtrain,
                   2,
                   valids,
                   min_data = 1,
                   learning_rate = 1,
                   early_stopping_rounds = 10)

```

If it crashes, try to play around `gpu_platform_id` and `gpu_device_id`.

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

Please have a read if you intend to use Intel Compilers. Intel Compilers vs gcc: https://sites.google.com/view/lauraepp/benchmarks/intel-compiler-may-2018 (tl;dr: do not waste all those hours trying it, the benefit is very small unless you are running a very CPU-limited machine such as a laptop)

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
