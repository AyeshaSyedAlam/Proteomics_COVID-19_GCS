################################################################################
## Immune-metabolic plasma proteomic signatures distinguish asymptomatic from
## symptomatic SARS-CoV-2 infection and predict post-COVID sequelae
##
## Analysis code
## Gutenberg COVID-19 Study (GCS) / Gutenberg Post-COVID Study (GPC)
##
## Contents
##   0  Setup: packages, paths, labels
##   1  Functions
##        1a  Penalised regression: E_Net1() and plot_enet()
##        1b  Ensemble mRMR selection frequency
##        1c  Selection stability (Nogueira)
##        1d  Helper functions for the downstream analyses
##   2  Data preparation
##        2a  Absolute peptide quantification: transformation, winsorisation
##        2b  Relative protein quantification and protein annotation
##        2c  Cohort data: GCS baseline, GCS follow-up, GHS examination 6
##        2d  Variable definitions
##        2e  COVID-19 variables: severity, sequelae, vaccination, testing
##        2f  Treatment data (GPC) and computer-assisted interviews (CATI)
##   3  Baseline table
##   4  Ensemble mRMR feature selection and selection stability
##   5  Ridge regression and the protein score
##   6  Protein signature: ridge coefficients (Figure 2A)
##   7  Protein score by infection severity (Figure 2, Supplementary Figure 3)
##   8  Pre- and post-infection protein dynamics (Supplementary Figures 4, 5)
##   9  Protein score and acute symptoms (Figure 3A, Figure 4)
##  10  Symptom prevalence threshold: sensitivity analysis (Suppl. Figure 6)
##  11  Protein score and post-COVID sequelae (Figure 3B)
##  12  Protein score and long-term symptoms (Supplementary Figure 9)
##  13  Protein score and cardiovascular risk factors
##  14  Per-protein statistics across severity groups (Suppl. Figure 10)
##  15  Analyses added during revision
##        15a  Serologically confirmed infections
##        15b  Baseline comparison: with vs without proteomic measurements
##        15c  Analysis 1: sensitivity to the method of infection ascertainment
##        15d  Analysis 2: objective severity definition (treatment)
##        15e  Analysis 3: extended confounder adjustment
##        15f  Forest plots for Analysis 3
##        15g  Export of the analysis samples
##
## Sections 1a to 5, 15c, 15d and 15e were written by Robin Martens; the
## hospitalisation model in 15d by Alexander Gieswinkel.
##
## Notes
##   - File paths and database identifiers are given as placeholders. The
##     data are not public; access is described in the data availability
##     statement of the manuscript.
##   - E_Net1_RM(), baseline13(), reg_tab(), reg_tab_list(), reg_tab_loop(),
##     TAB_miss() and plot_save() are wrappers from an in-house function
##     library that is not reproduced here. What each one computes is
##     documented where it is called, so that every model can be reproduced
##     with standard R functions. E_Net1(), the penalised regression function
##     that E_Net1_RM() wraps, is in section 1a.
##       reg_tab()       one regression model, as a table of effect
##                       estimates, confidence intervals and p-values
##       reg_tab_list()  the same model under several adjustment sets
##       reg_tab_loop()  the same model over a list of predictors
##       TAB_miss()      counts of missing values per variable
##     regression = 'robust.poisson' fits a Poisson model with a log link and
##     heteroscedasticity-consistent (HC0) standard errors, so that the
##     exponentiated coefficient is a prevalence ratio; regression = 'linear'
##     and 'logistic' fit lm() and glm(family = binomial).
##   - The protein score is written out in section 5 and read back in from
##     there on, so that the feature selection and the ridge fit do not have
##     to be repeated.
################################################################################


################################################################################
## 0  SETUP
################################################################################

## ---- packages ----------------------------------------------------------

library(readr)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(smplot2)
library(ggbeeswarm)
library(forplo)
library(sandwich)
library(openxlsx)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
library(stringr)
library(e1071)      # skewness, for the choice of transformation
library(mRMRe)      # ensemble mRMR feature selection
library(glmnet)     # ridge regression, via E_Net1()
library(ggdist)
library(pROC)
library(Hmisc)
library(doParallel)
library(survival)

## ---- in-house function library (not reproduced here) --------------------

dir_fun <- "<path to the in-house function library>"
suppressWarnings(eval(parse(file = file.path(dir_fun, "ade.fun.run.r"))))

## ---- paths -------------------------------------------------------------

dir_ms   <- "<path to the mass spectrometry data>"
dir_abs  <- "<path to the absolute protein values>"
dir_rel  <- "<path to the relative protein values>"
dir_proj <- "<path to the project data>"
dir_sql  <- "<path to the database read-in scripts>"
files    <- "<path to the exported data files>"

out <- "<path to the output folder>"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

f_prot   <- file.path(files, "data_ensemeble_feature_selection_proteins_and_peptides2026-01-22.csv")
f_scores <- file.path(files, "ensemble_feature_selection_scores_2025-08-06.csv")
f_covid  <- file.path(files, "Covid_Variables.csv")
f_sig    <- file.path(files, "relative_protein_signature_new_20260123.csv")
f_seq    <- file.path(files, "6_month_seq.csv")

stopifnot(all(file.exists(c(f_prot, f_scores, f_covid, f_sig, f_seq))))

## Cohort data frame expected in the workspace:
##   d   baseline examination of the GCS, with the serology

## ---- symptom and severity labels ---------------------------------------

symp_vars <- paste0("symp_ynu", 1:61)

symp_names <- c('Fatigue','Malaise','Fever','Chills','Cough','Sore throat',
                'Hoarseness','Runny nose','Headache','Shortness of breath',
                'Pain when breathing','Wheezing','Chest pain','Tachycardia',
                'Olfactory disturbances','Taste disturbances','Visual disturbances',
                'Conjunctivitis','Hearing problems','Tinnitus','Memory problems',
                'Confusion','Concentration problems','Sleep disturbances',
                'Less sleep','More sleep','Dizziness','Fainting spells',
                'Joint pain or swelling','Swollen ankles','Muscle pain',
                'Muscle stiffness','Body aches','Twiching of the limbs',
                'Weakness of the limbs','Tremor','Numbness or tingling',
                'Mood swings','Depression','Loss of interest or pleasure',
                'Anxiety','Abdominal pain','Diarrhea','Constipation',
                'Nausea or vomiting','Loss of appetite','Weight loss',
                'Skin rash','Hair loss','Seizures','Problems urinating',
                'Problems swallowing','Problems with balance',
                'Trouble walking or falling','Slow movements',
                'Painful menstrual bleeding','Erectile dysfunction',
                'Behavioural changes','Hallucinations','COVID toes',
                'Partial face or body paralysis')
symp_names <- gsub(" ", "_", symp_names)

severity_labels <- c("Asymptomatic", "Mild", "Moderate", "Severe")


################################################################################
## 1  FUNCTIONS
################################################################################

## -----------------------------------------------------------------------
## 1a  Penalised regression
##     Elastic net or ridge with the adjustment variables left
##     unpenalised, a cross-validated penalty parameter, and a ranking of
##     the selected variables by the penalty at which each enters the
##     model. plot_enet() draws the cross-validation curve and saves it
##     with plot_save().
## -----------------------------------------------------------------------

#/¯\___/¯\___/¯\____/¯\____/¯\__E__/¯\_____/¯\_____/¯\____/¯\____/¯\___/¯\___/¯\_
#_/¯\___/¯\___/¯\____/¯\____/¯\__L__/¯\_____/¯\_____/¯\____/¯\____/¯\___/¯\___/¯\_
#__/¯\___/¯\___/¯\____/¯\____/¯\__A__/¯\__N__/¯\_____/¯\____/¯\____/¯\___/¯\___/¯\_
#___/¯\___/¯\___/¯\____/¯\____/¯\__S__/¯\__E__/¯\__1__/¯\____/¯\____/¯\___/¯\___/¯\_
#\___/¯\___/¯\___/¯\____/¯\____/¯\__T__/¯\__T__/¯\_____/¯\____/¯\____/¯\___/¯\___/¯\_
#¯\___/¯\___/¯\___/¯\____/¯\____/¯\__I__/¯\_____/¯\_____/¯\____/¯\____/¯\___/¯\___/¯\_
#/¯\___/¯\___/¯\___/¯\____/¯\____/¯\__C__/¯\_____/¯\_____/¯\____/¯\____/¯\___/¯\___/¯\_


E_Net1 <- function(v_dep,
n_dep = NULL,
v_adj = c(), 
n_adj = c(), 
v_pen, 
n_pen, 
dat, 
t_dep = NULL, 
family='binomial',
lambda_eval = 'lambda.min',
alpha=0.5, 
nfolds=10, 
loocv = F, 
round = 3,
parallel = NULL,
seed = 4815, 
maxit = 10^6){
  
# Description
# v_ = variable name in data frame
# n_ = variable label we like to see in the output
# dep = dependent variable
# adj = adjustment variables
# pen = penalized variables (pool of variables we want to select the most important ones)
# dat = dataframe
# round = where to round the digits in the final table
# seed = sets the seed individually, i.e. if you want to compare different 10-fold CV 
#        outcomes to see, if there are stable and LOOCV can be omitted. 

# Features:
# - Performes elastic net with alpha (0.5 = 'typical elastic net', 1 = LASSO, 0 = RIDGE) as mixing parameter.
# - All adjustment and penalized variables are scaled and centered within the function
#   before the glm function is applied (like in Andreas function LASSO_FP4)
# - Calculates a robustness for each LR (= % of models between the model of first appearance of the 
#   corresponding variable and the model with the optimal lambda, where the variable is not 0)
# - Crossvalidation can be parallelized by setting the 'parallel' parameter to the wanted number of 
#   cores that should be used. 


# Ranking of the table ['What is the lambda ratio?']: 
# 
# Lambda is the weight for the additional penalty term and can therefore takes values from
# 0 = no penalty (i.e. this equals the 'regular' regression model) to (theoretically)
# infinity = maximum penalty (i.e. all coefficients are 0 and we call this the 'null model').
# [Note: There are applications with negative values for Lambda, but we ignore this here - see further:
# https://stats.stackexchange.com/questions/331264/understanding-negative-ridge-regression
# if I have time, I will try to dig a bit into this :)]
#   The function starts with the lowest lambda value for the null model and then slowly decreases lambda to
# eventually 0. On this way, the model gets more and more non-zero coefficents, i.e. the model includes with
# decreasing lambda an increasing amount of (the penalized) variables. 
# The 'optimal' lambda is then defined as the value that minimizes the loss (e.g. mean-sqaure error in the linear regression).
#   Now, we can determine for each variable the highest lambda, that this variable is in the model and then 
# devide it by the optimal lambda (-> this is the lambda ratio and its lowest possible value is 1).
# This can be used for a ranking of the penalized variables: 
#   If variable A has a higher lambda ratio than another variable B we can interpret this that A was included to 
#   the model befor B. 


  require(glmnet)
  require(pROC)
  require(Hmisc)
  require(doParallel)
  require(epade)
  require(survival)
  
  
  if(is.null(t_dep)){
    dat <- dat[, c(v_dep, v_adj, v_pen)]
    for (i in 1:ncol(dat)) {
      if(class(dat[,i]) == 'factor') dat[,i] <- as.numeric(dat[,i])-1
    }
    dat <- na.omit(dat)
    
    
    ### Exclude Variables with constant value, otherwise scale produces an error
    out_adj <- c()
    for (j in 1:length(v_adj)) {
      if(length(unique(dat[,v_adj[j]])) == 1) out_adj <- c(out_adj, j)
    }
    out_pen <- c()
    for (j in 1:length(v_pen)) {
      if(length(unique(dat[,v_pen[j]])) == 1) out_pen <- c(out_pen, j)
    }
    
    if(length(out_adj)>0) {
      v_adj <- v_adj[-out_adj]
      n_adj <- n_adj[-out_adj]
    }
    if(length(out_pen)>0) {
      v_pen <- v_pen[-out_pen]
      n_pen <- n_pen[-out_pen]
    } 
    
    
    dat <- dat[, c(v_dep, v_adj, v_pen)]
    
    
    y <- dat[,v_dep]
    
  } else {
    dat <- dat[which(dat[,t_dep] > 0), c(v_dep, t_dep, v_adj, v_pen)]
    
    for (i in 1:ncol(dat)) {
      # if((i > 2) & (i <= (length(v_adj)+2)) & (class(dat[,i]) == 'numeric')) dat[,i] <- dat[,i]/sd(dat[,i])
      if(class(dat[,i]) == 'factor') dat[,i] <- as.numeric(dat[,i])-1
    }
    
    dat <- na.omit(dat)
    
    ### Exclude Variables with constant value, otherwise scale produces an error
    out_adj <- c()
    for (j in 1:length(v_adj)) {
      if(length(unique(dat[,v_adj[j]])) == 1) out_adj <- c(out_adj, j)
    }
    out_pen <- c()
    for (j in 1:length(v_pen)) {
      if(length(unique(dat[,v_pen[j]])) == 1) out_pen <- c(out_pen, j)
    }
    
    if(length(out_adj)>0) {
      v_adj <- v_adj[-out_adj]
      n_adj <- n_adj[-out_adj]
    }
    if(length(out_pen)>0) {
      v_pen <- v_pen[-out_pen]
      n_pen <- n_pen[-out_pen]
    } 
    
    
    dat <- dat[, c(v_dep, t_dep, v_adj, v_pen)]
    
    
    y <- cbind(dat[,t_dep], dat[, v_dep]) # need matrix which is made by cbind
    colnames(y)<- c('time', 'status')
  } 
  
  
  BigX <- scale(as.matrix(dat[,c(v_adj, v_pen)]))
  pen_ration <- c(rep(0, length(v_adj)), rep(1, length(v_pen)))
  
  if(!is.null(v_adj)){
    if(length(pen_ration)==ncol(BigX)){
      cat('Penalty was set to 0 for adjusting variables.\n')
    }
  }
  
  if(loocv == T) nfolds <- nrow(BigX)
  
  slabels <- c(n_adj, n_pen)
  
  set.seed(seed)
  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ### Set criterium what lambda should minimize 

  family <- tolower(family)
  if(family %in% c('linear', 'gauss')) family <- 'gaussian'
  if(family %in% c('logistic', 'binary')) family <- 'binomial'
  if(family %in% c('tte')) family <- 'cox'
  
  mykeep<- TRUE
  if(family=="binomial" | family=="cox") typem="deviance"
  if(family=="gaussian") typem="mse" #mae or mse
  if(family=="mgaussian") typem="mse" #mae or mse
  if(family=="poisson" | family == 'quasipoisson') typem="default" 
  if(family=="cox") mykeep<- F
  
  cat(paste0('Running ', family, ' E-Net with ', nfolds, '-fold cross-validation.\n'))
  
  ### Parallelize crossvalidation
  if(is.null(parallel)) {
    crossval       <- cv.glmnet(x = BigX, y=y, family = family, grouped=T, nlambda = 1000, type.measure=typem, nfolds=nfolds, alpha = alpha, standardize = F, penalty.factor=pen_ration, keep=mykeep, maxit = maxit,
                                parallel = F)
  } else {
    detect <- parallel::detectCores()
    # At least 2 cores should not be used for the parallelization
    if(parallel > detect-2) parallel <- detect - 2
    registerDoParallel(parallel)
    crossval       <- cv.glmnet(x = BigX, y=y, family = family, grouped=T, nlambda = 1000, type.measure=typem, nfolds=nfolds, alpha = alpha, standardize = F, penalty.factor=pen_ration, keep=mykeep, maxit = maxit,
                                parallel = T)
    stopImplicitCluster()
  }
  
  lambda         <- crossval[[lambda_eval]] #optimal lambda, minimal shrinkage
  fit            <- glmnet(BigX, y=y, alpha = alpha, nlambda = 1000, family=family, standardize = F, penalty.factor=pen_ration, maxit = maxit)
  #fit2           <- glmnet(BigX, y=y, alpha = alpha, nlambda = 5000,   family=family, lambda.min.ratio = 0.00000001, standardize = F, penalty.factor=pen_ration)
  pred_0         <- predict(fit,  newx=BigX, s=lambda,  type="response")
  #print(coef(fit, s=lambda))
  if(sum(abs(coef(fit, s=lambda))[-1]>0)==0){
    return(F)
    stop('Nothing was selected.')
  }
  if(family=="cox") {
    library(foreach)
    events <- sum(y[,2])
    # normal (overoptimistic) C-index
    pred0 <- predict(fit, newx=BigX, s = lambda, type="response")
    AUC_0 <- as.numeric(as.vector((rcorr.cens(-pred0[,1],Surv(y[,1],y[,2]))[1])))
    # AUC_1 <- as.numeric(as.vector((rcorr.cens(-as.vector(crossval$fit.preval[,atds]),Surv(y[,1],y[,2]))[1])))
    
    # Cross-validation based C-index
    predCV <- rep(NA,nrow(BigX))
    cat(paste0('Cross-validation run for C-index: '))
    
    ### 
    if(loocv == F) {
      folds <- sample(rep(1:nfolds,ceiling(nrow(BigX)/nfolds))[1:nrow(BigX)])
      for (i in 1:10){
        cat(paste0(i, ', '))
        temp = glmnet(BigX[which(folds!=i),], y=y[which(folds!=i),], alpha = alpha, lambda = crossval$lambda, family="cox", standardize = F, penalty.factor=pen_ration)
        predCV[which(folds==i)] = predict(temp, newx=BigX[folds==i,],s=lambda, type="response")
      }
    } else {
      if(!is.null(parallel)) {
        # start <- Sys.time()
        detect <- parallel::detectCores()
        # At least 2 cores should not be used for the parallelization
        if(parallel > detect-2) parallel <- detect - 2
        registerDoParallel(parallel)
        preds <- foreach(i = 1:nfolds, .packages = 'glmnet', .combine = c) %dopar% {
          temp <- glmnet(BigX[-i,], y=y[-i,], alpha = alpha, lambda = crossval$lambda, family="cox", standardize = F, penalty.factor=pen_ration)
          predict(temp, newx=BigX[i,], s=lambda, type="response")
        }
        stopImplicitCluster()
        # end <- Sys.time()
        # end-start
        predCV <- preds
      } else {
        folds <- 1:nfolds
        for (i in 1:nfolds){
          if(i %% 25 == 0) cat(paste0(i, ', '))
          temp = glmnet(BigX[which(folds!=i),], y=y[which(folds!=i),], alpha = alpha, lambda = crossval$lambda, family="cox", standardize = F, penalty.factor=pen_ration)
          predCV[which(folds==i)] = predict(temp, newx=BigX[folds==i,],s=lambda, type="response")
        }
      }
    }
    
    cat('done! :)\n')
    AUC_cv <- as.numeric(as.vector((rcorr.cens(-predCV,Surv(y[,1],y[,2]))[1])))
    coef_name <- 'HR per SD'
    cat(paste0('Simple C-index: ', round(AUC_0, 4), '; ', nfolds, '-fold CV C-index: ', round(AUC_cv, 4), '\n'))
  }
  if(family=="binomial"){
    events <- sum(y)
    AUC_0  <- as.numeric(as.vector(roc(y,  as.vector(pred_0), direction='<')$auc))
    atds   <- which(crossval$lambda==crossval[[lambda_eval]])
    AUC_cv <- as.numeric(as.vector(roc(y,  as.vector(crossval$fit.preval[,atds]), direction='<')$auc))
    coef_name <- 'OR per SD'
    cat(paste0('Simple AUC: ', round(AUC_0, 4), '; CV AUC: ', round(AUC_cv, 4), '\n'))
  }
  if(family=="gaussian"){
	  events <- NA
    AUC_0  <- cor(y, pred_0)^2
    atds   <- which(crossval$lambda==crossval[[lambda_eval]])
    AUC_cv <- cor(y, as.vector(crossval$fit.preval[,atds]))^2
    coef_name <- 'Beta per SD'
    cat(paste0('Simple R^2: ', round(AUC_0, 4), '; CV R^2: ', round(AUC_cv, 4), '\n'))
  }
  if(family == 'poisson' | family == 'quasipoisson') {
    events <- NA
    coef_name <- 'IRR per SD'
    AUC_0  <- cor(y, pred_0)^2
    atds   <- which(crossval$lambda==crossval[[lambda_eval]])
    AUC_cv <- cor(y, as.vector(crossval$fit.preval[,atds]))^2
    cat(paste0('Simple pseudo R^2: ', round(AUC_0, 4), '; CV pseudo R^2: ', round(AUC_cv, 4), '\n'))
  }
  if(family!="cox"){
    # vlabels        <- c('Intercept^1', vlabels)
    slabels        <- c('Intercept',   slabels)
    # plabels        <- c('^1',          plabels)
    signs          <- sign(coef(fit, s=lambda))
    signs_txt      <- rep('', length(signs))
  }
  if(family=="cox"){
    # vlabels        <- vlabels[which(vlabels!='^1')]
    slabels        <- c(slabels)
    # plabels        <- c(plabels)
    signs  <- sign(coef(fit, s=lambda))
    signs_txt <- rep(length(signs))
  }
  signs_txt[which(signs==-1)] <- '-'
  signs_txt[which(signs==0)]  <- ''
  signs_txt[which(signs==1)]  <- '+'
  #########################################
  
  #########################################
  # Ranking
  LAM <- matrix(fit$lambda, nrow(fit$beta), ncol(fit$beta), byrow=T) # Row = Variables, Column = Lambdas
  LAM[which(fit$beta==0, arr.ind=F)] <- 0 # set Lambdas with beta coefficients 0 to 0 
  ranks <- apply(LAM, 1, max, na.rm=T) # maximum of every row of LAM
  ranks <- ranks/lambda # 'Lambda ratio' - first appear divided by best model
  
  # Trust in lambda ratio
  s.lambda <- which(fit$lambda == lambda)
  LAM2 <-LAM
  LAM2 <- ifelse(abs(LAM2) > 0, 1, 0)
  trust <- rep(0, nrow(LAM2))
  for(i in 1:length(trust)) {
    k <- 1
    while(LAM2[i,k] == 0 & k < s.lambda) {
      k <- k+1
    }
    st <- sum(LAM2[i,k:s.lambda])
    trust[i] <- st/(s.lambda-k+1)
  }

  
  #coefs2 <- coef(fit, s=fit$lambda[s.lambda], exact = T)
  coefs <- coef(fit, s=lambda)
  
  # There is no intercept in the Cox model as it is built into the baseline hazard
  # See: https://glmnet.stanford.edu/articles/Coxnet.html
  if(family!="cox") {
    ranks <- c(NA, ranks)
    trust <- c(0, trust)
  }
  
  #print(paste0('min lambda=', round(lambda,5)))
  sids <- which(abs(coefs)>0)
  M <- cbind(slabels[sids], signs_txt[sids], ranks[sids], trust[sids], coefs[sids])
  if(dim(M)[1]>1){
    M <- M[base:::order(as.numeric(M[,3]), decreasing = T),]
    M <- M[which(M[,1]!='Intercept'), , drop=FALSE]
    M2 <- M
    colnames(M) <- c('Variable', 'Direction','Lambda Ratio', 'Robustness of LR [%]',coef_name)
    
    
    M[,3] <- round_n.ade(as.numeric(as.character(M[,3])), round)
    M[,4] <- round_n.ade(as.numeric(as.character(M[,4]))*100, 2)
    if(family %in%c('cox', 'binomial', 'poisson')) M[,5]<- round_n.ade(exp(as.numeric(as.character(M[,5]))), round)
    if(family %in%c('gaussian', 'mgaussian')) M[,5]<- round_n.ade(as.numeric(as.character(M[,5])), round)
    M[which(M[,1]%in%n_adj),3]<- '--'
    
    selected <- v_pen[which(n_pen%in%M[,1])]
    
  }
  else
  {
    colnames(M) <- c('Variable', 'Direction','Lambda Ratio', 'Coefficient', coef_name)
    M2 <- M
  }
  #########################################
  
  selected <- M[which(M[,3]!='--'),1]
  
  cat(paste0('Selected ', length(selected), ' / ', sum(pen_ration), ' penalized variables.\n'))
  return(list(M, fit, ranks, AUC_0, AUC_cv, selected, slabels, signs_txt, M2, lambda, BigX, pen_ration, crossval, nrow(BigX), events, y, family, alpha, nfolds, loocv, n_dep, n_adj))
}


# results <- E_Net1('bmi', NULL, c('sex', 'age'), c('Sex', 'Age'), vars_prot, nams_prot, d, family = 'gaussian')


plot_enet <- function(results, plot_name, se = F, xlim = NULL, ylim = NULL, zoom = F){
  if(zoom) {
    xlim <- NULL
    ylim <- NULL
  }
  
  M              <- results[[1]]
  fit            <- results[[2]]
  # AUC_0          <- results[[4]]
  # AUC_cv         <- results[[5]]
  lambda         <- results[[10]]
  # BigX           <- results[[11]]
  # pen_ration     <- results[[12]]
  crossval <- cvobj <- results[[13]]
  # N              <- results[[14]]
  # events         <- results[[15]]
  family         <- results[[17]]
  
  lambda_se <- crossval$lambda.1se
  
  
  n_adj <- length(M[which(M[,3]=='--'),1])
  
  if(family != 'cox') sum_sel <- sum(coef(fit, s=lambda)[-1] != 0)
  if(family == 'cox') sum_sel <- sum(coef(fit, s=lambda) != 0)
  if(se) {
    if(family != 'cox') sum_sel_se <- sum(coef(fit, s=lambda_se)[-1] != 0)
    if(family == 'cox') sum_sel_se <- sum(coef(fit, s=lambda_se) != 0)
  }
  
  if(zoom) {
    pos1 <- which(crossval$lambda == crossval$lambda.min)
    pos2 <- which(crossval$lambda == crossval$lambda.1se)
    diff <- pos1-pos2
    limit <- (pos1-diff):(pos1+diff)
  } else {
    limit <- 1:length(crossval$lambda)
  }
  
  
  dev.new(noRStudioGD = T)
  sign.lambda <- 1
  xlab = expression(Log(lambda))
  xlab = paste("-", xlab, sep = "")
  if(is.null(ylim)) ylim <- range(cvobj$cvup[limit], cvobj$cvlo[limit])
  if(!is.null(ylim)) ylim <- ylim

  plot(x = sign.lambda * log(cvobj$lambda[limit]), y = cvobj$cvm[limit], 
       ylim = ylim, xlab = xlab, ylab = cvobj$name, 
       type = "n", xlim = xlim, las = 1)
  glmnet:::error.bars(sign.lambda * log(cvobj$lambda[limit]), cvobj$cvup[limit], 
             cvobj$cvlo[limit], width = 0.01, col = "darkgrey")
  points(sign.lambda * log(cvobj$lambda), cvobj$cvm, pch = 20, 
         col = "red")
  axis(side = 3, at = sign.lambda * log(cvobj$lambda), labels = paste(cvobj$nz), 
       tick = FALSE, line = 0)
  abline(v = sign.lambda * log(cvobj$lambda.min), lty = 3)
  abline(v = sign.lambda * log(cvobj$lambda.1se), lty = 3)
  
 
  
  abline(v=log(lambda), col=2)
  mtext(side=3, at=log(lambda), text=paste0(sum_sel, ' (', sum_sel-n_adj, ')'), line=0.1, col=2)
  mtext(side=1, at=log(lambda), text=round(lambda,3), col=2, line=0)
  if(se) {
    abline(v=log(lambda_se), col=4)
    mtext(side=3, at=log(lambda_se), text=paste0(sum_sel_se, ' (', sum_sel_se-n_adj, ')'), line=0.1, col=4)
    mtext(side=1, at=log(lambda_se), text=round(lambda_se,3), col=4, line=0)
  }
  if(!zoom) mtext(side=3, text='Number of variables', line=2.5)
  if(zoom) mtext(side=3, text='Number of variables [zoomed in]', line=2.5)
  box()
  plot_save(fname = plot_name, type=c('png'))
}


## -----------------------------------------------------------------------
## 1b  Ensemble mRMR selection frequency
## -----------------------------------------------------------------------

selection_frequency <- function(dd, iterations, feature_count = 50, seed = 2025) {
  set.seed(seed)
  ensemble <- mRMR.ensemble(data = dd, target_indices = c(1),
                            solution_count = iterations,
                            feature_count = feature_count)
  ensemble_list <- solutions(ensemble, with_fixed_features = T)[[1]]

  rank <- rep(0, length(c("asym_smod_and_hard", vars_prot_rel)))
  for (i in 1:iterations) {
    for (j in 1:length(c(vars_prot_rel))) {
      if (j %in% ensemble_list[, i]) rank[j] <- rank[j] + 1 / iterations
    }
  }
  list(rank = rank, ensemble = ensemble, ensemble_list = ensemble_list)
}


## -----------------------------------------------------------------------
## 1c  Selection stability
##     Nogueira S, Sechidis K, Brown G. On the stability of feature
##     selection algorithms. J Mach Learn Res 2018;18:1-54.
## -----------------------------------------------------------------------

getStability <- function(X, alpha = 0.05) {
  X <- na.omit(X)
  M <- nrow(X)
  d <- ncol(X)
  hatPF  <- colMeans(X)
  kbar   <- sum(hatPF)
  v_rand <- (kbar / d) * (1 - kbar / d)
  stability <- 1 - (M / (M - 1)) * mean(hatPF * (1 - hatPF)) / v_rand

  ki    <- rowSums(X)
  phi_i <- rep(0, M)
  for (i in 1:M) {
    phi_i[i] <- (1 / v_rand) * ((1 / d) * sum(X[i, ] * hatPF) -
                                (ki[i] * kbar) / d^2 -
                                (stability / 2) * ((2 * kbar * ki[i]) / d^2 -
                                                   ki[i] / d - kbar / d + 1))
  }
  phi_bar  <- mean(phi_i)
  var_stab <- (4 / M^2) * sum((phi_i - phi_bar)^2)

  z     <- qnorm(1 - alpha / 2)
  upper <- stability + z * sqrt(var_stab)
  lower <- stability - z * sqrt(var_stab)

  ## correct finite sample problems
  stability <- max(0, stability)
  lower     <- max(0, lower)
  upper     <- max(0, upper)

  list("stability" = stability, "variance" = var_stab,
       "lower" = lower, "upper" = upper)
}


## -----------------------------------------------------------------------
## 1d  Helper functions for the downstream analyses
## -----------------------------------------------------------------------

## post-infection protein scores with age, sex and severity
get_scores <- function(timepoint = "post", vars = c("age", "sex", "cov_severity_a")) {
  s <- read.csv2(f_scores)
  s <- subset(s, pre_post_infection == timepoint)
  s$cov_prot_score_ensemble <- as.numeric(s$cov_prot_score_ensemble)
  merge(s, d_cati[c("sid", vars)], by = "sid", all.x = TRUE)
}

## signature protein levels (post-infection), columns renamed to gene symbols
get_proteins <- function(timepoint = "post") {
  p <- read.csv2(f_prot)[, 1:53]
  if (!is.null(timepoint)) p <- subset(p, pre_post_infection == timepoint)
  colnames(p)[5:53] <- read.csv2(f_sig)$gene_symbol
  p
}

## recode symptom variables: 99 and any non-binary value to NA
recode_symptoms <- function(df, vars) {
  for (v in vars) {
    df[[v]][df[[v]] == 99] <- NA
    df[[v]][!df[[v]] %in% c(0, 1) & !is.na(df[[v]])] <- NA
    df[[v]] <- as.numeric(df[[v]])
  }
  df
}

## linear models of the protein score on a set of binary predictors
score_on <- function(df, predictors, covariates = c("age", "sex")) {
  res <- lapply(predictors, function(v) {
    m <- lm(reformulate(c(v, covariates), response = "cov_prot_score_ensemble"),
            data = df)
    co <- summary(m)$coefficients[2, ]
    data.frame(predictor = v, beta = co[1], se = co[2], p = co[4],
               n = nobs(m), row.names = NULL)
  })
  bind_rows(res) %>%
    mutate(LCI = beta - 1.96 * se,
           UCI = beta + 1.96 * se)
}

## prevalence ratio from Poisson regression with robust (HC0) standard errors
prevalence_ratio <- function(outcome, data, covariates = c("age", "sex"), label = outcome) {
  form <- reformulate(c("cov_prot_score_ensemble", covariates), response = outcome)
  m    <- glm(form, data = data, family = poisson(link = "log"))
  se   <- sqrt(diag(vcovHC(m, type = "HC0")))
  est  <- coef(m)[2]; s <- se[2]
  data.frame(outcome = label,
             PR  = exp(est),
             LCI = exp(est - 1.96 * s),
             UCI = exp(est + 1.96 * s),
             p   = 2 * (1 - pnorm(abs(est / s))),
             n   = nobs(m), row.names = NULL)
}


################################################################################
## 2  DATA PREPARATION
################################################################################

## -----------------------------------------------------------------------
## 2a  Absolute peptide quantification: transformation and winsorisation
## -----------------------------------------------------------------------

data <- read.csv2(file.path(dir_abs, "absolut_peptide_values.csv"), sep = ";")[, -1]

## For every peptide: choose the transformation that minimises skewness
## (none, log, square root or inverse), then winsorise at six median absolute
## deviations from the median, spacing the replaced values by the mean
## distance between neighbouring observations.
d1 <- data

for (i in 2:ncol(d1)) {

  ## ---- transformation ----------------------------------------------------
  a   <- NULL
  eps <- 0.01
  x        <- as.numeric(d1[, i])
  x_intern <- c(na.omit(as.numeric(d1[, i])))
  if (is.null(a) == T) a <- max(0, (-min(x_intern) + eps * sd(x_intern, na.rm = T)))

  x_untrans      <- x
  x_untrans_info <- skewness(x_untrans, na.rm = T)
  x_log          <- log(x + a)
  x_log_info     <- skewness(x_log, na.rm = T)
  x_sqrt         <- sqrt(x + a)
  x_sqrt_info    <- skewness(x_sqrt, na.rm = T)
  x_inv          <- -(1 / (x + a))
  x_inv_info     <- skewness(x_inv, na.rm = T)

  t <- which.min(abs(c(x_untrans_info, x_log_info, x_sqrt_info, x_inv_info)))

  if (t == 1) d1[, i] <- x
  if (t == 2) d1[, i] <- x_log
  if (t == 3) d1[, i] <- x_sqrt
  if (t == 4) d1[, i] <- x_inv

  ## ---- winsorisation -----------------------------------------------------
  x   <- d1[, i]
  sds <- mad(x, na.rm = T)
  mea <- median(x, na.rm = T)
  u_ids <- which(x > (mea + 6 * sds))
  d_ids <- which(x < (mea - 6 * sds))

  mean_dis <- mean(diff(sort(na.omit(x)), 1))
  x[u_ids] <- mea + 6 * sds + mean_dis * order(x[u_ids])
  x[d_ids] <- mea - 6 * sds - mean_dis * order(x[d_ids])

  d1[, i] <- x
}

d1$Run <- gsub("Covid", "Cov", d1$Run)

## sample metadata and tube codes
meta      <- read.delim(file.path(dir_ms, "2024-10-15_AF_50-0144_metadata.txt"))
tubecodes <- read_excel(file.path(dir_ms, "GCS_singleCell_tubecodes_lids.xlsx"))

meta <- merge(meta, tubecodes, by.x = "Sample.Name", by.y = "tubeCode", all.x = T)
meta$File.Name <- gsub(".d", "", meta$File.Name)

dat_all_norm_winz <- merge(meta, d1, by.x = "File.Name", by.y = "Run", all.x = T)
dat_all_norm_winz <- dat_all_norm_winz[which(!is.na(dat_all_norm_winz$lid)), ]


## -----------------------------------------------------------------------
## 2b  Relative protein quantification and protein annotation
## -----------------------------------------------------------------------

## normalised, winsorised and imputed relative protein values
dat_prot_rel <- read.csv2(file.path(dir_rel, "MassSpec_data_norm_winz_imp.csv"))[, -1]
dat_prot_rel$File.Name <- gsub(".d", "", dat_prot_rel$File.Name)

## protein group annotation from the DIA-NN report
Names <- read_tsv(file.path(dir_ms, "report.pg_matrix.tsv"))
Names$Protein.Ids <- gsub(";", ".", Names$Protein.Ids, fixed = TRUE)
Names$Protein.Ids <- gsub("-", ".", Names$Protein.Ids, fixed = TRUE)


## -----------------------------------------------------------------------
## 2c  Cohort data
## -----------------------------------------------------------------------

## ---- GCS baseline ------------------------------------------------------
suppressWarnings(eval(parse(file = file.path(dir_sql, "GCS_read_SQL_1_S.R"))))
d_BL <- read_SQL_BL_GCS(reload = T, user = "<database user>")

dat_all_norm_winz_GCC_BL <- merge(dat_all_norm_winz, d_BL,
                                  by.x = "lid", by.y = "lid.x", all = FALSE)
dat_all_norm_winz_GCC_BL$hem <- dat_all_norm_winz_GCC_BL$hemt

## ---- GCS follow-up -----------------------------------------------------
suppressWarnings(eval(parse(file = file.path(dir_sql, "GCS_read_SQL_FU_1_S.R"))))
d_FU <- read_SQL_FU_GCS(reload = T, user = "<database user>")

## carry the baseline test result forward
d_BL_copy <- d_BL[, c("sid", "cov_test_res_c")]
colnames(d_BL_copy) <- c("sid", "cov_test_res_c_bl")
d_FU <- merge(d_FU, d_BL_copy[c("sid", "cov_test_res_c_bl")], by = "sid", all.x = T)

dat_all_norm_winz_GCC_FU <- merge(dat_all_norm_winz, d_FU,
                                  by.x = "lid", by.y = "lid.x", all = FALSE)
dat_all_norm_winz_GCC_FU$hem <- dat_all_norm_winz_GCC_FU$hemt

## ---- Gutenberg Health Study, examination 6 (uninfected comparison) -----
suppressWarnings(eval(parse(file = file.path(dir_sql, "Read_SQL_GHS_C6.r"))))
d_C6 <- read_SQL_C6_data(reload = T, read_proteins = F, u_id = "<database user>")

dat_all_norm_winz_GHS <- merge(dat_all_norm_winz, d_C6,
                               by.x = "lid", by.y = "lid", all = FALSE)

dat_all_norm_winz_GHS$chf <- dat_all_norm_winz_GHS$hrzcapi_c004a_ynu
dat_all_norm_winz_GHS$chf[which(dat_all_norm_winz_GHS$chf == 99)] <- NA
dat_all_norm_winz_GHS$chf[which(dat_all_norm_winz_GHS$chf == 0)]  <- "no"
dat_all_norm_winz_GHS$chf <- as.factor(dat_all_norm_winz_GHS$chf)
dat_all_norm_winz_GHS$egfr <- dat_all_norm_winz_GHS$gfr_ckd_epi
dat_all_norm_winz_GHS$auto <- NA

## ---- Olink panel (subset with both platforms) --------------------------
olink <- read.csv2(file.path(dir_proj, "Olink_data_norm_winz_imp.csv"))


## -----------------------------------------------------------------------
## 2d  Variable definitions
## -----------------------------------------------------------------------

vars_cvrf <- c("sex", "age", "bmi", "adipos",
               "nic", "hyper", "dyslip",
               "fhmis", "mi",
               "afib", "hemo",
               "egfr", "stroke",
               "cad", "pad",
               "hem", "cld",
               "copd", "cancer",
               "cvd", "auto",
               "nic", "dyslip",
               "diab")

nams_cvrf <- c("Sex", "Age [y]", "BMI", "Obesity",
               "Smoking", "Arterial hypertension", "Dyslipidemia",
               "FH of MI/Stroke", "MI",
               "Atrial fibrillation", "Hemoglobin",
               "eGFR [ml/min/1.73m2]", "Stroke",
               "Coronary artery disease", "Peripheral artery disease",
               "Hematocrit", "Chronic liver disease",
               "Chronic obstructive pulmonary disease", "History of cancer",
               "Cardiovascular disease", "Autoimmune disease",
               "Smoking", "Dyslipidemia",
               "Diabetes")

vars_peptid   <- colnames(data)[-1]
vars_prot_rel <- colnames(dat_prot_rel)[-which(colnames(dat_prot_rel) == "File.Name")]

nams_prot_rel       <- gsub("_HUMAN", "", Names$Protein.Names[which(Names$Protein.Ids %in% vars_prot_rel)])
nams_genes_prot_rel <- Names$Genes[which(Names$Protein.Ids %in% vars_prot_rel)]

vars <- c("lid", "File.Name", vars_cvrf, vars_peptid)

## one data frame across the three examinations
dat_norm_winz <- rbind(dat_all_norm_winz_GCC_BL[, vars],
                       dat_all_norm_winz_GCC_FU[, vars],
                       dat_all_norm_winz_GHS[, vars])


## -----------------------------------------------------------------------
## 2e  COVID-19 variables
## -----------------------------------------------------------------------

dat_covid_vars <- read.csv2(file.path(dir_proj, "Covid_Variables.csv"), sep = ",")[, -1]

## ---- GCS baseline ------------------------------------------------------
dat_cov_1 <- merge(dat_all_norm_winz_GCC_BL, dat_covid_vars,
                   by.x = "sid", by.y = "ID", all.x = T)

## if the infection was not yet reported at baseline, severity and sequelae
## belong to the follow-up examination
dat_cov_1[which(dat_cov_1$cov19 == "no"),
          c("Severity.of.acute.infection",
            "Sequelae.nach.3.Monaten", "Sequelae.nach.6.Monaten")] <- NA

## the same if the infection date lies after the baseline examination
dat_cov_1[which(dat_cov_1$"Infektionsdatum..FU." > dat_cov_1$date),
          c("Severity.of.acute.infection",
            "Sequelae.nach.3.Monaten", "Sequelae.nach.6.Monaten")] <- NA

dat_cov_1$"cov_test_res_c_bl" <- NA

## ---- GCS follow-up -----------------------------------------------------
dat_cov_2 <- merge(dat_all_norm_winz_GCC_FU, dat_covid_vars,
                   by.x = "sid", by.y = "ID", all.x = T)

vars_cov_bl <- c(vars, "date", "cov19", "Geimpft",
                 "Infektionsdatum", "Severity.of.acute.infection",
                 "Sequelae.nach.3.Monaten", "Sequelae.nach.6.Monaten",
                 "cov_test_type_c", "cov_test_res_c", "cov_test_res_c_bl")
vars_cov_fu <- c(vars, "date", "cov19_fu", "Geimpft..FU.",
                 "Infektionsdatum..FU.", "Severity.of.acute.infection",
                 "Sequelae.nach.3.Monaten", "Sequelae.nach.6.Monaten",
                 "cov_test_type_c", "cov_test_res_c", "cov_test_res_c_bl")

## ---- GHS examination 6: no infection reported --------------------------
dat_cov_3 <- dat_all_norm_winz_GHS
dat_cov_3$"cov19"                       <- "no"
dat_cov_3$Geimpft                       <- "no"
dat_cov_3$Infektionsdatum               <- NA
dat_cov_3$"Severity.of.acute.infection" <- NA
dat_cov_3$"Sequelae.nach.3.Monaten"     <- NA
dat_cov_3$"Sequelae.nach.6.Monaten"     <- NA
dat_cov_3$"cov_test_type_c"             <- NA
dat_cov_3$"cov_test_res_c"              <- NA
dat_cov_3$"cov_test_res_c_bl"           <- NA

dat_cov_1 <- dat_cov_1[, c(vars_cov_bl)]
dat_cov_2 <- dat_cov_2[, c(vars_cov_fu)]
colnames(dat_cov_2) <- colnames(dat_cov_1)
dat_cov_3 <- dat_cov_3[, c(vars_cov_bl)]
colnames(dat_cov_3) <- colnames(dat_cov_1)

dat_cov_1$FU <- "GCS BL"
dat_cov_2$FU <- "GCS FU"
dat_cov_3$FU <- "GHS C6"

dat_cov <- rbind(dat_cov_1, dat_cov_2, dat_cov_3)

dat_cov$Severity.of.acute.infection <-
  factor(dat_cov$Severity.of.acute.infection,
         levels = c("asymptomatisch", "milder Verlauf",
                    "maessiger Verlauf", "schwerer Verlauf"))

## time between infection and blood draw
dat_cov$time_after_inf <- log(as.numeric(as.Date(dat_cov$date) -
                                         as.Date(dat_cov$Infektionsdatum)) + 1)

## join the relative protein values and the Olink panel
dat_cov_all <- merge(dat_cov, dat_prot_rel, by = "File.Name", all = FALSE)
dat_cov_all <- merge(dat_cov_all, olink, by.x = "lid", by.y = "LID", all.x = F)

## ---- severity contrasts ------------------------------------------------

dat_cov_all$asym <- NA
dat_cov_all$asym[which(dat_cov_all$Severity.of.acute.infection == "asymptomatisch")] <- "yes"
dat_cov_all$asym[which(dat_cov_all$Severity.of.acute.infection != "asymptomatisch")] <- "no"
dat_cov_all$asym <- as.factor(dat_cov_all$asym)
## asymptomatic = 0, symptomatic = 1
dat_cov_all$asym2 <- 2 - as.numeric(dat_cov_all$asym)

dat_cov_all$mild <- NA
dat_cov_all$mild[which(dat_cov_all$Severity.of.acute.infection == "milder Verlauf")] <- "yes"
dat_cov_all$mild[which(dat_cov_all$Severity.of.acute.infection != "milder Verlauf")] <- "no"
dat_cov_all$mild <- as.factor(dat_cov_all$mild)

dat_cov_all$maes <- NA
dat_cov_all$maes[which(dat_cov_all$Severity.of.acute.infection == "maessiger Verlauf")] <- "yes"
dat_cov_all$maes[which(dat_cov_all$Severity.of.acute.infection != "maessiger Verlauf")] <- "no"
dat_cov_all$maes <- as.factor(dat_cov_all$maes)

dat_cov_all$schw <- NA
dat_cov_all$schw[which(dat_cov_all$Severity.of.acute.infection == "schwerer Verlauf")] <- "yes"
dat_cov_all$schw[which(dat_cov_all$Severity.of.acute.infection != "schwerer Verlauf")] <- "no"
dat_cov_all$schw <- as.factor(dat_cov_all$schw)

dat_cov_all$sex_n <- as.numeric(dat_cov_all$sex) - 1

## English labels for the severity categories
dat_cov_all$Severity <- dat_cov_all$Severity.of.acute.infection
levels(dat_cov_all$Severity)[levels(dat_cov_all$Severity) == "asymptomatisch"]    <- "asymptomatic"
levels(dat_cov_all$Severity)[levels(dat_cov_all$Severity) == "milder Verlauf"]    <- "mild"
levels(dat_cov_all$Severity)[levels(dat_cov_all$Severity) == "maessiger Verlauf"] <- "moderate"
levels(dat_cov_all$Severity)[levels(dat_cov_all$Severity) == "schwerer Verlauf"]  <- "severe"

## asymptomatic versus every symptomatic category
dat_cov_all$asym_smod_and_hard <- NA
dat_cov_all$asym_smod_and_hard[which(dat_cov_all$Severity.of.acute.infection == "asymptomatisch")]    <- "yes"
dat_cov_all$asym_smod_and_hard[which(dat_cov_all$Severity.of.acute.infection == "milder Verlauf")]    <- "no"
dat_cov_all$asym_smod_and_hard[which(dat_cov_all$Severity.of.acute.infection == "maessiger Verlauf")] <- "no"
dat_cov_all$asym_smod_and_hard[which(dat_cov_all$Severity.of.acute.infection == "schwerer Verlauf")]  <- "no"
dat_cov_all$asym_smod_and_hard <- 2 - as.numeric(as.factor(dat_cov_all$asym_smod_and_hard))

## mild versus moderate and severe
dat_cov_all$mild_vs_mod_and_hard <- NA
dat_cov_all$mild_vs_mod_and_hard[which(dat_cov_all$Severity.of.acute.infection == "mild")]             <- "yes"
dat_cov_all$mild_vs_mod_and_hard[which(dat_cov_all$Severity.of.acute.infection == "maessiger Verlauf")] <- "no"
dat_cov_all$mild_vs_mod_and_hard[which(dat_cov_all$Severity.of.acute.infection == "schwerer Verlauf")]  <- "no"
dat_cov_all$mild_vs_mod_and_hard <- 2 - as.numeric(as.factor(dat_cov_all$mild_vs_mod_and_hard))

## pairwise contrasts between all severity categories
dat_cov_all$asymvsmild <- NA
dat_cov_all$asymvsmild[which(dat_cov_all$Severity.of.acute.infection == "asymptomatisch")] <- "asymptomatic"
dat_cov_all$asymvsmild[which(dat_cov_all$Severity.of.acute.infection == "milder Verlauf")] <- "mild"
dat_cov_all$asymvsmild <- as.factor(dat_cov_all$asymvsmild)

dat_cov_all$asymvsmaes <- NA
dat_cov_all$asymvsmaes[which(dat_cov_all$Severity.of.acute.infection == "asymptomatisch")]    <- "asymptomatic"
dat_cov_all$asymvsmaes[which(dat_cov_all$Severity.of.acute.infection == "maessiger Verlauf")] <- "moderate"
dat_cov_all$asymvsmaes <- as.factor(dat_cov_all$asymvsmaes)

dat_cov_all$asymvshard <- NA
dat_cov_all$asymvshard[which(dat_cov_all$Severity.of.acute.infection == "asymptomatisch")]   <- "asymptomatic"
dat_cov_all$asymvshard[which(dat_cov_all$Severity.of.acute.infection == "schwerer Verlauf")] <- "hard"
dat_cov_all$asymvshard <- as.factor(dat_cov_all$asymvshard)

dat_cov_all$mildvsmaes <- NA
dat_cov_all$mildvsmaes[which(dat_cov_all$Severity.of.acute.infection == "milder Verlauf")]    <- "mild"
dat_cov_all$mildvsmaes[which(dat_cov_all$Severity.of.acute.infection == "maessiger Verlauf")] <- "moderate"
dat_cov_all$mildvsmaes <- as.factor(dat_cov_all$mildvsmaes)

dat_cov_all$mildvshard <- NA
dat_cov_all$mildvshard[which(dat_cov_all$Severity.of.acute.infection == "milder Verlauf")]   <- "mild"
dat_cov_all$mildvshard[which(dat_cov_all$Severity.of.acute.infection == "schwerer Verlauf")] <- "hard"
dat_cov_all$mildvshard <- as.factor(dat_cov_all$mildvshard)

dat_cov_all$maesvshard <- NA
dat_cov_all$maesvshard[which(dat_cov_all$Severity.of.acute.infection == "maessiger Verlauf")] <- "moderate"
dat_cov_all$maesvshard[which(dat_cov_all$Severity.of.acute.infection == "schwerer Verlauf")]  <- "hard"
dat_cov_all$maesvshard <- as.factor(dat_cov_all$maesvshard)

## ---- method of infection ascertainment ---------------------------------
## Participants without a reported infection are counted as negative; only
## five participants had a negative test result on record.
dat_cov_all$pcr_qPCR <- NA
dat_cov_all$pcr_qPCR[which(dat_cov_all$cov_test_type_c == 1)] <- "yes"
dat_cov_all$pcr_qPCR[which(dat_cov_all$cov_test_res_bl == 2 &
                           dat_cov_all$cov_test_res == 2)]    <- "yes"
dat_cov_all$pcr_qPCR[which(dat_cov_all$cov19 == "no")]        <- "no"

table(dat_cov_all$pcr_qPCR, dat_cov_all$cov19)

dat_cov_all$"Geimpft" <- factor(dat_cov_all$"Geimpft")


## -----------------------------------------------------------------------
## 2f  Treatment data and computer-assisted interviews
## -----------------------------------------------------------------------

## quarantine, outpatient, inpatient and intensive care treatment, from the
## Gutenberg Post-COVID Study screening interview
suppressWarnings(eval(parse(file = file.path(dir_sql, "GPS_read_SQL_1_S.r"))))
d_GPC_raw <- read_SQL_GPS(user = "<database user>", odbc = "<ODBC data source>")

channel <- odbcConnect("<ODBC data source>", uid = "<database user>")
info6   <- sqlQuery(channel, "select * from <symptom and treatment table>")
d_GPC   <- merge(d_GPC_raw, info6, by = "record_id", all.x = T,
                 suffixes = c("", "_screening"))

dat_cov_all <- merge(dat_cov_all,
                     d_GPC[, c("sid", "treat_treat_c___0", "treat_treat_c___1",
                               "treat_treat_c___2", "treat_inten_ynu")],
                     by = "sid", all.x = T)

## computer-assisted telephone interviews: severity category and symptoms
source(file.path(dir_sql, "GCS_read_SQL_CATI_try.R"))
d_cati <- read_SQL_CATI_GCS(reload = T, user = "<database user>")


################################################################################
## 3  BASELINE TABLE
################################################################################

## baseline13() is an in-house wrapper that tabulates the variables by group:
## mean and standard deviation for continuous variables, percentage and count
## for categorical ones, with a test for trend.

## one record per participant: follow-up where available, otherwise baseline
dat_BL_Tab <- dat_cov_all[which(dat_cov_all$lid %in%
                                  c(dat_all_norm_winz_GCC_BL$lid[
                                      which(!(dat_all_norm_winz_GCC_BL$lid_fu %in%
                                              dat_all_norm_winz_GCC_FU$lid))],
                                    dat_all_norm_winz_GCC_FU$lid)), ]

## post-infection samples only
dat_BL_Tab <- dat_BL_Tab[which(dat_BL_Tab$cov19 == "yes"), ]

dat_BL_Tab$Infektionsdatum_bin <- is.na(dat_BL_Tab$Infektionsdatum)

vars_tab <- c("age", "bmi",
              "nic", "hyper", "dyslip",
              "fhmis", "mi",
              "afib", "hemo",
              "egfr", "stroke",
              "cad", "pad",
              "hem", "cld",
              "copd", "cancer",
              "cvd", "auto",
              "diab", "adipos",
              "Severity.of.acute.infection",
              "Infektionsdatum_bin",
              "Sequelae.nach.3.Monaten", "Sequelae.nach.6.Monaten", "Geimpft")

nams_tab <- c("Age [y]", "BMI",
              "Smoking", "Arterial hypertension", "Dyslipidemia",
              "FH of MI/Stroke", "MI",
              "Atrial fibrillation", "Hemoglobin",
              "eGFR [ml/min/1.73m2]", "Stroke",
              "Coronary artery disease", "Peripheral artery disease",
              "Hematocrit", "Chronic liver disease",
              "Chronic obstructive pulmonary disease", "History of cancer",
              "Cardiovascular disease", "Autoimmune disease",
              "Diabetes", "Obesity",
              "Severity of infection",
              "Date of infection known",
              "Sequelae after 3 months", "Sequelae after 6 months",
              "Vaccination status")

cuts    <- c("age", "hyper", "afib", "Severity.of.acute.infection")
cutlabs <- c("Sociodemographics", "Traditional cardiovascular risk factors",
             "Clinical profile", "COVID-19 related characteristics")

BT1 <- baseline13(vars_tab, nams_tab, data = dat_BL_Tab, style = 3, group = "sex",
                  cuts = cuts, cutlabs = cutlabs, pv = F)


################################################################################
## 4  ENSEMBLE mRMR FEATURE SELECTION AND SELECTION STABILITY
################################################################################

set.thread.count(2)

## feature pool: all relative protein values; target: asymptomatic versus
## symptomatic infection
df <- dat_cov_all[which(!is.na(dat_cov_all$"asym_smod_and_hard")),
                  c("asym_smod_and_hard", vars_prot_rel)]
dd <- mRMR.data(data = df)


## ---- stability of a single mRMR run ------------------------------------
iterations <- 100
sel_100       <- selection_frequency(dd, iterations)
rank_100      <- sel_100$rank
ensemble_list_100 <- sel_100$ensemble_list

n_features <- 649
n_runs     <- ncol(ensemble_list_100)

selection_matrix <- matrix(FALSE, nrow = n_features, ncol = n_runs)
for (j in 1:n_runs) {
  selected <- ensemble_list_100[, j]
  selection_matrix[selected - 1, j] <- TRUE
}

stab_base_method <- getStability(t(selection_matrix))

## ---- stability of the ensemble itself ----------------------------------
## 100 repetitions of the 100-iteration ensemble; a protein counts as
## selected when it appears in more than half of the iterations
rep       <- 100
ranks_100 <- list()
for (k in 1:rep) {
  ensemble <- mRMR.ensemble(data = dd, target_indices = c(1),
                            solution_count = iterations, feature_count = 50)
  ensemble_list <- solutions(ensemble, with_fixed_features = T)[[1]]
  rank <- rep(0, length(c("asym_smod_and_hard", vars_prot_rel)))
  for (i in 1:iterations) {
    for (j in 1:length(c(vars_prot_rel))) {
      if (j %in% ensemble_list[, i]) rank[j] <- rank[j] + 1 / iterations
    }
  }
  ranks_100[[k]] <- rank
}

selection_matrix <- do.call(cbind, lapply(ranks_100, function(x) x > 0.5))
stab_ensemble    <- getStability(t(selection_matrix))

## ---- selection frequency across numbers of iterations ------------------
rank_5    <- selection_frequency(dd,    5)$rank
rank_20   <- selection_frequency(dd,   20)$rank
rank_50   <- selection_frequency(dd,   50)$rank
rank_200  <- selection_frequency(dd,  200)$rank
rank_500  <- selection_frequency(dd,  500)$rank
rank_1000 <- selection_frequency(dd, 1000)$rank

## unique gene symbols; the first occurrence keeps the plain symbol
protein_names <- make.unique(nams_genes_prot_rel, sep = "_")
protein_names[which(protein_names == "PROS1_1")] <- "PROS1_j"
protein_names[which(protein_names == "PROS1")]   <- "PROS1_1"
protein_names[which(protein_names == "PROS1_j")] <- "PROS1"

results <- cbind(rank_5[-1], rank_20[-1], rank_50[-1], rank_100[-1],
                 rank_200[-1], rank_500[-1], rank_1000[-1])
colnames(results) <- c("rank_5", "rank_20", "rank_50", "rank_100",
                       "rank_200", "rank_500", "rank_1000")
rownames(results) <- protein_names

result_long <- as.data.frame(results)
result_long$protein <- protein_names
result_long$feature <- seq_len(nrow(result_long))

result_long <- result_long %>%
  pivot_longer(cols = starts_with("rank_"),
               names_to = "bootstraps", values_to = "stability") %>%
  mutate(bootstraps = as.numeric(sub("rank_", "", bootstraps))) %>%
  mutate(stability = pmin(stability, 1))

all_bootstraps <- c(5, 20, 50, 100, 200, 500, 1000)

first_appearance <- result_long %>%
  filter(!is.na(stability)) %>%
  group_by(protein) %>%
  summarise(first_bootstrap = min(bootstraps), .groups = "drop")

zeros_before <- tidyr::expand_grid(first_appearance, bootstraps = all_bootstraps) %>%
  filter(bootstraps < first_bootstrap) %>%
  transmute(protein, bootstraps, stability = 0)

result_long <- bind_rows(result_long %>% dplyr::select(protein, bootstraps, stability),
                         zeros_before) %>%
  arrange(protein, bootstraps)

## proteins that reach the 50% threshold and later fall below it
drops <- result_long %>%
  arrange(protein, bootstraps) %>%
  group_by(protein) %>%
  mutate(previous_max  = lag(cummax(stability), default = 0),
         dropped_below = previous_max >= 0.5 & stability < 0.5) %>%
  filter(dropped_below) %>%
  slice(1) %>%
  ungroup()

labels_end <- result_long %>%
  filter(protein %in% drops$protein, bootstraps == 1000)

stability_path <- ggplot(result_long, aes(x = bootstraps, y = stability, group = protein)) +
  geom_line(alpha = 0.2) +
  geom_line(data = result_long %>% filter(protein %in% drops$protein),
            aes(group = protein), color = "red") +
  geom_text(data = labels_end, aes(x = 1000, y = stability, label = protein),
            color = "red", hjust = -0.05, size = 3) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  annotate("label", x = 1000, y = 1, hjust = 0, vjust = 1, fill = "white",
           label = paste0("Nogueira stability\n", "Ensemble (100 iterations)\n",
                          round(stab_ensemble$stability, 3), "\n", "95% CI: ",
                          round(stab_ensemble$lower, 3), "-",
                          round(stab_ensemble$upper, 3), "\n\n",
                          "Nogueira stability\n", "Base Method\n",
                          round(stab_base_method$stability, 3), "\n", "95% CI: ",
                          round(stab_base_method$lower, 3), "-",
                          round(stab_base_method$upper, 3))) +
  scale_x_continuous(trans = "log10", breaks = all_bootstraps, limits = c(5, 1000),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Number of bootstrap iterations", y = "Selection frequency",
       title = "Stability Path") +
  coord_cartesian(clip = "off") +
  theme_minimal() +
  theme(plot.margin = margin(t = 5.5, r = 160, b = 5.5, l = 5.5))

ggsave(file.path(out, "stability_path.png"), stability_path,
       width = 12, height = 8)


################################################################################
## 5  RIDGE REGRESSION AND THE PROTEIN SCORE
################################################################################

## proteins selected in more than half of the 100 ensemble iterations
sel      <- selection_frequency(dd, 100)
rank     <- sel$rank
vars_new <- featureNames(sel$ensemble)[which(rank > 0.5)]
nams_new <- nams_prot_rel[c(which(rank > 0.5) - 1)]

## Ridge regression (alpha = 0) of the selected proteins on asymptomatic
## versus symptomatic infection, with age and sex left unpenalised, and
## 10-fold cross-validation for the penalty parameter.
##
## E_Net1_RM() is a wrapper around E_Net1(), in section 1a.
Elastic_Net_covid_prot <- E_Net1_RM(v_dep = "asym_smod_and_hard",
                                    n_dep = "asym_smod_and_hard",
                                    v_pen = vars_new, n_pen = nams_new,
                                    v_adj = c("sex", "age"),
                                    n_adj = c("Sex", "Age"),
                                    dat = dat_cov_all, family = "logistic",
                                    nfolds = 10, parallel = 12, alpha = 0)

vars_sel_prot <- vars_new
nams_sel_prot <- nams_new

## discrimination before and after cross-validation
TAB <- data.frame(simpleAUC = c(Elastic_Net_covid_prot[[5]]),
                  cvAUC     = c(Elastic_Net_covid_prot[[6]]),
                  selected  = c(length(Elastic_Net_covid_prot[[7]])))
print(TAB)

## coefficient table of the selected proteins
print(Elastic_Net_covid_prot[[2]])

## cross-validation curve
plot_enet(Elastic_Net_covid_prot[-1], "Ensemble_proteins_severity_asym_vs_mod_and_hard")

## ---- protein score ------------------------------------------------------
## weighted sum of the standardised protein levels, using the ridge
## coefficients and dropping the intercept and the two covariates
cov_prot_score_ensemble <- as.matrix(coef(Elastic_Net_covid_prot[[3]],
                                          s = Elastic_Net_covid_prot[[11]]))
cov_prot_score_ensemble <- as.matrix(cov_prot_score_ensemble[-c(1:3), ])

dat_cov_all$cov_prot_score_ensemble <-
  scale(scale(as.matrix(dat_cov_all[, c(rownames(cov_prot_score_ensemble))])) %*%
        cov_prot_score_ensemble)

## ---- distribution of the score ------------------------------------------

score_plot <- function(df, x, filename) {
  p <- ggplot(df, aes(y = cov_prot_score_ensemble, x = .data[[x]], fill = .data[[x]])) +
    stat_slab(aes(thickness = after_stat(pdf * n)), scale = 0.7) +
    stat_dotsinterval(side = "bottom", scale = 0.5, slab_size = NA,
                      dotsize = 0.75, thickness = 4) +
    scale_fill_brewer(palette = "Set2") +
    theme_light()
  ggsave(file.path(out, filename), p, width = 10, height = 8)
  p
}

score_plot(dat_cov_all, "cov19", "score_by_infection.png")
score_plot(dat_cov_all[which(!is.na(dat_cov_all$Severity.of.acute.infection)), ],
           "Severity.of.acute.infection", "score_by_severity.png")
score_plot(dat_cov_all[which(!is.na(dat_cov_all$Sequelae.nach.3.Monaten)), ],
           "Sequelae.nach.3.Monaten", "score_by_sequelae_3m.png")
score_plot(dat_cov_all[which(!is.na(dat_cov_all$Sequelae.nach.6.Monaten)), ],
           "Sequelae.nach.6.Monaten", "score_by_sequelae_6m.png")

## The score is written out here and read back in for the analyses reported
## in the main figures and for the revision analyses below, so that the
## selection and the ridge fit do not have to be repeated.
write.csv2(dat_cov_all[, c("lid", "sid", "cov_prot_score_ensemble")],
           file.path(out, "ensemble_feature_selection_scores.csv"))


## ---- reload of the published score -------------------------------------
## the first submission have to be rerun for the revision.

set.thread.count(2)

dat_cov_all$asym_smod_and_hard <- NA
dat_cov_all$asym_smod_and_hard[which(dat_cov_all$Severity.of.acute.infection == "asymptomatisch")] <- "yes"
# dat_cov_all$asym_smod_and_hard[which(dat_cov_all$Severity.of.acute.infection == 'milder Verlauf')] <- 'no'
dat_cov_all$asym_smod_and_hard[which(dat_cov_all$Severity.of.acute.infection == "maessiger Verlauf")] <- "no"
dat_cov_all$asym_smod_and_hard[which(dat_cov_all$Severity.of.acute.infection == "schwerer Verlauf")]  <- "no"
dat_cov_all$asym_smod_and_hard <- 2 - as.numeric(as.factor(dat_cov_all$asym_smod_and_hard))

score <- read.csv2(file.path(dir_proj, "ensemble_feature_selection_scores_2025-08-26.csv"))
dat_cov_all <- merge(dat_cov_all, score[, c("lid", "cov_prot_score_ensemble")],
                     by = "lid", all.x = T)


################################################################################
## 6  PROTEIN SIGNATURE: RIDGE COEFFICIENTS (Figure 2A)
################################################################################

sig <- read.csv2(f_sig)
sig$Log_Odds_Ratio <- log(sig$OR.per.SD)
sig <- sig[order(sig$Log_Odds_Ratio, decreasing = FALSE), ]

signature_panel <- ggplot(sig, aes(x = Log_Odds_Ratio,
                                   y = factor(gene_name, levels = sig$gene_name),
                                   fill = Direction)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("+" = "#c98fa4", "-" = "#8eb56e")) +
  xlab("Log Odds Ratio") + ylab("Gene name") +
  theme_classic() +
  theme(axis.text  = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.position = "none",
        panel.grid = element_blank())

ggsave(file.path(out, "signature_log_odds_ratios.png"), signature_panel,
       width = 13, height = 11)


################################################################################
## 7  PROTEIN SCORE BY INFECTION SEVERITY
##    post-infection: Supplementary Figure 3A; pre-infection: Figure 3B
################################################################################

plot_score_by_severity <- function(timepoint, adjusted, title, fills, filename) {

  dat <- na.omit(get_scores(timepoint))
  dat$sex <- as.factor(dat$sex)

  if (adjusted) {
    dat$y    <- residuals(lm(cov_prot_score_ensemble ~ age + sex, data = dat))
    y_label  <- "Protein Score (adjusted for age and sex)"
  } else {
    dat$y    <- dat$cov_prot_score_ensemble
    y_label  <- "Protein Score for Symptomatic Infection"
  }

  dat$cov_severity_a <- factor(dat$cov_severity_a, levels = 1:4, labels = severity_labels)

  p <- ggplot(dat, aes(x = cov_severity_a, y = y,
                       color = cov_severity_a, fill = cov_severity_a)) +
    sm_raincloud() +
    xlab("Course of SARS-CoV-2 Infection") + ylab(y_label) +
    scale_fill_manual(values = fills) +
    scale_color_manual(values = fills) +
    stat_summary(fun.data = function(x) data.frame(y = Inf, label = paste0("n = ", length(x))),
                 geom = "text", hjust = 0.5, vjust = 1.5, col = "#262222", size = 5) +
    ylim(-3.5, 5.5) +
    theme_bw() +
    theme(axis.text  = element_text(size = 18),
          axis.title = element_text(size = 18),
          plot.title = element_text(size = 20),
          panel.grid.minor = element_blank()) +
    stat_compare_means(comparisons = list(c("Asymptomatic", "Mild"),
                                          c("Asymptomatic", "Moderate"),
                                          c("Asymptomatic", "Severe")),
                       method = "t.test", position = "identity", vjust = -0.2, size = 5) +
    ggtitle(title) +
    guides(fill = "none", color = "none")

  ggsave(file.path(out, filename), p, width = 12, height = 8)
  p
}

post_cols <- c("#cf5389", "#6379c2", "#a581c9", "#94ad5c")
pre_cols  <- c("#9c3e86", "#1e8fa6", "#714694", "#499c3e")

plot_score_by_severity("post", FALSE, "Post-infection samples", post_cols,
                       "post_infection_severity.png")
plot_score_by_severity("post", TRUE,  "Post-infection samples", post_cols,
                       "post_infection_severity_adjusted.png")
plot_score_by_severity("pre",  FALSE, "Pre-infection samples",  pre_cols,
                       "pre_infection_severity.png")
plot_score_by_severity("pre",  TRUE,  "Pre-infection samples",  pre_cols,
                       "pre_infection_severity_adjusted.png")


################################################################################
## 8  PRE- AND POST-INFECTION PROTEIN DYNAMICS
##    paired t-tests in the 189 individuals sampled at both timepoints
##    Supplementary Figure 4
################################################################################

prot_all <- get_proteins(timepoint = NULL)
gene_names <- read.csv2(f_sig)$gene_symbol

paired_sids <- intersect(prot_all$sid[prot_all$pre_post_infection == "pre"],
                         prot_all$sid[prot_all$pre_post_infection == "post"])
prot_paired <- prot_all %>% filter(sid %in% paired_sids)

pre_df  <- prot_paired %>% filter(pre_post_infection == "pre")  %>% arrange(sid)
post_df <- prot_paired %>% filter(pre_post_infection == "post") %>% arrange(sid)
stopifnot(all(pre_df$sid == post_df$sid))

paired_res <- data.frame(
  gene_name = gene_names,
  p_value   = sapply(gene_names, function(g) t.test(pre_df[[g]], post_df[[g]], paired = TRUE)$p.value),
  row.names = NULL
)
paired_res$fdr <- p.adjust(paired_res$p_value, method = "fdr")
paired_res <- paired_res[order(paired_res$p_value), ]

sig_proteins <- paired_res %>% filter(fdr < 0.05)
cat("Proteins tested:", nrow(paired_res),
    "| FDR < 0.05:", nrow(sig_proteins), "\n")

plot_df <- prot_paired %>%
  select(sid, pre_post_infection, all_of(sig_proteins$gene_name)) %>%
  pivot_longer(cols = all_of(sig_proteins$gene_name),
               names_to = "gene_name", values_to = "expression") %>%
  mutate(pre_post_infection = factor(pre_post_infection, levels = c("pre", "post")))

fdr_lab <- sig_proteins %>% mutate(lab = paste0("FDR p = ", signif(fdr, 2)))

pre_post_proteins <- ggplot(plot_df, aes(x = pre_post_infection, y = expression,
                                         fill = pre_post_infection, color = pre_post_infection)) +
  geom_violin(trim = FALSE, alpha = 1, size = 1, show.legend = FALSE) +
  geom_boxplot(width = 0.25, outlier.shape = NA, size = 1) +
  geom_beeswarm(dodge.width = 1, alpha = 0.5, cex = 1.5, priority = "ascending") +
  scale_fill_manual(values = c("pre" = "#f5dbb3", "post" = "#f5cbdf")) +
  scale_color_manual(values = c("pre" = "#db8704", "post" = "#d11d73")) +
  scale_x_discrete(labels = c("pre" = "Pre-infection", "post" = "Post-infection")) +
  xlab("") + ylab("Protein Expression") +
  facet_wrap(~ gene_name, scales = "free_y") +
  geom_text(data = fdr_lab, aes(x = 1.5, y = Inf, label = lab),
            inherit.aes = FALSE, vjust = 1.8, size = 5, colour = "black") +
  theme_bw() +
  theme(axis.text = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 15),
        strip.text = element_text(size = 14, face = "bold"),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

ggsave(file.path(out, "pre_post_proteins.png"), pre_post_proteins,
       width = 18, height = 7, dpi = 150)


################################################################################
## 9  PROTEIN SCORE AND ACUTE SYMPTOMS
##    4a  protein score ~ symptom (Figure 3A)
##    4b  individual proteins ~ symptom (Figure 4, heatmap)
################################################################################

## ---- 4a  protein score ~ symptom ---------------------------------------

symp <- merge(get_scores("post", c("age", "sex")),
              d_cati[c("sid", symp_vars)], by = "sid", all.x = TRUE)
names(symp)[match(symp_vars, names(symp))] <- symp_names
symp <- recode_symptoms(symp, symp_names)

## symptoms reported by at least 5% of the post-infection sample
counts <- sapply(symp[symp_names], function(x) sum(x == 1, na.rm = TRUE))
reduced_symps <- symp_names[counts >= 24]

res_symptoms <- score_on(symp, reduced_symps)

forplo(res_symptoms[, c("beta", "LCI", "UCI")],
       row.labels = gsub("_", " ", res_symptoms$predictor),
       em = "\u03b2", linreg = TRUE,
       col = "black", ci.edge = FALSE, left.bar = FALSE,
       left.align = TRUE, sort = TRUE, size = 1.5,
       margin.left = 10, margin.right = 25, column.spacing = 4,
       save = TRUE, save.path = out,
       save.name = "forplo_symptoms_score", save.type = "png",
       save.width = 10, save.height = 10)

## ---- 4b  individual proteins ~ symptom ---------------------------------

prot_post <- get_proteins("post")
df_prot_symp <- merge(prot_post[, c("sid", gene_names)],
                      symp[, c("sid", "age", "sex", reduced_symps)],
                      by = "sid")

estmat  <- matrix(NA, nrow = length(gene_names), ncol = length(reduced_symps),
                  dimnames = list(gene_names, reduced_symps))
pvalmat <- estmat

for (i in seq_along(gene_names)) {
  for (j in seq_along(reduced_symps)) {
    m  <- glm(reformulate(c(paste0("scale(`", gene_names[i], "`)"), "age", "sex"),
                          response = reduced_symps[j]),
              family = binomial, data = df_prot_symp)
    co <- summary(m)$coefficients[2, ]
    estmat[i, j]  <- co[1]
    pvalmat[i, j] <- co[4]
  }
}

## signed -log10(p), capped at 10 and set to zero below p = 0.05
signed_sig <- pmin(-log10(pvalmat), 10)
signed_sig[signed_sig < 1.3] <- 0
signed_sig <- signed_sig * sign(estmat)
colnames(signed_sig) <- gsub("_", " ", colnames(signed_sig))

heatmap_matrix <- t(signed_sig)

## label proteins and symptoms with many associations
sig_threshold <- 1.3
proteins_per_symptom <- colSums(abs(signed_sig) > sig_threshold, na.rm = TRUE)
symptoms_per_protein <- rowSums(abs(signed_sig) > sig_threshold, na.rm = TRUE)
cut_sym <- 10; cut_prot <- 10

row_hl <- proteins_per_symptom[rownames(heatmap_matrix)] >= cut_sym
col_hl <- symptoms_per_protein[colnames(heatmap_matrix)] >= cut_prot

row_labels <- ifelse(row_hl,
                     paste0(rownames(heatmap_matrix), " (",
                            proteins_per_symptom[rownames(heatmap_matrix)], ")"),
                     rownames(heatmap_matrix))
col_labels <- ifelse(col_hl,
                     paste0(colnames(heatmap_matrix), " (",
                            symptoms_per_protein[colnames(heatmap_matrix)], ")"),
                     colnames(heatmap_matrix))

ht <- Heatmap(
  heatmap_matrix,
  name = "-log10(p)",
  col  = colorRamp2(c(-3, 0, 3), c("#0439b5", "white", "#027302")),
  rect_gp = gpar(col = "grey", lwd = 0.5),
  row_labels = row_labels, column_labels = col_labels,
  row_names_gp = gpar(fontsize = 8, fontface = ifelse(row_hl, "bold", "plain"),
                      col = ifelse(row_hl, "#2c8f40", "black")),
  column_names_gp = gpar(fontsize = 8, fontface = ifelse(col_hl, "bold", "plain"),
                         col = ifelse(col_hl, "#2c8f40", "black"), rot = 60),
  row_gap = unit(2, "mm"), column_gap = unit(2, "mm")
)

png(file.path(out, "heatmap_symptoms_proteins.png"),
    width = 2400, height = 1800, res = 300)
draw(ht, heatmap_legend_side = "right")
dev.off()


################################################################################
## 10  SYMPTOM PREVALENCE THRESHOLD: SENSITIVITY ANALYSIS
##    5% (>= 24 individuals) vs 10% (>= 50 individuals)
##    Supplementary Figure 6
################################################################################

run_threshold <- function(min_n) {
  score_on(symp, symp_names[counts >= min_n]) %>%
    mutate(n_cases = counts[predictor])
}

thr5  <- 24
thr10 <- ceiling(0.10 * nrow(symp))

res5  <- run_threshold(thr5)
res10 <- run_threshold(thr10)

cat("Symptoms at 5% threshold: ", nrow(res5),
    "| at 10% threshold:", nrow(res10), "\n")
cat("p < 0.05 at 5%: ", sum(res5$p  < 0.05), "of", nrow(res5),  "\n")
cat("p < 0.05 at 10%:", sum(res10$p < 0.05), "of", nrow(res10), "\n")

dropped <- setdiff(res5$predictor, res10$predictor)
cat("Symptoms dropped at the 10% threshold:", length(dropped),
    "| of these with p < 0.05:", sum(res5$p[res5$predictor %in% dropped] < 0.05), "\n")

plot_dat <- res5 %>%
  mutate(retained = ifelse(n_cases >= thr10, 1, 2),
         label = gsub("_", " ", predictor)) %>%
  arrange(desc(beta))

forplo(plot_dat[, c("beta", "LCI", "UCI")],
       row.labels = plot_dat$label,
       em = "\u03b2", linreg = TRUE,
       ci.edge = FALSE, left.bar = FALSE, right.bar = FALSE,
       left.align = TRUE, sort = FALSE, size = 2,
       margin.left = 12, margin.right = 22, column.spacing = 4,
       fill.by = plot_dat$retained,
       fill.colors = c("black", "#c1121f"),
       fill.labs = c("Retained at 10% threshold",
                     "Excluded at 10% threshold, included at 5%"),
       legend = TRUE, legend.spacing = 2, legend.vadj = 1, legend.hadj = -2.5,
       save = TRUE, save.path = out,
       save.name = "forplo_symptom_threshold_sensitivity",
       save.type = "png", save.width = 10, save.height = 10)


################################################################################
## 11  PROTEIN SCORE AND POST-COVID SEQUELAE (Figure 3B)
##    prevalence ratios from Poisson regression with robust standard errors
################################################################################

covid_data <- read.csv(f_covid)

pcs <- merge(get_scores("post", c("age", "sex")),
             covid_data[c("ID", "Sequelae.nach.3.Monaten", "Sequelae.nach.6.Monaten")],
             by.x = "sid", by.y = "ID", all.x = TRUE)

pcs$pcs3 <- ifelse(pcs$Sequelae.nach.3.Monaten == "yes", 1,
                   ifelse(pcs$Sequelae.nach.3.Monaten == "no", 0, NA))
pcs$pcs6 <- ifelse(pcs$Sequelae.nach.6.Monaten == "yes", 1,
                   ifelse(pcs$Sequelae.nach.6.Monaten == "no", 0, NA))

table(pcs$pcs3, useNA = "ifany")
table(pcs$pcs6, useNA = "ifany")

res_pcs <- rbind(
  prevalence_ratio("pcs3", pcs, label = "Post-COVID sequelae > 3 months"),
  prevalence_ratio("pcs6", pcs, label = "Post-COVID sequelae > 6 months")
)
print(res_pcs)

forplo(res_pcs[, c("PR", "LCI", "UCI")],
       row.labels = res_pcs$outcome,
       em = "PR", linreg = FALSE,
       col = "#80644f", ci.edge = FALSE,
       left.bar = FALSE, right.bar = FALSE, left.align = TRUE,
       xlim = c(0.9, 1.5), margin.left = 12, margin.right = 10,
       legend = FALSE,
       save = TRUE, save.path = out,
       save.name = "forplo_post_covid_sequelae", save.type = "png",
       save.width = 7, save.height = 4)


################################################################################
## 12  PROTEIN SCORE AND LONG-TERM SYMPTOMS (> 6 months)
##    Supplementary Figure 9
################################################################################

seq_data <- read.csv2(f_seq)
long_symp <- merge(get_scores("post", c("age", "sex")), seq_data,
                   by.x = "sid", by.y = "X", all.x = TRUE)
names(long_symp)[9:69] <- symp_names
long_symp <- recode_symptoms(long_symp, symp_names)

long_counts <- sapply(long_symp[symp_names], function(x) sum(x == 1, na.rm = TRUE))
long_reduced <- symp_names[long_counts >= 24]

res_long <- lapply(long_reduced, function(s) {
  m  <- glm(reformulate(c("cov_prot_score_ensemble", "age", "sex"), response = s),
            data = long_symp, family = binomial)
  co <- summary(m)$coefficients[2, ]
  data.frame(symptom = s, OR = exp(co[1]),
             LCI = exp(co[1] - 1.96 * co[2]),
             UCI = exp(co[1] + 1.96 * co[2]),
             p = co[4], row.names = NULL)
}) %>% bind_rows()

forplo(res_long[, c("OR", "LCI", "UCI")],
       row.labels = gsub("_", " ", res_long$symptom),
       col = "black", ci.edge = FALSE, left.bar = FALSE,
       left.align = TRUE, sort = TRUE, size = 1.5,
       margin.left = 10, margin.right = 25, column.spacing = 4,
       save = TRUE, save.path = out,
       save.name = "forplo_long_term_symptoms", save.type = "png",
       save.width = 10, save.height = 10)


################################################################################
## 13  PROTEIN SCORE AND CARDIOVASCULAR RISK FACTORS
##    (the analysis of lipid markers measured before infection, Figure 5,
##     was performed separately and is not part of this script)
################################################################################

cvrf_vars <- c("hyper", "diab", "dyslip", "cancer", "mi", "adipos", "nic",
               "stroke", "afib", "cad", "chf", "ckd", "cld", "copd", "cvd",
               "pad", "pulmd", "vte")

cvrf <- merge(get_scores("post", c("age", "sex")),
              d[c("sid", cvrf_vars)], by = "sid", all.x = TRUE)
cvrf[cvrf_vars] <- lapply(cvrf[cvrf_vars], function(x) {
  x[x == 99] <- NA
  ifelse(x == "yes", 1, ifelse(is.na(x), NA, 0))
})

res_cvrf <- score_on(cvrf, cvrf_vars)

cvrf_labels <- c("Arterial Hypertension", "Diabetes Mellitus Type 2", "Dyslipidemia",
                 "History of Cancer", "History of Myocardial Infarction", "Obesity",
                 "Smoking", "History of Stroke", "Atrial Fibrillation",
                 "Coronary Artery Disease", "Chronic Heart Failure",
                 "Chronic Kidney Disease", "Chronic Liver Disease",
                 "Chronic Obstructive Pulmonary Disease", "Cardiovascular Disease",
                 "Peripheral Artery Disease", "Pulmonary Disease",
                 "History of Venous Thromboembolism")

forplo(res_cvrf[, c("beta", "LCI", "UCI")],
       row.labels = cvrf_labels,
       em = "\u03b2", linreg = TRUE,
       col = "black", ci.edge = FALSE, left.bar = FALSE,
       left.align = TRUE, sort = TRUE,
       margin.left = 10, margin.right = 25, column.spacing = 4,
       save = TRUE, save.path = out,
       save.name = "forplo_score_cvrf", save.type = "png",
       save.width = 12, save.height = 6)


################################################################################
## 14  PER-PROTEIN STATISTICS ACROSS SEVERITY GROUPS
##    group means relative to asymptomatic (Supplementary Figure 5)
##    per-protein table (Supplementary Figure 10)
################################################################################

ens <- merge(get_proteins("post"), d_cati[c("sid", "cov_severity_a")],
             by = "sid", all.x = TRUE)

dat <- ens %>%
  filter(!is.na(cov_severity_a)) %>%
  mutate(severity = factor(cov_severity_a, levels = 1:4, labels = severity_labels))

X  <- apply(as.matrix(dat[, gene_names]), 2, as.numeric)
Xz <- scale(X)

cat("N with severity category:", nrow(dat), "\n"); print(table(dat$severity))

## ---- group means, expressed as the difference from asymptomatic --------

grp_mean <- as.data.frame(Xz) %>%
  mutate(severity = dat$severity) %>%
  group_by(severity) %>%
  summarise(across(all_of(gene_names), \(x) mean(x, na.rm = TRUE)), .groups = "drop")

asymp_row <- grp_mean %>% filter(severity == "Asymptomatic") %>%
  select(all_of(gene_names)) %>% as.numeric()

mat_diff <- t(sweep(as.matrix(grp_mean[, gene_names]), 2, asymp_row, "-"))
colnames(mat_diff) <- as.character(grp_mean$severity)
mat_diff <- mat_diff[, colnames(mat_diff) != "Asymptomatic", drop = FALSE]

ht_diff <- Heatmap(
  mat_diff, name = "Z-score",
  col = colorRamp2(c(-0.8, 0, 0.8), c("#1d3557", "white", "#e63946")),
  cluster_columns = FALSE, cluster_rows = TRUE,
  column_names_rot = 45, column_names_gp = gpar(fontsize = 9),
  column_title = "Course of SARS-CoV-2 infection", column_title_side = "bottom",
  row_title = "Signature proteins", row_names_gp = gpar(fontsize = 8),
  width = unit(3.6, "cm"),
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.text(sprintf("%.2f", mat_diff[i, j]), x, y, gp = gpar(fontsize = 6.5, col = "grey20"))
  }
)

png(file.path(out, "heatmap_group_means_vs_asymptomatic.png"),
    width = 1400, height = 2600, res = 300)
draw(ht_diff, heatmap_legend_side = "right")
dev.off()

## ---- per-protein comparisons with BH adjustment -------------------------

long_df <- dat %>%
  mutate(severity_group = factor(case_when(
    cov_severity_a == 1 ~ "Asymptomatic",
    cov_severity_a == 2 ~ "Mild",
    cov_severity_a %in% c(3, 4) ~ "Moderate/Severe"),
    levels = c("Asymptomatic", "Mild", "Moderate/Severe"))) %>%
  pivot_longer(cols = all_of(gene_names), names_to = "gene", values_to = "expression") %>%
  group_by(gene) %>%
  mutate(expression_z = as.numeric(scale(expression))) %>%
  ungroup()

stats_tab <- long_df %>%
  group_by(gene) %>%
  summarise(
    n_asymp  = sum(severity_group == "Asymptomatic"    & !is.na(expression)),
    n_mild   = sum(severity_group == "Mild"            & !is.na(expression)),
    n_modsev = sum(severity_group == "Moderate/Severe" & !is.na(expression)),
    smd_mild   = mean(expression_z[severity_group == "Mild"], na.rm = TRUE) -
                 mean(expression_z[severity_group == "Asymptomatic"], na.rm = TRUE),
    smd_modsev = mean(expression_z[severity_group == "Moderate/Severe"], na.rm = TRUE) -
                 mean(expression_z[severity_group == "Asymptomatic"], na.rm = TRUE),
    p_mild   = t.test(expression[severity_group == "Asymptomatic"],
                      expression[severity_group == "Mild"])$p.value,
    p_modsev = t.test(expression[severity_group == "Asymptomatic"],
                      expression[severity_group == "Moderate/Severe"])$p.value,
    .groups = "drop") %>%
  mutate(padj_mild   = p.adjust(p_mild,   method = "BH"),
         padj_modsev = p.adjust(p_modsev, method = "BH"),
         min_padj    = pmin(padj_mild, padj_modsev, na.rm = TRUE))

cat("Proteins with BH-adjusted p < 0.05 in at least one comparison:",
    sum(stats_tab$min_padj < 0.05, na.rm = TRUE), "\n")

write.csv(stats_tab, file.path(out, "per_protein_statistics.csv"), row.names = FALSE)


################################################################################
## 15  ANALYSES ADDED DURING REVISION
################################################################################

## -----------------------------------------------------------------------
## 15a  Serologically confirmed infections
## -----------------------------------------------------------------------

##      nucleocapsid antibodies above the assay cut-off
##      (Abbott IgG >= 1.4 RLU or Roche >= 0.8 U/ml, baseline or follow-up)

sero <- merge(get_scores("post"), d[c("sid", "igg", "roche", "igg_fu", "roche_fu")],
              by = "sid", all.x = TRUE)

sero$any_measured <- !is.na(sero$igg) | !is.na(sero$roche) |
                     !is.na(sero$igg_fu) | !is.na(sero$roche_fu)

sero$sero_pos <- (sero$igg      >= 1.4) %in% TRUE |
                 (sero$roche    >= 0.8) %in% TRUE |
                 (sero$igg_fu   >= 1.4) %in% TRUE |
                 (sero$roche_fu >= 0.8) %in% TRUE

table(sero$any_measured)
table(sero$sero_pos, useNA = "ifany")
table(sero$cov_severity_a[sero$sero_pos], useNA = "ifany")

sero_sub <- subset(sero, sero_pos & !is.na(cov_severity_a))
sero_sub$symptomatic <- ifelse(sero_sub$cov_severity_a == 1, 0, 1)

cat("N in serologically confirmed subgroup:", nrow(sero_sub), "\n")
print(table(sero_sub$cov_severity_a))

summary(lm(cov_prot_score_ensemble ~ symptomatic + age + sex, data = sero_sub))
summary(lm(cov_prot_score_ensemble ~ factor(cov_severity_a) + age + sex, data = sero_sub))


## -----------------------------------------------------------------------
## 15b  Baseline characteristics: individuals with vs without proteomics
## -----------------------------------------------------------------------

prot_ids <- unique(read.csv2(f_scores)$sid)

base_vars <- c("age", "sex", "bmi", "nic", "hyper", "diab", "dyslip",
               "adipos", "cancer", "mi", "stroke", "cvd")

base <- d[!duplicated(d$sid), c("sid", base_vars)]
base$group <- ifelse(base$sid %in% prot_ids, "proteomics", "no proteomics")
table(base$group)

bin_vars <- c("nic", "hyper", "diab", "dyslip", "adipos", "cancer", "mi", "stroke", "cvd")
base[bin_vars] <- lapply(base[bin_vars], function(x) {
  x[x %in% c(99, "unknown", "")] <- NA
  ifelse(x %in% c("yes", 1, "1"), 1, ifelse(is.na(x), NA, 0))
})
base$female <- ifelse(base$sex == "Women", 1, 0)

row_cont <- function(v, label) {
  a <- base[[v]][base$group == "proteomics"]
  b <- base[[v]][base$group == "no proteomics"]
  data.frame(Characteristic = label,
             Proteomics    = sprintf("%.1f \u00b1 %.1f", mean(a, na.rm = TRUE), sd(a, na.rm = TRUE)),
             No_proteomics = sprintf("%.1f \u00b1 %.1f", mean(b, na.rm = TRUE), sd(b, na.rm = TRUE)),
             p = signif(t.test(a, b)$p.value, 3))
}

row_bin <- function(v, label) {
  pct <- function(g) sprintf("%.1f (%d)",
                             100 * mean(base[[v]][base$group == g], na.rm = TRUE),
                             sum(base[[v]][base$group == g] == 1, na.rm = TRUE))
  data.frame(Characteristic = label,
             Proteomics = pct("proteomics"),
             No_proteomics = pct("no proteomics"),
             p = signif(chisq.test(table(base$group, base[[v]]))$p.value, 3))
}

baseline_table <- rbind(
  row_cont("age", "Age [years]"),
  row_bin("female", "Female sex"),
  row_cont("bmi", "BMI [kg/m2]"),
  row_bin("nic",    "Smoking"),
  row_bin("hyper",  "Arterial hypertension"),
  row_bin("dyslip", "Dyslipidemia"),
  row_bin("adipos", "Obesity"),
  row_bin("diab",   "Diabetes mellitus type 2"),
  row_bin("cvd",    "Cardiovascular disease"),
  row_bin("cancer", "History of cancer"),
  row_bin("mi",     "History of myocardial infarction"),
  row_bin("stroke", "History of stroke")
)

print(baseline_table)
write.xlsx(baseline_table, file.path(out, "baseline_comparison.xlsx"), rowNames = FALSE)


## -----------------------------------------------------------------------
## 15c  Analysis 1: sensitivity to the method of infection ascertainment
## -----------------------------------------------------------------------

## Restricting the sample to PCR-confirmed infections leaves too few
## asymptomatic infections for the comparison to be repeated.

table(dat_cov_all[which(dat_cov_all$pcr_qPCR == "yes"), "Severity.of.acute.infection"])


## -----------------------------------------------------------------------
## 15d  Analysis 2: objective severity definition
## -----------------------------------------------------------------------

## Treatment during the acute infection, from the screening interview of the
## Gutenberg Post-COVID Study.

table(dat_cov_all$treat_treat_c___0)   # quarantine
table(dat_cov_all$treat_treat_c___1)   # outpatient treatment
table(dat_cov_all$treat_treat_c___2)   # inpatient treatment
table(dat_cov_all$treat_inten_ynu)     # intensive care

## ---- protein score and hospitalisation ---------------------------------

dat <- read_excel(file.path(dir_proj, "revision_samples.xlsx"), sheet = 1)
dat <- merge(dat, d_GPC[, c("sid", "treat_treat_c___0", "treat_treat_c___1",
                            "treat_treat_c___2", "treat_inten_ynu")],
             by = "sid", all.x = T)

dat$hospitalised <- NA
dat$hospitalised[which(dat$treat_treat_c___2 == T | dat$treat_inten_ynu == T)] <- 1
dat$hospitalised[which(dat$treat_treat_c___2 == F)] <- 0

## logistic regression of hospitalisation on the protein score, age and sex
TAB <- reg_tab("hospitalised", "hospitalised",
               c("cov_prot_score_ensemble", "age", "sex"),
               c("Protein Score", "Age", "Sex"),
               data = dat, regression = "logistic")
print(TAB)


## -----------------------------------------------------------------------
## 15e  Analysis 3: extended confounder adjustment
## -----------------------------------------------------------------------

## post-infection samples, joined to the severity category and the symptoms
## of the acute infection from the telephone interview
dat_post <- dat_cov_all[which(dat_cov_all$pre_post_infection == "post"), ]

covid_symp_vars <- paste0("symp_ynu", 1:61)

dat_post <- merge(dat_post, d_cati[c("sid", "cov_severity_a", covid_symp_vars)],
                  by = "sid", all.x = TRUE)

covid_symp_nams <- c("Fatigue", "Malaise", "Fever", "Chills", "Cough", "Sore throat",
  "Hoarseness", "Runny nose", "Headache", "Shortness of breath",
  "Pain when breathing", "Wheezing", "Chest pain", "Tachycardia",
  "Olfactory disturbances", "Taste disturbances", "Visual disturbances",
  "Conjunctivitis", "Hearing problems", "Tinnitus", "Memory problems",
  "Confusion", "Concentration problems", "Sleep disturbances", "Less sleep",
  "More sleep", "Dizziness", "Fainting spells", "Joint pain or swelling",
  "Swollen ankles", "Muscle pain", "Muscle stiffness", "Body aches",
  "Twiching of the limbs", "Weakness of the limbs", "Tremor",
  "Numbness or tingling", "Mood swings", "Depression",
  "Loss of interest or pleasure", "Anxiety", "Abdominal pain", "Diarrhea",
  "Constipation", "Nausea or vomiting", "Loss of appetite", "Weight loss",
  "Skin rash", "Hair loss", "Seizures", "Problems urinating",
  "Problems swallowing", "Problems with balance",
  "Trouble walking or falling", "Slow movements",
  "Painful menstrual bleeding", "Erectile dysfunction", "Behavioural changes",
  "Hallucinations", "COVID toes", "Partial face or body paralysis")

vars_pred <- "cov_prot_score_ensemble"
nams_pred <- "Protein Score"

## ---- the three adjustment sets -----------------------------------------
vars_adj <- c("sex", "age",
              "bmi", "nic",
              ## comorbidities
              "hyper", "diab", "dyslip", "adipos", "cancer", "mi", "stroke",
              "cvd", "cad", "cld", "copd", "pad", "afib",
              ## vaccination
              "Geimpft")

nams_adj <- c("Sex", "Age",
              "BMI", "Current or former smoker",
              ## comorbidities
              "Hypertension", "Diabetes", "Dyslipidemia", "Obesity", "Cancer",
              "Myocardial infarction", "Stroke", "Cardiovascular disease",
              "Coronary artery disease", "Chronic liver disease",
              "Chronic obstructive pulmonary disease",
              "Peripheral artery disease", "Atrial fibrillation",
              ## vaccination
              "Vaccination status")

## Model 1: age and sex
## Model 2: Model 1 plus BMI and smoking
## Model 3: Model 2 plus comorbidities and vaccination before infection
vars_adj_list <- list(vars_adj[1:2], vars_adj[1:4], vars_adj)
nams_adj_list <- list(nams_adj[1:2], nams_adj[1:4], nams_adj)

## symptoms reported by fewer than 24 participants (5% of the post-infection
## sample) are not analysed
symptom_counts_min <- 24

out_symp <- c()
for (i in 1:length(covid_symp_vars)) {
  dat_post[which(dat_post[, covid_symp_vars[i]] == 99), covid_symp_vars[i]] <- NA
  dat_post[, covid_symp_vars[i]] <- factor(dat_post[, covid_symp_vars[i]],
                                           levels = c(0, 1),
                                           labels = c("no", "yes"))

  if (length(which(dat_post[, covid_symp_vars[i]] == "yes")) < symptom_counts_min) {
    out_symp <- c(out_symp, i)
  }
}
if (length(out_symp) > 0) {
  covid_symp_vars <- covid_symp_vars[-out_symp]
  covid_symp_nams <- covid_symp_nams[-out_symp]
}

## missing values per symptom
TAB1 <- TAB_miss(covid_symp_vars, covid_symp_nams, 494, dat_post)
print(TAB1)


## -----------------------------------------------------------------------
## Protein score by severity category, three adjustment sets
##   The score is regressed on the covariates and the residuals are
##   plotted, so that the group comparison is adjusted.
## -----------------------------------------------------------------------

form_adj_1 <- as.formula(paste0("cov_prot_score_ensemble ~ ",
                                paste0(vars_adj_list[[1]], collapse = " + ")))
model_adj_1 <- lm(form_adj_1, data = dat_post)
dat_post$resid_adjusted_1 <- residuals(model_adj_1)

form_adj_2 <- as.formula(paste0("cov_prot_score_ensemble ~ ",
                                paste0(vars_adj_list[[2]], collapse = " + ")))
model_adj_2 <- lm(form_adj_2, data = dat_post, na.action = na.exclude)
dat_post$resid_adjusted_2 <- residuals(model_adj_2)

form_adj_3 <- as.formula(paste0("cov_prot_score_ensemble ~ ",
                                paste0(vars_adj_list[[3]][-which(vars_adj_list[[3]] == "cld")],
                                       collapse = " + ")))
model_adj_3 <- lm(form_adj_3, data = dat_post, na.action = na.exclude)
dat_post$resid_adjusted_3 <- residuals(model_adj_3)

dat_post$cov_severity_a <- factor(dat_post$cov_severity_a,
                                  levels = c(1, 2, 3, 4),
                                  labels = c("Asymptomatic", "Mild",
                                             "Moderate", "Severe"))
dat_post_plot <- dat_post[!is.na(dat_post$cov_severity_a), ]

plot_severity_adjusted <- function(resid_var, y_label, filename) {

  p <- ggplot(dat_post_plot, aes(x = cov_severity_a, y = .data[[resid_var]],
                                 color = cov_severity_a, fill = cov_severity_a)) +
    sm_raincloud() +
    xlab("Course of SARS-CoV-2 Infection") +
    ylab(y_label) +
    scale_fill_manual(values  = c("#cf5389", "#6379c2", "#a581c9", "#94ad5c")) +
    scale_color_manual(values = c("#cf5389", "#6379c2", "#a581c9", "#94ad5c")) +
    stat_summary(fun.data = function(x) data.frame(y = Inf, label = paste0("n = ", length(x))),
                 geom = "text", hjust = 0.5, vjust = 1.5, col = "#262222", size = 5) +
    ylim(-3.5, 7) +
    theme_bw() +
    theme(axis.text.x  = element_text(size = 18),
          axis.title.x = element_text(size = 18, margin = margin(t = 20)),
          axis.title.y = element_text(size = 18),
          axis.text.y  = element_text(size = 18),
          plot.title   = element_text(size = 20),
          panel.grid.minor = element_blank()) +
    stat_compare_means(comparisons = list(c("Asymptomatic", "Mild"),
                                          c("Asymptomatic", "Moderate"),
                                          c("Asymptomatic", "Severe")),
                       method = "t.test", position = "identity",
                       vjust = -0.2, size = 5) +
    labs(fill = "Group", color = "Group") +
    ggtitle("Post-infection samples") +
    guides(fill = "none")

  ggsave(plot = p, filename = file.path(out, filename), width = 12, height = 12)
  p
}

plot_severity_adjusted("resid_adjusted_1",
                       "Protein Score (adjusted for age and sex)",
                       "severity_adjusted_model1.png")
plot_severity_adjusted("resid_adjusted_2",
                       "Protein Score (adjusted for age, sex, BMI, smoking)",
                       "severity_adjusted_model2.png")
plot_severity_adjusted("resid_adjusted_3",
                       paste("Protein Score (adjusted for age, sex, BMI, smoking,",
                             "comorbidities, vaccination before infection)"),
                       "severity_adjusted_model3.png")


## -----------------------------------------------------------------------
## Post-COVID sequelae, three adjustment sets
##   prevalence ratios from Poisson regression with robust standard errors
## -----------------------------------------------------------------------

sequ <- read_excel(file.path(dir_proj, "pcs_sample_494.xlsx"))
dat_post <- merge(dat_post, sequ, by = "sid")

sequelae <- c("pcs3", "pcs6")
labs     <- c("Post-COVID sequelae > 3 months", "Post-COVID sequelae > 6 months")

vars_adj_list2 <- list(c(vars_pred, vars_adj[1:2]),
                       c(vars_pred, vars_adj[1:4]),
                       c(vars_pred, vars_adj))
nams_adj_list2 <- list(c(nams_pred, nams_adj[1:2]),
                       c(nams_pred, nams_adj[1:4]),
                       c(nams_pred, nams_adj))

TAB1 <- reg_tab_list(sequelae[1], labs[1], vars_adj_list2, nams_adj_list2,
                     data = dat_post, regression = "robust.poisson")
TAB2 <- reg_tab_list(sequelae[2], labs[2], vars_adj_list2, nams_adj_list2,
                     data = dat_post, regression = "robust.poisson")

## the row of the protein score from each of the three models, side by side
TAB <- matrix("", ncol = 24, nrow = 2)
TAB[1, ] <- c("", TAB1[2, 2:5],   TAB1[22:24, 1],
              "", TAB1[2, 7:10],  TAB1[22:24, 6],
              "", TAB1[2, 12:15], TAB1[22:24, 11])
TAB[2, ] <- c("", TAB2[2, 2:5],   TAB2[22:24, 1],
              "", TAB2[2, 7:10],  TAB2[22:24, 6],
              "", TAB2[2, 12:15], TAB2[22:24, 11])
colnames(TAB) <- c("Model 1:", "PR", "L 95% CI", "U 95% CI", "p-value", "AUC", "N", "Events",
                   "Model 2:", "PR", "L 95% CI", "U 95% CI", "p-value", "AUC", "N", "Events",
                   "Model 3:", "PR", "L 95% CI", "U 95% CI", "p-value", "AUC", "N", "Events")
rownames(TAB) <- labs
TAB_pcs_pr <- TAB
print(TAB_pcs_pr)


## -----------------------------------------------------------------------
## Symptoms of the acute infection, three adjustment sets
##   linear regression of the protein score on each symptom
## -----------------------------------------------------------------------

result_1 <- reg_tab_loop(vars_dep = vars_pred, nams_dep = nams_pred,
                         vars_pred = covid_symp_vars, nams_pred = covid_symp_nams,
                         vars_adj = c(vars_adj_list[[1]]), nams_adj = c(nams_adj_list[[1]]),
                         data = dat_post, regression = "linear",
                         label_n = F, sort_p = F)
result_2 <- reg_tab_loop(vars_dep = vars_pred, nams_dep = nams_pred,
                         vars_pred = covid_symp_vars, nams_pred = covid_symp_nams,
                         vars_adj = c(vars_adj_list[[2]]), nams_adj = c(nams_adj_list[[2]]),
                         data = dat_post, regression = "linear",
                         label_n = F, sort_p = F)
result_3 <- reg_tab_loop(vars_dep = vars_pred, nams_dep = nams_pred,
                         vars_pred = covid_symp_vars, nams_pred = covid_symp_nams,
                         vars_adj = c(vars_adj_list[[3]]), nams_adj = c(nams_adj_list[[3]]),
                         data = dat_post, regression = "linear",
                         label_n = F, sort_p = F)

TAB_symp <- cbind("", result_1, "", result_2, "", result_3)
colnames(TAB_symp) <- c("Model 1:", "beta", "L 95% CI", "U 95% CI", "p-value", "R2", "N",
                        "Model 2:", "beta", "L 95% CI", "U 95% CI", "p-value", "R2", "N",
                        "Model 3:", "beta", "L 95% CI", "U 95% CI", "p-value", "R2", "N")
print(TAB_symp)


## -----------------------------------------------------------------------
## 15f  Forest plots for Analysis 3
##   Input: the estimates from section 15e, in
##   forest_symptoms_data.csv and forest_pcs_data.csv
## -----------------------------------------------------------------------

adjust <- c("Model 1" = "adjusted for age and sex",
            "Model 2" = "adjusted for age, sex, BMI and smoking",
            "Model 3" = "adjusted for age, sex, BMI, smoking, comorbidities and vaccination")

## shorter versions used as row labels in the post-COVID plot
adjust_short <- c("Model 1" = "Model 1: age, sex",
                  "Model 2" = "Model 2: Model 1 + BMI, smoking",
                  "Model 3" = "Model 3: Model 2 + comorbidities, vaccination")

cols <- c("Model 1" = "#bd4f7f", "Model 2" = "#3a6e9c", "Model 3" = "#70913f")

## ---- acute symptoms: one plot per adjustment set -----------------------

sym <- read.csv("forest_symptoms_data.csv")

## keep the ordering of the first model for all three plots
first <- sym[sym$model == "Model 1", ]
ord   <- first$symptom[order(-first$beta)]

for (m in names(adjust)) {

  d <- sym[sym$model == m, ]
  d <- d[match(ord, d$symptom), ]

  forplo(d[, c("beta", "lci", "uci")],
         em          = "β",
         linreg      = TRUE,
         row.labels  = d$symptom,
         col         = cols[[m]],
         ci.sep      = ";",
         ci.edge     = FALSE,
         left.bar    = FALSE,
         right.bar   = FALSE,
         horiz.bar   = FALSE,
         left.align  = TRUE,
         sort        = FALSE,
         size        = 1.5,
         margin.left = 12,
         margin.right = 22,
         column.spacing = 4,
         legend      = FALSE,
         font        = "Arial",
         save        = TRUE,
         save.path   = out,
         save.name   = paste0("forest_symptoms_", tolower(sub(" ", "", m))),
         save.type   = "png",
         save.width  = 10,
         save.height = 10)

  cat(m, "-", adjust[[m]], "\n")
}

## ---- post-COVID sequelae: adjustment sets as row labels ----------------

pcs <- read.csv("forest_pcs_data.csv")

pcs$outcome <- factor(pcs$outcome, levels = unique(pcs$outcome))
pcs <- pcs[order(pcs$outcome, pcs$model), ]

n_o <- nlevels(pcs$outcome)

forplo(pcs[, c("pr", "lci", "uci")],
       em          = "PR",
       linreg      = FALSE,
       groups      = rep(1:n_o, each = 3),
       grouplabs   = levels(pcs$outcome),
       row.labels  = as.character(adjust_short[pcs$model]),
       fill.by     = rep(c(1, 2, 3), n_o),
       fill.colors = unname(cols),
       xlim        = c(0.9, 1.5),
       ci.sep      = ";",
       ci.edge     = FALSE,
       left.bar    = FALSE,
       right.bar   = FALSE,
       horiz.bar   = FALSE,
       left.align  = TRUE,
       sort        = FALSE,
       margin.left = 16,
       margin.right = 22,
       legend      = FALSE,
       font        = "Arial",
       save        = TRUE,
       save.path   = out,
       save.name   = "forest_pcs_extended_adjustment",
       save.type   = "png",
       save.width  = 11,
       save.height = 5)


## -----------------------------------------------------------------------
## 15g  Export of the analysis samples
## -----------------------------------------------------------------------

##      (shared with the statisticians for the additional revision analyses)

post_export <- merge(get_scores("post"),
                     covid_data[c("ID", "Sequelae.nach.3.Monaten", "Sequelae.nach.6.Monaten")],
                     by.x = "sid", by.y = "ID", all.x = TRUE)
post_export$symptomatic <- ifelse(post_export$cov_severity_a == 1, 0,
                                  ifelse(post_export$cov_severity_a %in% 2:4, 1, NA))
post_export$pcs_3m <- ifelse(post_export$Sequelae.nach.3.Monaten == "yes", 1, 0)
post_export$pcs_6m <- ifelse(post_export$Sequelae.nach.6.Monaten == "yes", 1, 0)
post_export <- post_export[, c("sid", "cov_prot_score_ensemble", "age", "sex",
                               "cov_severity_a", "symptomatic", "pcs_3m", "pcs_6m")]

pre_export <- get_scores("pre")
pre_export$symptomatic <- ifelse(pre_export$cov_severity_a == 1, 0,
                                 ifelse(pre_export$cov_severity_a %in% 2:4, 1, NA))
pre_export <- pre_export[, c("sid", "cov_prot_score_ensemble", "age", "sex",
                             "cov_severity_a", "symptomatic")]

write.xlsx(list(post_infection = post_export, pre_infection = pre_export),
           file.path(out, "analysis_samples.xlsx"), rowNames = FALSE)


################################################################################
## End of script
################################################################################
