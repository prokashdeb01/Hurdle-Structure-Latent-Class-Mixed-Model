# Call Library
library(lcmm)
library(dplyr)
library(tidyr)
library(tidyverse)
library(haven)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(tidyr)
library(zoo)
library(MASS)
library(plm)
library(pglm)
library(pscl)
library(geepack)
library(statmod)
library(maxLik)
library(parallel)
library(cowplot)


# Import Data
rm(list = ls())
pbm <- read_csv("your data path.csv")


############################### First Hurdle ###############################################

pdata <- pdata.frame(pbm, index = c("new_household_code", "mon"))
pprobit <- pglm(snp ~ Household_Income+Household_Size+
                  male_under_45+female_under_45+male_less_highschool+female_less_highschool+
                  children_yes+african_american+other+hispanic_yes,
                data = pdata,
                family = binomial('probit'),
                model = "pooling", method = "bfgs", print.level = 3, R = 5)

########################################## Second Hurdle ####################################
# Subset the participants data
pbm_1.1 <- subset(pbm, snp==0) %>% 
  dplyr::select(new_household_code, mon,projection_factor, pound, price_hat, price_final_pound, 
                resid_stage1, Household_Income,Household_Size,male_over_45, male_under_45, 
                female_over_45, female_under_45, children_yes, children_no, married, 
                other_than_married, male_less_highschool, male_morethan_highschool, 
                female_less_highschool, female_morethan_highschool, caucasian, african_american, 
                other, hispanic_yes, hispanic_no)

set.seed(123)

pbm_1.1 <- as.data.frame(pbm_1.1)
pbm_1.1$mon <- pbm_1.1$mon/10


m.1 <- hlme(fixed = pound~mon+I(mon^2)+price_final_pound+resid_stage1+Household_Income+Household_Size,
            subject = "new_household_code",
            ng=1,
            maxiter = 1000,
            data = pbm_1.1,
            convB = 0.0001,
            convL = 0.0001,
            convG = 0.0001,
            # partialH = TRUE,
            verbose = TRUE,
            nproc = detectCores()-3)

# Determine the number of cores to use
n_cores <- detectCores() - 3
# Create a cluster for parallel processing
cl <- makeCluster(n_cores)

m.4e <- gridsearch(
  hlme(
    fixed = pound ~ mon+I(mon^2)+price_final_pound+resid_stage1+Household_Income+Household_Size,
    mixture = ~ mon+I(mon^2)+price_final_pound+resid_stage1+Household_Income+Household_Size,
    classmb = ~ male_under_45 + female_under_45 +
      male_less_highschool + female_less_highschool +
      children_yes + african_american + other + hispanic_yes,
    subject = "new_household_code",
    ng = 4,
    nwg = FALSE,
    partialH = TRUE,
    maxiter = 10000,
    data = pbm_1.1,
    convB = 0.001,
    convL = 0.001,
    convG = 0.001,
    verbose = TRUE,
    # nproc = detectCores()-3
  ),
  rep = 100,         # Number of random starts
  maxiter = 50,    # Max iterations for each start
  minit = m.1,      # Initial model with ng = 1
  cl = cl           # Cluster for parallel processing
)


################################ Combined log-likelihood function ###############################

# Initial guesses setup
initial_guess_probit <- coef(pprobit)
initial_guess_probit

initial_guess_lcmm <- estimates(m.4e, cholesky = T)

initial_guess_combined <- c(initial_guess_probit, initial_guess_lcmm)
initial_guess_combined


# Log-likelihood for the probit model
log_likelihood_probit <- function(params_probit) {
  X <- model.matrix(~ price_final_pound + resid_stage1 + Household_Income + Household_Size +
                      male_under_45 + female_under_45 +
                      children_yes + male_less_highschool + female_less_highschool +
                      african_american + other + hispanic_yes,
                    data = pdata)
  
  eta <- X %*% params_probit
  p <- pnorm(eta)  # Probit link function
  y <- pdata$snp
  log_likelihood <- sum(y * log(p) + (1 - y) * log(1 - p))
  
  return(log_likelihood)
}


# Log-likelihood for the LCMM model
log_likelihood_lcmm <- function(params_lcmm) {
  # Ensure params_lcmm is numeric
  # params_lcmm <- as.numeric(params_lcmm)
  pbm_1.1 <- pbm_1.1 %>% arrange(new_household_code, mon)
  # Y0: Numeric vector of outcomes
  Y0 <- as.numeric(pbm_1.1$pound)
  
  # X0: Numeric matrix of covariates
  X0 <- model.matrix(~mon + I(mon^2) + price_final_pound + resid_stage1 + Household_Income + Household_Size +
                       male_under_45 + female_under_45 + children_yes + male_less_highschool + 
                       female_less_highschool + african_american + other + hispanic_yes,
                     data = pbm_1.1)
  
  # Other necessary parameters for loglikhlme
  nmes0 <- as.vector(table(pbm_1.1$new_household_code)) #number of mesures for each subject (length ns0 or dom ns0*ny0)
  ns0 <- length(nmes0) # Number of subjects (24 row per subject)
  ng0 <- 4  # Number of latent classes
  nom.X0 <- colnames(X0)
  nvar.exp <- length(nom.X0)
  nv0 <- nvar.exp #number of covariates
  prior2 <- as.integer(rep(0,ns0))
  prior0 <- prior2 #the prior latent class (length ns0)
  pprior0 <- matrix(1, ns0, ng0) #the prior probabilty of each latent class; 1081 is the number of unique household
  idprob0 <- m.4e[["idprob0"]] #indicator of presence in the class membership submodel (length nv0)
  idea0 <- m.4e[["idea0"]] #indicator of presence in the random part of the longitudinal submodel (length nv0)
  idg0 <- m.4e[["idg0"]] #indicator of presence in the fixed part of the longitudinal submodel (length nv0)
  z.X0 <- strsplit(colnames(X0),split=":",fixed=TRUE)
  z.X0 <- lapply(z.X0,sort)
  idcor0 <- rep(0, length(z.X0)) #indicator of presence in the correlation part of the longitudinal submodel (length nv0)
  nobs0 <- length(Y0)     # Number of observations
  nea0 <- sum(idea0==1)  # Number of random effects
  idiag0 <- as.integer(0) #indicator of diagonal variance matrix of the random effects
  nwg0 <- as.integer(0) #number of parameters for proportional random effects over latent classes
  ncor0 <- 0 #number of parameters for the correlation
  npm0 <- length(m.4e[["best"]])  # Total number of parameters
  fix0 <- rep(0,npm0) #indicator of non estimated parameter (length npm0+nfix0)
  nfix0 <- sum(fix0) # number of non estimated parameters
  bfix0 <- 0 #vector of non estimated parameters
  
  
  # Call the Fortran function to calculate log-likelihood
  log_likelihood <- loglikhlme(b = params_lcmm, Y0 = Y0, X0 = X0, ns0 = ns0, ng0 = ng0, 
                               nobs0 = nobs0, nea0 = nea0, npm0 = npm0, 
                               prior0=prior0, pprior0=pprior0, 
                               idprob0=idprob0, idea0=idea0, idg0=idg0, idcor0=idcor0, nv0=nv0,nmes0=nmes0, 
                               idiag0=idiag0, nwg0=nwg0, ncor0=ncor0, fix0=fix0,
                               nfix0=nfix0, bfix0=bfix0)
  
  return(log_likelihood)
}


# Combined log-likelihood function
combined_log_likelihood <- function(params) {
  m <- length(initial_guess_probit)
  
  params_probit <- params[1:m]
  params_lcmm <- params[(m+1):length(params)]
  
  log_likelihood_probit_val <- log_likelihood_probit(params_probit)
  log_likelihood_lcmm_val <- log_likelihood_lcmm(params_lcmm)
  
  return((log_likelihood_probit_val + log_likelihood_lcmm_val))
}


# Fixed Parameters and MaxLik
library(numDeriv)
library(future.apply)

# Making fixed parameter for 4 class model
print(names(initial_guess_combined))
fixed_params <- rep(FALSE, 69)
fixed_positions <- c(1,3, 14:16, 41:52, 57:60, 69)
fixed_params[fixed_positions] <- TRUE
fixed_params

# Calculate numerical gradient
plan(multisession, workers = parallel::detectCores() - 3)

combined_gradient_function <- function(params) {
  grad(func = combined_log_likelihood, x = params, method = "simple")  # faster than "complex"
}

parameters_DHLCMM <- maxLik(logLik = combined_log_likelihood, 
                            start = initial_guess_combined, 
                            method = "BFGSR",
                            grad = combined_gradient_function,
                            control = list(printLevel = 1, tol = 1e-4, gradtol=1e-4),
                            fixed = fixed_params)

summary(parameters_DHLCMM)













