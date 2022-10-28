## Main div rate calculation script

## load the packages:

library(TESS)

CDSTree <- read.tree("../TimeTrees/GTRG/3TaxCalTrees/CDS_CLCK_medianTT.tree")
intronTree <- read.tree("../TimeTrees/GTRG/3TaxCalTrees/intron_CLCK_medianTT.tree")
UTRTree <- read.tree("../TimeTrees/GTRG/3TaxCalTrees/UTR_CLCK_medianTT.tree")
lnc_RNATree <- read.tree("../TimeTrees/GTRG/3TaxCalTrees/lnc_RNA_CLCK_medianTT.tree")
otherTree <- read.tree("../TimeTrees/GTRG/3TaxCalTrees/other_CLCK_medianTT.tree")
ALL50KTree <- read.tree("../TimeTrees/GTRG/3TaxCalTrees/ALL50K_CLCK_medianTT.tree")
unannotatedTree <- read.tree("../TimeTrees/GTRG/3TaxCalTrees/unannotated_CLCK_medianTT.tree")
pseudogeneTree <- read.tree("../TimeTrees/GTRG/3TaxCalTrees/pseudogene_CLCK_medianTT.tree")




CDS_times <- as.numeric( branching.times(CDSTree) )


#Setup the priors

CDS_prior_delta <- function(x) { dexp(x,rate=10.0,log=TRUE) }
CDS_prior_tau <- function(x) { dexp(x,rate=10.0,log=TRUE) }
CDS_priorsConstBD <- c("diversification"=CDS_prior_delta,
                   "turnover"=CDS_prior_tau)


##Likelihood function

CDS_likelihoodConstBD <- function(params) {
  
  speciation <- params[1] + params[2]
  extinction <- params[2]
  
  lnl <- tess.likelihood(CDS_times,
                         lambda = speciation,
                         mu = extinction,
                         samplingProbability = 0.68,
                         log = TRUE)
  
  return (lnl)
  
}
                         

## MCMC for CDS

CDS_samplesConstBD <- tess.mcmc(likelihoodFunction = CDS_likelihoodConstBD,
                            priors = CDS_priorsConstBD,
                            parameters = runif(2,0,1),
                            logTransforms = c(TRUE,TRUE),
                            delta = c(1,1),
                            iterations = 10000,
                            burnin = 1000,
                            thinning = 10,
                            adaptive = TRUE,
                            verbose = TRUE)


CDS_summary <- summary(CDS_samplesConstBD)

plot(CDS_samplesConstBD)


## varying rates

CDS_prior_delta <- function(x) { dexp(x,rate=0.1,log=TRUE) }
CDS_prior_lambda <- function(x) { dexp(x,rate=10.0,log=TRUE) }
CDS_prior_alpha <- function(x) { dexp(x,rate=0.1,log=TRUE) }
CDS_priorsDecrBD <- c("turnover"=CDS_prior_delta,
                  "initial speciation"=CDS_prior_lambda,
                  "speciation decay"=CDS_prior_alpha)
            
#likelihood function

CDS_likelihoodDecrBD <- function(params) {
  speciation <- function(t) params[1] + params[2] * exp(-params[3]*t)
  extinction <- function(t) params[1]
  lnl <- tess.likelihood(CDS_times,
                         lambda = speciation,
                         mu = extinction,
                         samplingProbability = 0.68,
                         log = TRUE)
  
  return (lnl)
  
}

CDS_samplesDecrBD <- tess.mcmc(likelihoodFunction = CDS_likelihoodDecrBD,
                           priors = CDS_priorsDecrBD,
                           parameters = runif(3,0,1),
                           logTransforms = c(TRUE,TRUE,TRUE),
                           delta = c(1,1,1),
                           iterations = 10000,
                           burnin = 1000,
                           thinning = 10,
                           adaptive = TRUE,
                           verbose = TRUE)



plot(CDS_samplesDecrBD)







#model evalution

marginalLikelihoodConstBD <- tess.steppingStoneSampling(
  likelihoodFunction = likelihoodConstBD,
  priors = priorsConstBD,
  parameters = runif(2,0,1),
  logTransforms = c(TRUE,TRUE),
  iterations = 1000,
  burnin = 100,
  K = 50)


marginalLikelihoodDecrBD <- tess.steppingStoneSampling(
  likelihoodFunction = likelihoodDecrBD,
  priors = priorsDecrBD,
  parameters = runif(3,0,1),
  logTransforms = c(TRUE,TRUE,TRUE),
  iterations = 1000,
  burnin = 100,
  K = 50)

# compare the bayesion factors

candidateModels <- c("ConstBD"=marginalLikelihoodConstBD,
                     "DecrBD"=marginalLikelihoodDecrBD)


# Make all possible combinations of the models.
marginalLikelihoodGrid <- expand.grid(M0=names(candidateModels),
                                      M1=names(candidateModels))

marginalLikelihoodGrid$BF <- 2 * (candidateModels[marginalLikelihoodGrid$M0] -
                                    candidateModels[marginalLikelihoodGrid$M1])


marginalLikelihoodGrid <- marginalLikelihoodGrid[order(marginalLikelihoodGrid$BF,
                                                       decreasing=TRUE),]


marginalLikelihoodGrid




