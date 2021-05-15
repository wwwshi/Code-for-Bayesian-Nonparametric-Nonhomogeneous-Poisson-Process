rm(list = ls(all = TRUE))
library(Rcpp)
library(spatstat)

sourceCpp("PP_Rcpp_update_z.cpp")
source("function.R")

# load earthquake data
load("data_train.Rdata")
data_train <- as.vector(quakecounts$counts)
n.grid <- length(data_train)

#function for Collapsed Gibbs Sampler
NPPMFM <- function(data, ini.K, hyper, n.MCMC.iter, n.grid){
  # Initialize z
  ini.z <- c(sample(1:ini.K, size = ini.K, replace = FALSE),
             sample(1:ini.K, size = n.grid - ini.K, replace = TRUE))
  
  # Initialization
  cur.sample <- Initialization(hyper, data, ini.K, ini.z)
  
  save.result <- list()
  save.result$z <- array(dim = c(n.MCMC.iter, n.grid))
  save.result$K <- rep(NA, n.MCMC.iter)
  save.result$lambda <- save.result$N.k <- vector("list", n.MCMC.iter)
  
  for(i.iter in 1:n.MCMC.iter)
  {
    # Update z
    cur.sample <- UpdateZ(cur.sample, data, hyper, n.grid)
    
    # Update lambda
    cur.sample$lambda <- UpdateLambda(cur.sample, data, hyper)
    
    save.result$z[i.iter, ] <- cur.sample$z
    save.result$K[i.iter] <- cur.sample$K
    save.result$N.k[[i.iter]] <- cur.sample$N.k
    save.result$lambda[[i.iter]] <- cur.sample$lambda
    
    cat("iteration:", i.iter, "\n", cur.sample$N.k, "\n")
  }
  return(save.result)
}

## Dahl's method to summarize the samples from the MCMC
GetDahl <- function(MFMfit, n.burn.in)
{
  ################################################################
  
  ## Input: MFMfit = the result from NPPMFM ##
  ##        n.burn.in = the number of burn-in iterations ##
  
  ## Output: 
  ##         Dahl.res = estimated output ##
  
  #################################################################
  z <- MFMfit$z[-(1:n.burn.in),]
  n.iter <- dim(z)[1]
  n.grid <- dim(z)[2]
  sum.matrix <- matrix(0, nrow = n.grid, ncol = n.grid)
  for(i in 1:n.iter){
    membership.matrix <- outer(z[i,], z[i,], "==")
    sum.matrix <- sum.matrix + membership.matrix
  }
  membership.avg <- sum.matrix / n.iter
  min.sq.error = n.grid^2
  for(i in 1:n.iter){
    membership.matrix <- outer(z[i,], z[i,], "==")
    sq.error <- sum((membership.matrix - membership.avg)^2)
    if(sq.error < min.sq.error){
      min.sq.error <- sq.error
      Dahl.index <- i
    }
  }
  iter.index <- n.burn.in + Dahl.index
  Dahl.res <- list(z = MFMfit$z[iter.index,], K = MFMfit$K[iter.index], 
                   N.k = MFMfit$N.k[[iter.index]], lambda = MFMfit$lambda[[iter.index]])
  attr(Dahl.res, "iterIndex") <- iter.index
  attr(Dahl.res, "burnin") <- n.burn.in
  return(Dahl.res)
}

# hyper-parameter values setting #############################
hyper <- list()

# prior for lambda
hyper$a <- 1 #shape
hyper$b <- 1 #rate

# DP prior
hyper$g <- 1

# Poisson
hyper$l <- 1

# VN
hyper$VN <- fn_calc_VN(n.grid, hyper$g, hyper$l)

set.seed(1234)
ini.K <- 5
n.MCMC.iter <- 20000
n.burnin <- 5000
#memory.limit(size = 3000000)
res_train <- NPPMFM(data_train, ini.K, hyper, n.MCMC.iter, n.grid)
save(res_train, file = paste0("tt_", n.grid, "_", "results", ".Rdata"))
Dahl_train <- GetDahl(res_train, n.burnin)
save(Dahl_train, file = paste0("tt_", n.grid, "_", "Dahl", ".Rdata"))

# MAE Comparison
load("tt_7292_Dahl.Rdata")
load("tt_7292_results.Rdata")

load("data_test.Rdata")
data_test <- as.vector(quakecounts$counts)

Dahl.lambda.matrix <- Dahl_train$lambda[Dahl_train$z]
lambda.matrix <- sapply(5001:20000, function(i) res_train$lambda[[i]][res_train$z[i, ]])
lambda.mean <- rowMeans(lambda.matrix)

Dahl.MAE <- mean(abs(data_test - Dahl.lambda.matrix/4))
mean.MAE <- mean(abs(data_test - lambda.mean/4))

d <- as.matrix(data_train)
poisson_mix_model <- stepFlexmix(d ~ 1, k = 1:15, nrep = 3, model = FLXMCmvpois(), 
                                 control = list(tolerance = 1e-8, iter = 1000))
bestfit <- getModel(poisson_mix_model, which = "BIC")
pois.mix_intensity <- unname(parameters(bestfit)[clusters(bestfit)])
pois.mix.MAE <- mean(abs(data_test - pois.mix_intensity/4))
all.MAE <- data.frame(mean.MAE, Dahl.MAE, pois.mix.MAE)
colnames(all.MAE) <- c("MFM-NHPP Mean", "MFM-NHPP Dahl", "Poisson Mixture")
save(all.MAE, file = "tt_MAE.Rdata")
