# Initialize the samples #################################
Initialization <- function(hyper, data, ini.K, ini.z){
  ini.sample <- NULL
  ini.sample$K <- as.integer(ini.K)
  ini.sample$z <- as.integer(ini.z)
  #ini.sample$N.k <- table(ini.z)
  ini.sample$N.k <- rep(NA, ini.sample$K)  #size of cluster
  for(k in 1:ini.sample$K) ini.sample$N.k[k] <- sum(ini.sample$z == k)
  ini.sample$lambda <- rgamma(ini.K, shape = hyper$a, rate = hyper$b)
  
  return(ini.sample)
}

UpdateZ <- function(this.sample, data, hyper, N.grid){
  #In R, I define N.k as a length K vector
  #But in Rcpp, I define N.k as a N.grid vector which is its max dim.  
  #So I fill the rest of it as 0. The same for lambda.
  
  N.k <- c(this.sample$N.k, rep(as.integer(0), (N.grid - this.sample$K)))
  lambda <- c(this.sample$lambda, rep(as.integer(0), (N.grid - this.sample$K)))
  
  fn_z_update(data, hyper$g, hyper$a, hyper$b, hyper$VN, this.sample$z, this.sample$K, N.k, lambda, N.grid)
  
  this.sample$N.k <- N.k[1:this.sample$K]
  this.sample$lambda <- lambda[1:this.sample$K]
  
  return(this.sample)
}

UpdateLambda <- function(this.sample, data, hyper){
  lambda <- sapply(1:this.sample$K, function(k){
    rgamma(1, shape = hyper$a + sum(data[this.sample$z == k]), rate = hyper$b + this.sample$N.k[k])
  })
  
  return(lambda)
}
