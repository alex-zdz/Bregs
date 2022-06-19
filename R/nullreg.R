#' Title
#'
#' @param X Design matrix
#' @param Y Respones variable
#' @param nburn   MCMC Burnin
#' @param nthin   MCMC Thinning
#' @param n_save  MCMC Samples saved
#' @param prior.mean.beta Prior mean of the regression coefficient
#' @param prior.var.beta  Prior variance of the regression coefficient
#' @param prior.sig2      Prior for the variance term
#' @param sig2.samp       Starting value for sigma^2
#' @param beta.samp       Starting value for beta
#' @param verbose         Print progress
#'
#' @return                Posterior beta and sigma samples
#' @export
#'
#'
#'
#'
#' @examples # AR-1 model:
#'set.seed(1234)
#' n     = 100
#' alpha = 0
#' beta  = 0.9
#' #beta  = 1.2
#' sig   = 1
#' z     = rep(0,n)
#' z[1]  = alpha/(1-beta)
#' e     = rnorm(n,0,sig)
#' for (t in 2:n)
#'   z[t] = alpha+beta*z[t-1]+e[t]
#'
#' par(mfrow=c(2,2))
#' ts.plot(z,main="Simulated AR(1) data")
#'
#' Y = z[2:n]
#' x <- z[1:(n-1)]
#' X <- matrix(x,length(x),1)
#' nullreg(
#' X,Y,
#' nburn = 0 ,nthin = 1,n_save = 1e3,
#' prior.mean.beta=NULL, prior.var.beta=NULL, prior.sig2=NULL,
#' sig2.samp= NULL, beta.samp = NULL,
#' verbose=TRUE
#' )
#'
nullreg<-function(X,Y,
                  nburn,nthin,n_save,
                  prior.mean.beta=NULL, prior.var.beta=NULL, prior.sig2=NULL,
                  sig2.samp= NULL, beta.samp = NULL,intercept = FALSE,
                  verbose=TRUE){



  ## extract dimensions from design matrix
  #need to include check if X already includes an intercept

  if(length(dim(X)) == 3){
    stop("Not yet implemented for panel data", call.=FALSE)
    I<-dim(X)[1]
    Ttot<-dim(X)[2]
    p<-dim(X)[3]
  }else{
    I<-dim(X)[1]
    p<-dim(X)[2]
  }


  ## Number of MCMC iterations
  n_sample <-  nburn+ nthin*n_save

  #### Priors
  if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
  if(is.null(prior.var.beta))  prior.var.beta  <- diag(rep(1, p))
  if(is.null(prior.sig2)) prior.sig2 <- c(1, 0.01)

  data.precision.beta <- crossprod(X)
  if(length(prior.var.beta)==1){
    prior.precision.beta <- 1 / prior.var.beta
  }else{
    prior.precision.beta <- diag(1/prior.var.beta)
  }



  ###########################
  ## Initial parameter values
  ###########################
  if(intercept){
  mod.glm   <- glm(Y ~ -1 + X)
  }else{
  mod.glm   <- glm(Y ~  X)
  }
  beta0     <- mod.glm$coefficients[1]
  beta.mean <- mod.glm$coefficients[-1]
  beta.sd <- sqrt(diag(summary(mod.glm)$cov.scaled))[-1]
  beta.samp <- rnorm(n=length(beta.mean), mean=beta.mean, sd=beta.sd)

  res.temp <- Y - beta0 - X %*% beta.samp
  res.sd <- sd(res.temp, na.rm=TRUE)
  sig2.samp <- res.sd/2


  ###########################
  ## Results storage
  ###########################

  beta_out <- vector("list", n_save)
  sig2_out <- rep(NA, n_save)

  ###########################
  ## pre-calculated values:
  ###########################
  pre1 <- crossprod(X) + prior.precision.beta
  pre2 <- crossprod(X, Y) #+ prior.precision.beta %*% prior.mean.beta
  pre3 <- prior.sig2[1] + I/2 + p/2

  ###########################
  ## Run the model
  ###########################
  ## Start timer
  if(verbose) {
    cat("Generating", n_save, "post burnin and thinned (if requested) samples.\n", sep = " ")
    progressBar <- txtProgressBar(style = 3)
    percentage.points<-round((1:100/100)*n_sample)
  } else {
    percentage.points<-round((1:100/100)*n_sample)
  }


  ###########################
  ## Start Gibbs
  ###########################

  for (iter in 1:n_sample){

  # Sample beta:
  Sigma <- solve(pre1/sig2.samp)
  mu <- Sigma %*% (pre2/sig2.samp)
  beta.samp <- mvtnorm::rmvnorm(1, mu, Sigma)

  # Sample sigma^2:
  tmp <- prior.sig2[2] + .5*(crossprod(y - X%*%beta.samp) +
                    t(beta.samp-prior.mean.beta) %*% prior.precision.beta %*% (beta.samp-prior.mean.beta))
  sig2.samp <- 1/rgamma(1, pre3, rate=tmp)



  if(iter > nburn && ((iter - nburn) %% nthin) == 0){
    iter_aux <- (iter- nburn)/nthin
    beta_out[[iter_aux]] <- beta.samp
    sig2_out[iter_aux]<-sig2.samp
  }

  }

  return(list(beta_out = beta_out,sig2_out=sig2_out))

}















