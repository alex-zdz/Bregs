nullreg<-function(X,Y,
                  nburn,nthin,n_save,
                  prior.mean.beta=NULL, prior.var.beta=NULL, prior.sig2=NULL,
                  sig2.samp= NULL, beta.samp = NULL,
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
  mod.glm   <- glm(Y ~ X)
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















