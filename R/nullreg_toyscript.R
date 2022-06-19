
set.seed(1234)
n     = 100
alpha = 0
beta  = 0.9
#beta  = 1.2
sig   = 1
z     = rep(0,n)
z[1]  = alpha/(1-beta)
e     = rnorm(n,0,sig)
for (t in 2:n)
  z[t] = alpha+beta*z[t-1]+e[t]

par(mfrow=c(2,2))
ts.plot(z,main="Simulated AR(1) data")

y = z[2:n]
x <- z[1:(n-1)]
#X = cbind(1,x)
X <- matrix(x,length(x),1)

p <- ncol(X)
Y <- y

nullreg(
        X,Y,
        nburn = 0 ,nthin = 1,n_save = 1e3,
        prior.mean.beta=NULL, prior.var.beta=NULL, prior.sig2=NULL,
        sig2.samp= NULL, beta.samp = NULL,
        verbose=TRUE
)

X;Y;
nburn = 0 ;nthin = 1;n_save = 1e3;
prior.mean.beta=NULL; prior.var.beta=NULL; prior.sig2=NULL;
sig2.samp= NULL; beta.samp = NULL;
verbose=TRUE



mu <- 5
sig2 <- 2
set.seed(28)
y <- rnorm(110,mu, sqrt(sig2))
X <- matrix(1,110,1)
Y <- y
