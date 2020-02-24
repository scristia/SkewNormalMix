## Skew normal mixture model
library(devtools)
load_all()
library(sn)
library(mvtnorm)
library(msm)
library(MASS)
library(gtools)
library(truncnorm)

simulated data
omega <- c(4, 1)
omega2 <- omega^2
alpha <- c(-3, 0)
mu <- c(0, 4)

xx <- c(rsn(5000, mu[1], omega[1], alpha[1]), rsn(8000, mu[2], omega[2], alpha[2]))
xx <- xx[sample.int(8000)]
pdf("~/fig/tmp.pdf")
par(bg="white")
plot(density(xx), type="l")
dev.off()
n <- length(xx)

##transformations
delta <- alpha/sqrt(1+alpha^2)
Ey <- mu+omega2*delta*sqrt(2/3.1415)
psi <- omega*delta
sigma2 <- omega2*(1-delta^2)

K = 2
nsim=10000
burnin <- 1:500

set.seed(4321)
res = skewnormal.gibbs(xx, K=K, nsim=nsim)
mus <- colMeans(res$MU[-burnin, ])
omegas <- colMeans(res$OMEGA[-burnin, ])
alphas <- colMeans(res$ALPHA[-burnin, ])
etas <- colMeans(res$ETA[-burnin, ])

set.seed(4321)
res <- skewNormalCpp(r=xx, K=K, nsim=nsim)
mus <- colMeans(res$MU[-burnin, ])
omegas <- colMeans(res$OMEGA[-burnin, ])
alphas <- colMeans(res$ALPHA[-burnin, ])
etas <- colMeans(res$ETA[-burnin, ])

snmixture <- list("mu"=res$MU[-burnin,], "omega"=res$OMEGA[-burnin,],
                  "alpha"=res$ALPHA[-burnin,], "P"=res$ETA[-burnin,])
loglik <- loglik.snmix(xx, snmixture, K)
bic <- -2*loglik + (4*K-1)*log(length(xx))


## PLOTTING
lim <- range(xx, na.rm=TRUE)
pdf("~/fig/tmp.pdf")
par(las=1, mfrow=c(1,1), mar=c(4, 4, 4, 4), bg="white")
hist(xx, breaks = 500, col='lightgray', border='gray', freq=FALSE, main="",
     xlim=lim)
y2 <- seq(min(xx), max(xx), len=5000)
#post.dens <- pis[1]*dsn(y2, mus[1], omegas[1], alphas[1] ) + pis[2]*dsn(y2, mus[2], omegas[2], alphas[2])
#lines(y2, post.dens,lwd=2)
#mx <- max(post.dens)
for(k in 1:K) lines(y2, etas[k]*dsn(y2, mus[k], omegas[k], alphas[k] ), col="gray40", lty=2, lwd=2)
lines(y2, rowSums(sapply(1:K, function(x) etas[x]*dsn(y2, mus[x], omegas[x], alphas[x]))), col="skyblue3", lwd=2)
dev.off()

## check vs truth
mu; mus
omega; omegas
alpha; alphas
c(8/13, 5/13); etas

library(rbenchmark)
benchmark(skewNormalCpp(r=xx, K=2, nsim=200),
          skewnormal.gibbs(xx, K=2, nsim=200, thin=1),
          replications=3)
## over 10x faster

### S4 objects with Rcpp
library(inline)
src <- '
S4 foo(x) ; foo.slot(".Data") = "bar" ; foo.slot("x")=100; return(foo);
'
fun <- cxxfunction(signature(x="any"), src,
                   plugin="Rcpp")
setClass( "S4ex", contains = "character",
         representation( x = "numeric" ) )
x <- new( "S4ex", "bla", x = 10 )
fun(x)
str(fun(x))
