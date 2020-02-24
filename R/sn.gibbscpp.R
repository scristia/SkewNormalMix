skewNormalCpp <- function(r, K, nsim, burnin=500) {

    ##mu <- rep(mean(xx), K)
    xx <- r
    alpha0 <- rep(0, K) ## skewness parameter
    ##alpha0 <- c(-3, 0)
    omega0 <- rep(mad(xx), K) ## scale parameter
    omega20 <- omega0^2

    pars <- kmeans(xx, centers=K, nstart=15)
    mu <- sort(pars$centers)
    S <- rep(NA, length(xx))
    for(i in 1:K) S[pars$cluster == order(pars$center)[i]] <- i
    S <- S-1L
    nn <- pars$size[order(pars$centers)]

    eta0 <- rep(1/K, K) ## intitial mixing params
    mat <- .Call("skewnormal_mix", r, K=K, S=S, centers=mu, alpha=alpha0,
                 omega2=omega20, eta=eta0, nsim)

    return(mat)
}

dsnmix <-
    function(r, mixture, K, log=FALSE) {
        if(K == 1){
            mix.mu <- mean(mixture[["mu"]])
            mix.alpha <- mean(mixture[["alpha"]])
            mix.omega <- mean(mixture[["omega"]])
            p <- mean(mixture[["P"]])
        }
        else{
            mix.mu <- colMeans(mixture[["mu"]])
            mix.alpha <- colMeans(mixture[["alpha"]])
            mix.omega <- colMeans(mixture[["omega"]])
            p <- colMeans(mixture[["P"]])
        }
        ## Find likelihood for each component
        lik.comp <- function(r, comp) {
            p[comp]*dsn(r, mix.mu[comp], mix.omega[comp], mix.alpha[comp])
        }
        ## Create likelihood array with likelihoods of from each component
        liks <- sapply(1:K, lik.comp, r=r)

        d <- rowSums(liks, na.rm=TRUE)
        if(log) {
            d <- log(d)
        }
        return(d)
    }



loglik.snmix <-
    function(r, mixture, K, burnin=1) {
        loglike <- dsnmix(r, mixture, K, log=TRUE)
        return(sum(loglike))
    }
