rssmsn <-
  function(
    n,           # sample size
    mu= NULL,          # location parameter
    sigma2= NULL,      # scale parameter
    lambda= NULL,      # skewness parameter
    nu= NULL,          # degree freedom
    family="skew.t.t"  # "skew.t.t", "skew.generalized.laplace.normal", "skew.slash.normal"
  ) {
    if(family=="skew.t.t")
    {
      tau1 <- rgamma(n,nu[1]/2,nu[1]/2)
      tau2 <- rgamma(n,nu[2]/2,nu[2]/2)
    }

    if(family=="skew.generalized.laplace.normal")
    {
      tau1 <- rinvgamma(n, shape=nu, 0.5)
      tau2 <- 1
    }

    if(family=="skew.slash.normal")
    {
      tau1 <- rbeta(n, nu/2, 1)
      tau2 <- 1
    }

    Z1     <- rnorm(n)
    Z2     <- rnorm(n)
    f      <- tau1/tau2
    y      <- mu + sqrt(sigma2)*((tau1^(-0.5)*f^0.5)/sqrt(f + lambda^2)*Z2 + (lambda*tau1^(-0.5))/sqrt(f + lambda^2)*abs(Z1))

    return(
      y
    )
  }

#yrSTT  <- rssmsn(n=1000,mu=-4,sigma2=1,lambda=1,nu=c(3,4),"skew.t.t");hist(yrSTT)
#yrSGLN <- rssmsn(n=1000,mu=-4,sigma2=1,lambda=1,nu=3,"skew.generalized.laplace.normal");hist(yrSGLN)
#yrSSN  <- rssmsn(n=1000,mu=-4,sigma2=1,lambda=1,nu=3,"skew.slash.normal");hist(yrSSN)


dssmsn <-
  function(
    x,           # evaluation points
    mu= NULL,          # location parameter
    sigma2= NULL,      # scale parameter
    lambda= NULL,      # skewness parameter
    nu= NULL,          # degree freedom
    family="skew.t.t"  # "skew.t.t", "skew.generalized.laplace.normal", "skew.slash.normal"
  ) {
    if(family=="skew.t.t")
    {
      u    <- (x - mu)/sqrt(sigma2)
      eta  <- lambda*u
      dens <- (2/sqrt(sigma2))*dt(u, df = nu[1])*pt(eta,nu[2])
    }

    if(family=="skew.generalized.laplace.normal")
    {
      u    <- (x - mu)/sqrt(sigma2)
      fGL  <- (1/(sqrt(2*pi)*2^(nu-1)*gamma(nu)*sqrt(sigma2)))*abs(u)^(nu-0.5)*besselK(abs(u), 0.5-nu)
      dens <- 2*fGL*pnorm(lambda*u)
    }

    if(family=="skew.slash.normal")
    {
      u    <- (x - mu)/sqrt(sigma2)
      eta  <- lambda*u
      if(u == 0)
      {
        fS <- nu/(sqrt(sigma2)*sqrt(2*pi)*(nu + 1))
      }else{
        fS <- (nu*gamma((nu + 1)/2)*2^(nu/2 - 1))/(sqrt(sigma2)*sqrt(pi)) * pgamma(u^2/2,(nu+1)/2,1)/abs(u)^(nu+1)
      }
      dens <- 2*fS*pnorm(eta)
    }

    return(
      dens
    )
  }


#dSTT  <- dssmsn(0.5,mu=-4,sigma2=1,lambda=1,nu=c(3,4),"skew.t.t"); dSTT
#dSGLN <- dssmsn(0.5,mu=-4,sigma2=1,lambda=1,nu=3,"skew.generalized.laplace.normal"); dSGLN
#dSSN  <- dssmsn(0.5,mu=-4,sigma2=1,lambda=1,nu=3,"skew.slash.normal"); dSSN
