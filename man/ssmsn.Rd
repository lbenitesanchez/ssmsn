\name{ssmsn}
\alias{ssmsn}
\alias{dssmsn}
\alias{rssmsn}

\title{Scale-Shape Mixtures of Skew-Normal Distributions}

\description{It provides the density and random number generator.
}
\usage{
dssmsn(x, mu= NULL,sigma2= NULL,lambda= NULL,nu= NULL,family="skew.t.t")
rssmsn(n,mu= NULL,sigma2= NULL,lambda= NULL,nu= NULL,family="skew.t.t")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{
vector of observations.
}
\item{n}{
numbers of observations.
}
  \item{mu}{
location parameter.
}
  \item{sigma2}{
scale parameter.
}
  \item{lambda}{
skewness parameter.
}
  \item{nu}{
degree freedom
}

  \item{family}{
 distribution family to be used in fitting ("skew.t.t", "skew.generalized.laplace.normal, "skew.slash.normal")
}

}

\details{
As discussed in Jamalizadeh and Lin (2016) the scale-shape mixture of skew-normal (SSMSN) distribution admits the following conditioning-type stochasctic representation

\deqn{Y=\mu + \sigma \tau_1^{-1/2}[Z_1 | (Z_2 < \lambda f^{-1/2} Z_1)],}

where \eqn{f = \tau_1}/\eqn{\tau_2} and \eqn{(Z_1,Z_2)} and (\eqn{\tau_1,\tau_2}) are independent. Alternatively the SSMSN distribution can be generated via the convolution-type stochastic representation, given by

\deqn{Y=\mu + \sigma \left(\frac{\tau_1^{-1/2} f^{1/2}}{\sqrt{f + \lambda^2}}Z_2 + \frac{\lambda \tau_1^{-1/2}}{\sqrt{f + \lambda^2}}|Z_1|\right).}

}


\value{
\code{dssmsn} gives the density, \code{rssmsn} generates a random sample.

The length of the result is determined by n for \code{rssmsn}, and is the maximum of the lengths of the numerical arguments for the other functions \code{dssmsn}.
}

\references{
Jamalizadeh, Ahad and Lin, Tsung-I (2016). A general class of scale-shape mixtures of skew-normal distributions: properties and estimation. Computational Statistics, 1-24.
}

\author{Rocio Maehara \email{rmaeharaa@gmail.com} and Luis Benites \email{lbenitesanchez@gmail.com}
}




\examples{

rSTT  <- rssmsn(n=1000,mu=-4,sigma2=1,lambda=1,nu=c(3,4),"skew.t.t");hist(rSTT)
rSGLN <- rssmsn(n=1000,mu=-4,sigma2=1,lambda=1,nu=3,"skew.generalized.laplace.normal");hist(rSGLN)
rSSN  <- rssmsn(n=1000,mu=-4,sigma2=1,lambda=1,nu=3,"skew.slash.normal");hist(rSSN)

dSTT  <- dssmsn(0.5,mu=-4,sigma2=1,lambda=1,nu=c(3,4),"skew.t.t")
dSGLN <- dssmsn(0.5,mu=-4,sigma2=1,lambda=1,nu=3,"skew.generalized.laplace.normal")
dSSN  <- dssmsn(0.5,mu=-4,sigma2=1,lambda=1,nu=3,"skew.slash.normal")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Shape-scale mixture of skew-normal distributions}
\keyword{ssmsn}
\keyword{models}
