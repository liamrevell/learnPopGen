\name{clt}
\alias{clt}
\title{Illustrates concept of central limit theorem}
\usage{
	clt(nvar=1, nobs=1000, df=c("normal","uniform","exponential"), theta=1,
		breaks="Sturges")
}
\arguments{
	\item{nvar}{number of random variables to add.}
	\item{nobs}{total number of observations.}
	\item{df}{distribution functions of individual random variables to add.}
	\item{theta}{parameter of the distribution functions: variance in the case of \code{df = "normal"}, maximum value in the case of \code{"uniform"}, and rate in the case of \code{"exponential"}.}
	\item{breaks}{breaks (see \code{\link{hist}}).}
}
\description{
	Illustrates CLT.
}
\value{
	Creates a plot or animation.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\examples{
	clt(df="exponential")
}
\keyword{statistics}
