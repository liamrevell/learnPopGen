\name{clt}
\alias{clt}
\title{Illustrates the concept of the Central Limit Theorem}
\usage{
	clt(nvar=1, nobs=1000, df=c("normal","uniform","exponential"), theta=1,
		breaks="Sturges")
}
\arguments{
	\item{nvar}{number of random variables to sum (1 or more).}
	\item{nobs}{total number of observations (per random variable).}
	\item{df}{distribution functions of individual random variables to sum.}
	\item{theta}{parameter of the distribution functions: variance in the case of \code{df = "normal"}, maximum value in the case of \code{"uniform"} (the minimum value will be assumed to be 0), and rate in the case of \code{"exponential"}.}
	\item{breaks}{breaks (see \code{\link{hist}}).}
}
\description{
	The Central Limit Theorem tells us that when independent random variables are added together, the distribution of their sum tends towards a normal distribution, regardless of the shape of their individual distributions. This function attempts to illustrate this concept by allowing the user to visualize the sum of an arbitrary number of different independent random variables with different underlying distributions.
}
\value{
	Creates a plot.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\examples{
	clt(nvar=1,df="exponential")
	clt(nvar=10,df="exponential")
	clt(nvar=40,df="exponential")
}
\keyword{statistics}
