\name{clt}
\alias{clt}
\title{Illustrates the concept of the Central Limit Theorem}
\usage{
	clt(nvar=1, nobs=1000, df=c("normal", "uniform", "exponential", "binomial"),
		theta=if(df%in%c("normal", "uniform", "exponential")) 1 else 0.5,
		breaks="Sturges", show=c("sum", "mean"))
}
\arguments{
	\item{nvar}{number of random variables to sum (1 or more).}
	\item{nobs}{total number of observations (per random variable).}
	\item{df}{distribution functions of individual random variables to sum.}
	\item{theta}{parameter of the distribution functions: variance in the case of \code{df = "normal"}, maximum value in the case of \code{"uniform"} (the minimum value will be assumed to be 0), and rate in the case of \code{"exponential"}.}
	\item{breaks}{breaks (see \code{\link{hist}}).}
	\item{show}{whether to show the row-wise \emph{sum} of the independent random variables (\code{show="sum"}), or their mean (\code{show="mean"}).}
}
\description{
	The \emph{Central Limit Theorem} tells us that when independent random variables are added together, the distribution of their sum tends towards a normal distribution, regardless of the shape of their individual distributions. This function attempts to illustrate this concept by allowing the user to visualize the sum of an arbitrary number of different independent random variables with different underlying distributions.
}
\details{
	The central limit theorem (CLT) establishes that when independent random variables are added together their (normalized) sum will tend towards a normal distribution, regardless of the distribution of the original random variables. That is to say if we were to generate a set of \code{nvar} (say) independent, uniform, random variables, normalize each one to have the same variance, and then sum or average the variables by observation, this sum or average will tend towards a normal distribution as \code{nvar} is increased.
}
\value{
	Creates a plot showing the observation-wise distribution of the sum or average of the independent random variables.
	
	The distribution of the observation-wise sum or average and the underlying data are also returned invisibly to the user in the form of an object of class \code{"clt"}.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{phenotype.freq}}
}
\examples{
	clt(nvar=1,df="exponential")
	clt(nvar=10,df="exponential")
	clt(nvar=40,df="exponential")
}
\keyword{statistics}
