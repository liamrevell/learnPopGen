\name{coalescent.plot}
\alias{coalescent.plot}
\title{Computes phenotypic distribution and its change through time due to natural selection on a polygenic trait}
\usage{
	coalescent.plot(n=10,ngen=20,colors=NULL,...)
}
\arguments{
	\item{n}{number of haploid individuals or gene copies.}
	\item{ngen}{number of generations.}
	\item{colors}{colors to use for plotting individuals and lines.}
	\item{...}{optional arguments.}
}
\description{
	Creates a coalescent genealogy plot.
}
\value{
	Creates a plot or animation.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\examples{
	coalescent.plot()
}
\keyword{population genetics}
