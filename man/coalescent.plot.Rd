\name{coalescent.plot}
\alias{coalescent.plot}
\title{Creates a (usually animated) simulation of gene coalescence within a population}
\usage{
	coalescent.plot(n=10, ngen=20, colors=NULL, ...)
}
\arguments{
	\item{n}{number of haploid individuals or gene copies.}
	\item{ngen}{number of generations.}
	\item{colors}{colors to use for plotting individuals and lines.}
	\item{...}{optional arguments.}
}
\description{
	Coalescence or coalescent theory is a model for genetic drift within a population in which we envision gene copies merging or "coalescing" into ancestors in the past. This function generates a(n) (optionally animated) visualization of this process of coalescence within a population.
}
\value{
	Creates a plot or animation.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{drift.selection}}, \code{\link{genetic.drift}}
}
\examples{
	coalescent.plot()
	coalescent.plot(n=20,ngen=30,col.order="alternating")
}
\keyword{population genetics}
