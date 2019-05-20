\name{drift.selection}
\alias{drift.selection}
\title{Simulation of genetic drift & natural selection at a biallelic locus}
\usage{
	drift.selection(p0=0.5, Ne=100, w=c(1,1,1), ngen=400, nrep=10, colors=NULL, ...)
}
\arguments{
	\item{p0}{starting frequency for the \emph{A} allele.}
	\item{Ne}{effective population size.}
	\item{w}{fitnesses of the three genotypes: \emph{AA}, \emph{Aa}, and \emph{aa}.}
	\item{nrep}{number of replicate simulations.}
	\item{ngen}{total time, in number of generations, for the simulation.}
	\item{colors}{colors to use for plotting.}
	\item{...}{optional arguments. Presently the only arguments are \code{type} (e.g., \code{"l"}, \code{"s"}) and \code{lwd} in the \code{plot} method.}
}
\description{
	Simulates drift and natural selection at a single biallelic locus within one or various populations.
}
\value{
	The function creates a plot and returns an object of class \code{"drift.selelction"} consisting of list containing the allele frequency through time for each simulation. This object can be printed or plotted using corresponding methods. (See examples.)
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
 	\code{\link{genetic.drift}}, \code{\link{selection}}
}
\examples{
	drift.selection()
	p<-drift.selection(p0=0.01,Ne=100,w=c(1,0.9,0.8),ngen=200,nrep=5)
	print(p)
	plot(p)
}
\keyword{population genetics}
\keyword{drift}
