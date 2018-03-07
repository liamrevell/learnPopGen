\name{drift.selection}
\alias{drift.selection}
\title{Genetic drift & selection simulation}
\usage{
	drift.selection(p0=0.5,Ne=100,w=c(1,1,1),ngen=400,nrep=10,colors=NULL,...)
}
\arguments{
	\item{p0}{starting frequency for the A allele.}
	\item{Ne}{effective population size.}
	\item{w}{fitnesses of the three genotypes: \emph{AA}, \emph{Aa}, and \emph{aa}.}
	\item{nrep}{number of replicate simulations.}
	\item{ngen}{total time, in number of generations, for the simulation.}
	\item{colors}{colors to use for plotting.}
	\item{...}{optional arguments.}
}
\description{
	Simulates drift and natural selection.
}
\value{
	The function creates a plot and returns a list containing the allele frequency through time for each simulation.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
 	\code{\link{genetic.drift}}, \code{\link{selection}}
}
\examples{
	drift.selection()
	drift.selection(p0=0.01,Ne=100,w=c(1,0.9,0.8),ngen=200,nrep=5)
}
\keyword{population genetics}
\keyword{drift}
