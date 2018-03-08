\name{msd}
\alias{msd}
\title{Migration, drift, and selection}
\usage{
	msd(p0=c(0.5,0.5),Ne=c(100,100),w=list(c(1,1,1),c(1,1,1)),
		m=c(0.01,0.01),ngen=400,colors=c("red","blue"),...)
}
\arguments{
	\item{p0}{starting frequency for the A allele in each of two populations.}
	\item{Ne}{effective population size for each of two populations.}
	\item{w}{fitnesses of the three genotypes (\emph{AA}, \emph{Aa}, and \emph{aa}, in that order) in each of the two populations.}
	\item{m}{rate of migration \emph{from} each of the two populations to the other.}
	\item{ngen}{total time, in number of generations, for the simulation.}
	\item{colors}{colors to use for plotting.}
	\item{...}{optional arguments.}
}
\description{
	Simulates migration, natural selection, and genetic drift. Selection can be in opposite directions in the two populations experiencing gene flow.
}
\value{
	The function creates a plot and returns a list containing the allele frequency through time for each of the two simulated populations.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
 	\code{\link{drift.selection}}
}
\examples{
	msd()
	msd(p0=c(0.25,0.75),w=list(c(1,0.9,0.8),c(0.8,0.9,1)))
	msd(p0=c(0.25,0.75),w=list(c(1,0.9,0.8),c(0.8,0.9,1)),
		m=c(0.1,0.1))
}
\keyword{population genetics}
\keyword{drift}
\keyword{migration}
