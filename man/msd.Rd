\name{msd}
\alias{msd}
\title{Migration, drift, and selection}
\usage{
	msd(p0=c(0.5,0.5), Ne=c(100,100), w=list(c(1,1,1),c(1,1,1)), 
		m=c(0.01,0.01), ngen=400, colors=c("red","blue"), ...)
}
\arguments{
	\item{p0}{starting frequency for the \emph{A} allele in each of two populations.}
	\item{Ne}{effective population size for each of two populations.}
	\item{w}{fitnesses of the three genotypes (\emph{AA}, \emph{Aa}, and \emph{aa}, in that order) in each of the two populations. \code{w} should take the form of a list of two vectors.}
	\item{m}{rates of migration \emph{from} the first population to the second, and from the second population to the first, in that order. This value is best interpreted as the \emph{probability} that an individual born in population 1 will migrate to population 2 before reproduction, and \emph{vice versa}.}
	\item{ngen}{total time, in number of generations, for the simulation.}
	\item{colors}{colors to use for plotting.}
	\item{...}{optional arguments. Presently, the only optional argument for \code{msd} is \code{show.legend} (which defaults to \code{TRUE}). The \code{plot} method adds the additional optional arguments of \code{colors} (a vector of colors for the two simulated populations), \code{lwd}, and \code{type} (e.g., \code{"l"} or \code{"s"}).}
}
\description{
	Simulates migration, natural selection, and genetic drift. Selection can be in opposite directions in the two populations experiencing gene flow.
}
\value{
	The function creates a plot and invisibly returns a list containing the allele frequency through time for each of the two simulated populations.
	
	The returned object is of class \code{"msd"} and can be printed or re-plotted using corresponding \code{print} or \code{plot} methods. (See examples.)
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
 	\code{\link{drift.selection}}
}
\examples{
	msd()
	msd(p0=c(0.25,0.75),w=list(c(1,0.9,0.8),c(0.8,0.9,1)))
	object<-msd(p0=c(0.75,0.25),w=list(c(1,0.9,0.8),
		c(0.8,0.9,1)),m=c(0.1,0.1))
	print(object)
	plot(object,colors=c("black","grey"),lwd=4)
}
\keyword{population genetics}
\keyword{drift}
\keyword{migration}
