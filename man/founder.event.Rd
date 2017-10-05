\name{founder.event}
\alias{founder.event}
\title{Simulation of a founder event}
\usage{
	founder.event(p0=0.5,Ne=1000,Nf=10,ttime=100,etime=50,show="p")
}
\arguments{
	\item{p0}{Starting frequency for the A allele.}
	\item{Ne}{Effective population size at the start of the simulation and after the founding event.}
	\item{Nf}{Size of the founding population.}
	\item{ttime}{Total time of the simulation.}
	\item{etime}{Time for the founding event.}
	\item{show}{Two different options for plotting. \code{"p"} shows the frequency of the A allele through time; \code{"var"} shows the genetic variation in the population, calculated as p*(1-p). The default is \code{show="p"}.}
}
\description{
	This function simulates genetic drift with a founding event at time \code{etime}.
}
\value{
	The function creates one of two different plots, depending on the value of \code{show}.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
  \code{\link{genetic.drift}}
}
\examples{
	founder.event()
	founder.event(show="variation")
}
\keyword{population genetics}
\keyword{drift}
