\name{founder.event}
\alias{founder.event}
\title{Simulation of a founder event or population bottleneck}
\usage{
	founder.event(p0=0.5, Ne=1000, Nf=10, ttime=100, etime=50, show="p", ...)
}
\arguments{
	\item{p0}{Starting frequency for the \emph{A} allele.}
	\item{Ne}{Effective population size at the start of the simulation and after the founding event.}
	\item{Nf}{Size of the founding population.}
	\item{ttime}{Total time of the simulation.}
	\item{etime}{Time for the founding event. Can either be a single generation, or a sequence of generations (e.g., \code{etime=40:50}) for a prolonged founder event. Note that to simulate a prolonged bottleneck the user \emph{must} supply the sequence of generations during which the population is to be bottlenecked and \emph{not} merely the start and end times of the bottleneck.}
	\item{show}{Two different options for plotting. \code{"p"} shows the frequency of the A allele through time; \code{"var"} shows the genetic variation in the population, calculated as p*(1-p). The default is \code{show="p"}.}
	\item{...}{optional arguments. Presently, the only optional argument in \code{founder.event} is \code{ltype} which specifies the line type and defaults to \code{"s"}.}
}
\description{
	This function simulates genetic drift with a founding event / population bottleneck at time \code{etime}.
}
\value{
	The function creates one of two different plots, depending on the value of \code{show}.
	
	The function also invisibly returns an object of class \code{"founder.event"} which can be printed or plotted using corresponding \code{print} and \code{plot} methods. (See examples.)
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
  \code{\link{drift.selection}}, \code{\link{genetic.drift}}
}
\examples{
	founder.event()
	p<-founder.event(show="variation")
	plot(p,show="p",ltype="l")
}
\keyword{population genetics}
\keyword{drift}
