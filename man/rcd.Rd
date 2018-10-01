\name{rcd}
\alias{rcd}
\title{Simulation of reproductive character displacement in an ecological community}
\usage{
	rcd(nsp=3,nindivs=c(700,400,100),w_t=10,gen=c(500,500),figs="on",pf=100)
}
\arguments{
	\item{nsp}{Number of species in the simulation. If \code{figs="on"}, \code{nsp} must be 1, 2, or 3.}
	\item{nindivs}{A vector of length \code{nsp} containing the integer number of individuals in each species of the simulation.}
	\item{w_t}{Shape parameter of the Gaussian selection surface for the male signalling trait.}
	\item{gen}{Vector containing the number of allopatric generations followed by the number of sympatric generations for simulation.}
	\item{figs}{Either \code{"on"} if plotting is turned on, or \code{"off"} to suppress plotting.}
	\item{pf}{Print frequency for the simulation status to screen.}
}
\description{
	This function conducts individual-based, genetically explicit numerical simulation of reproductive character displacement in an ecological community. The model is one of multiple species (with fixed relative abundance) competing to utilize the same signal space. There is both stabilizing selection on the signal trait for detectability, as well as (in multi-species simulations) countervailing selection for divergence due to the costs of erroneous mismating attempts.
}
\value{
	The function returns a list containing the mean male signal trait and the mean female preference over time. It also (optionally) plots these.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
 	\code{\link{genetic.drift}}, \code{\link{selection}}
}
\examples{
	\dontrun{
		obj<-rcd(nsp=2,nindivs=c(500,500))
	}
}
\keyword{population genetics}
\keyword{drift}
\keyword{sexual selection}
\keyword{character displacement}
