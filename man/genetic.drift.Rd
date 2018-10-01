\name{genetic.drift}
\alias{genetic.drift}
\title{Genetic drift simulation}
\usage{
	genetic.drift(p0=0.5, Ne=20, nrep=10, time=100, show="p", pause=0.1, ...)
}
\arguments{
	\item{p0}{Starting frequency for the A allele.}
	\item{Ne}{Effective population size.}
	\item{nrep}{Number of replicate simulations.}
	\item{time}{Total time, in number of generations, for the simulation.}
	\item{show}{Various options for plotting. \code{"p"} shows the frequency of the A allele through time; \code{"genotypes"} creates an animated histogram with the frequencies of each of the three genotypes through time; \code{"fixed"} shows the fraction of populations that have fixed for each allele, a or A; \code{"heterozygosity"} plots the mean heterozygosity and the expected heterozygosity through time. The default is \code{show="p"}.}
	\item{pause}{Pause between generations. \code{pause=0.01} (for instance) might smooth animation.}
	\item{...}{optional arguments.}
}
\description{
	This function simulates genetic drift at a biallelic genetic locus with no selection and no mutation in a sexually reproducing diploid population or set of populations. It is essentially redundant with \code{\link{drift.selection}}, but in which there is no difference in relative fitness among genotypes; however, it also allows the user to visualize heterozygosity or genetic variation through time - options that are not yet implemented in \code{\link{drift.selection}}.
}
\value{
	The function creates one of several possible plots, depending on the value of \code{show}.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
 	\code{\link{drift.selection}}, \code{\link{founder.event}}, \code{\link{selection}}
}
\examples{
	genetic.drift()
	genetic.drift(p0=0.7,show="heterozygosity")
}
\keyword{population genetics}
\keyword{drift}
