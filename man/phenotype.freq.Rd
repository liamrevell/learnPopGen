\name{phenotype.freq}
\alias{phenotype.freq}
\alias{phenotype.selection}
\title{Computes phenotypic distribution and its change through time due to natural selection on a polygenic trait}
\usage{
	phenotype.freq(nloci=6, p=NULL, effect=1/nloci)
	phenotype.selection(nloci=6, p=NULL, effect=1/nloci, beta=0.1, ngen=20, ...)
}
\arguments{
	\item{nloci}{number of loci. For simplicity all loci are assumed to be biallelic.}
	\item{p}{allele frequency, \emph{p}, for each locus, in a vector. If not supplied, initially frequencies will be assumed to be 0.5 at all loci.}
	\item{effect}{additive effect of an allele substitution. For simplicity, this is assumed to be the same at all loci.}
	\item{beta}{selection gradient.}
	\item{ngen}{number of generations to analyze.}
	\item{...}{optional arguments.}
}
\description{
	\code{phenotype.freq} computes the phenotypic trait distribution for a polygenic trait. Can be used to demonstrate that the phenotypic distribution of a polygenic trait will tend to normality as the number of loci is increased, regardless of the allele frequencies at each locus.
	
	\code{phenotype.selection} computes the change in the phenotypic trait distribution through time under natural selection. Can be used to show that natural selection on a polygenic trait can move the value of the trait well beyond its original distribution in the population.
}
\value{
	Creates a plot or animation.
	
	\code{phenotype.freq} also invisibly returns an object of class \code{"phenotype.freq"} that can be printed or re-plotted using \code{print} and \code{plot} methods corresponding to the object type. (See examples.)
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{clt}}, \code{\link{selection}}
}
\examples{
	\dontrun{
		phenotype.freq(n=4)
		object<-phenotype.freq(nloci=6,p=runif(n=6),effect=1/6)
		print(object)
		plot(object)
		phenotype.selection(ngen=100)
	}
}
\keyword{population genetics}
