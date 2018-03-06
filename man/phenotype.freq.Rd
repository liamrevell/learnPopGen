\name{phenotype.freq}
\alias{phenotype.freq}
\title{Computes phenotypic distribution for polygenic trait}
\usage{
	phenotype.freq(nloci=10,p=NULL,effect=1)
}
\arguments{
	\item{nloci}{number of loci. For simplicity all loci are assumed to be biallelic.}
	\item{p}{allele frequency \emph{p} for each locus.}
	\item{effect}{additive effect of an allele substitution. For simplicity, this is assumed to be the same at all loci.}
}
\description{
	Computes the phenotypic trait distribution for a polygenic trait. Can be used to demonstrate that the phenotypic distribution of a polygenic trait will tend to normality as the number of loci is increased, regardless of the allele frequencies at each locus.
}
\value{
	Creates a plot
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\examples{
	phenotype.freq(n=4)
	phenotype.freq(nloci=6,p=runif(n=6),effect=1/6)
}
\keyword{population genetics}
