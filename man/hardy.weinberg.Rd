\name{hardy.weinberg}
\alias{hardy.weinberg}
\alias{multilocus.hw}
\title{Computes Hardy-Weinberg frequencies for a multiallelic locus or across multiple loci}
\usage{
	hardy.weinberg(p=c(0.5,0.5),alleles=c("A","a"),as.matrix=FALSE)
	multilocus.hw(nloci=2,p=NULL)
}
\arguments{
	\item{p}{allele frequencies. In the case of \code{multilocus.hw} the frequencies of the 'dominant' allele at each locus.}
	\item{alleles}{names of the alleles.}
	\item{as.matrix}{logical argument indicating whether to return the result in the form of a matrix (if \code{TRUE}) or a vector.}
	\item{nloci}{for \code{multilocus.hw} the number of loci.}
}
\description{
	\code{hardy.weinberg} computes Hardy-Weinberg frequencies for a multiallelic locus using arbitrary allele frequencies.
	
	\code{multilocus.hw} computes multilocus Hardy-Weinberg frequencies for a set of biallelic loci.
}
\value{
	Returns a matrix or vector.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\examples{
	hardy.weinberg()
	hardy.weinberg(p=c(0.4,0.3,0.2,0.1),alleles=letters[1:4])
}
\keyword{population genetics}
