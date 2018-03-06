\name{hardy.weinberg}
\alias{hardy.weinberg}
\title{Computes Hardy-Weinberg frequencies for a multiallelic locus}
\usage{
	hardy.weinberg(p=c(0.5,0.5),alleles=c("A","a"),as.matrix=FALSE)
}
\arguments{
	\item{p}{allele frequencies.}
	\item{alleles}{names of the alleles.}
	\item{as.matrix}{logical argument indicating whether to return the result in the form of a matrix (if \code{TRUE}) or a vector.}
}
\description{
	Computes Hardy-Weinberg frequencies for a multiallelic locus using arbitrary allele frequencies.
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
