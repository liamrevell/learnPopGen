\name{sexratio}
\alias{sexratio}
\title{Hypothetical analysis of frequency dependent selection on a sex determining genetic locus}
\usage{
	sexratio(p0=0.01,time=100,show="p",pause=0)
}
\arguments{
	\item{p0}{Starting frequency for the A allele. Individuals with AA genotypes are male, while individuals with Aa genotypes are male or female with equal probability.}
	\item{time}{Number of generations to run the analysis.}
	\item{show}{Two different options for plotting. \code{"p"} shows the frequency of A through time; \code{"fitness"} plots the mean fitness through time. The default is \code{show="p"}.}
	\item{pause}{Pause between generations. \code{pause=0.01} (for instance) might smooth animation.}
}
\description{
	This function performs numerical analysis of a frequency dependent selection model of a hypothetical diploid sexually reproducing population in which sex is determined by the genotype at a biallelic genetic locus. Genotype AA are male, genotype aa are female, and genotype Aa might be male or female with equal probability.
}
\value{
	The function creates one of two possible plots, depending on the value of \code{show}. Numerical analysis of this model shows how frequency dependent selection should favor alleles that tend to produce the rarer sex in the population.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
  \code{\link{sexratio}}
}
\examples{
	sexratio(time=100)
}
\keyword{population genetics}
\keyword{evolutionary theory}
