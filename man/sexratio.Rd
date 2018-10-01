\name{sexratio}
\alias{sexratio}
\title{Hypothetical analysis of frequency dependent selection on a sex determining genetic locus}
\usage{
	sexratio(p0=0.01,time=40,show="p",pause=0,sex.Aa=c(0.5,0.5))
}
\arguments{
	\item{p0}{Starting frequency for the \emph{A} allele. Individuals with \emph{AA} genotypes are male, while individuals with \emph{Aa} genotypes are male or female with probability given by \code{sex.Aa}.}
	\item{time}{Number of generations to run the analysis.}
	\item{show}{Two different options for plotting. \code{"p"} shows the frequency of \emph{A} through time; \code{"fitness"} plots the mean fitness through time. The default is \code{show="p"}.}
	\item{pause}{Pause between generations. \code{pause=0.01} (for instance) might smooth animation.}
	\item{sex.Aa}{Probability that individuals with genotype \emph{Aa} are male or female, respectively.}
}
\description{
	This function performs numerical analysis of a frequency dependent selection model of a hypothetical diploid sexually reproducing population in which sex is determined by the genotype at a biallelic genetic locus. Genotype \emph{AA} are male, genotype \emph{aa} are female, and genotype \emph{Aa} might be male or female with probabilities that can be specified by the user. (Users may find, for instance, that setting \code{sex.Aa=c(1,0)} will result in evolution towards an \emph{XY} sex determination system; whereas \code{sex.Aa=c(0,1)} will evolve towards a \emph{ZW} system.)
}
\value{
	The function creates one of two possible plots, depending on the value of \code{show}. Numerical analysis of this model shows how frequency dependent selection should favor alleles that tend to produce the rarer sex in the population.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
  \code{\link{sexratio}}
}
\examples{
	sexratio()
	sexratio(p0=0.001,show="sex-ratio",pause=0.01)
	sexratio(p0=0.001,show="fitness")
	sexratio(p0=0.001,show="sex-ratio",sex.Aa=c(1,0))
}
\keyword{population genetics}
\keyword{evolutionary theory}
