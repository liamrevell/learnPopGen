\name{freqdep}
\alias{freqdep}
\title{Numerical analysis of a frequency dependent seletion model}
\usage{
	freqdep(p0=0.01,s=0,time=100,show="p",pause=0)
}
\arguments{
	\item{p0}{Starting frequency for the A allele.}
	\item{s}{Parameter that determines the strength of selection against heterozygotes when they are common.}
	\item{time}{Number of generations to run the analysis.}
	\item{show}{Various options for plotting. \code{"p"} shows the frequency of A through time; \code{"surface"} plots the mean fitness as a function of p; \code{"deltap"} shows the change in p as a function of p; \code{"cobweb"} creates a cobweb plot showing p(t) by p(t+1). The default is \code{show="p"}.}
	\item{pause}{Pause between generations. \code{pause=0.01} (for instance) might smooth animation.}
}
\description{
	This function performs numerical analysis of a frequency dependent selection model based on Rice (\emph{Evolutionary Theory}). The fitnesses of the three genotypes in the model are as follows, where f(\emph{XX}) denotes the frequency of the \emph{XX} genotype: w(AA)=1-3*f(Aa)+3*f(aa); w(Aa)=1-s*f(Aa); w(aa)=1-3*f(Aa)+3*f(AA).
}
\value{
	The function creates one of several possible plots, depending on the value of \code{show}.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
  \code{\link{sexratio}}
}
\examples{
	freqdep(time=100)
	freqdep(s=1.5,time=100)
	freqdep(s=2,show="cobweb",time=100)
}
\keyword{population genetics}
\keyword{evolutionary theory}
