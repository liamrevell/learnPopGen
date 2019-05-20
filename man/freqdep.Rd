\name{freqdep}
\alias{freqdep}
\title{Numerical analysis of a frequency dependent seletion model}
\usage{
	freqdep(p0=0.01, s=0, time=100, show="p", pause=0, ...)
}
\arguments{
	\item{p0}{Starting frequency for the \emph{A} allele.}
	\item{s}{Parameter that determines the strength of selection against heterozygotes when they are common.}
	\item{time}{Number of generations to run the analysis.}
	\item{show}{Various options for plotting. \code{"p"} shows the frequency of A through time; \code{"surface"} plots the mean fitness as a function of p; \code{"deltap"} shows the change in p as a function of p; \code{"cobweb"} creates a cobweb plot showing p(t) by p(t+1). The default is \code{show="p"}.}
	\item{pause}{Pause between generations. \code{pause=0.01} (for instance) might smooth animation.}
	\item{...}{optional arguments. Presently, the only optional argument in \code{freqdep} is \code{color}, which can be used to change the color of the lines of the plot. The \code{plot} method can also accept the optional arguments \code{type} (e.g., \code{"l"} or \code{"s"}) and \code{lwd}.}
}
\description{
	This function performs numerical analysis of a frequency dependent selection model based on Rice (2004; \emph{Evolutionary Theory: Mathematical & Conceptual Foundations}). The fitnesses of the three genotypes in the model are as follows, where f(\emph{XX}) denotes the frequency of the \emph{XX} genotype: w(\emph{AA})=1-3*f(\emph{Aa})+3*f(\emph{aa}); w(\emph{Aa})=1-s*f(\emph{Aa}); and w(\emph{aa})=1-3*f(\emph{Aa})+3*f(\emph{AA}). As shown in Rice (2004), though entirely deterministic, the model can exhibit chaotic behavior under some values for \code{s}.
}
\value{
	The function creates one of several possible plots, depending on the value of \code{show}.
	
	The function also invisibly returns an object of class \code{"freqdep"} containing the frequency of the allele A through time, if this was calculated by the selected method. This can be printed or plotted using the corresponding methods. (See examples.)
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
  \code{\link{sexratio}}
}
\examples{
	freqdep(time=100)
	freqdep(s=1.5,time=100)
	p<-freqdep(s=2,show="cobweb",time=100)
	plot(p)
}
\keyword{population genetics}
\keyword{evolutionary theory}
