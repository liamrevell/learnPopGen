\name{selection}
\alias{selection}
\title{Numerical analysis of biallelic locus frequency independent selection}
\usage{
	selection(p0=0.01,w=c(1.0,0.9,0.8),time=100,show="p",pause=0,...)
}
\arguments{
	\item{p0}{Starting frequency for the \emph{A} allele.}
	\item{w}{Fitnesses for the three genotypes in the following order: \emph{AA}, \emph{Aa}, \emph{aa}.}
	\item{time}{Number of generations to run the analysis.}
	\item{show}{Various options for plotting. \code{"p"} shows the frequency of \emph{A} through time; \code{"surface"} plots the mean fitness as a function of \emph{p}; \code{"deltap"} shows the change in \emph{p} as a function of \emph{p}; \code{"cobweb"} creates a cobweb plot showing \emph{p}(\emph{t}) by \emph{p}(\emph{t}+1). The default is \code{show="p"}.}
	\item{pause}{Pause between generations. \code{pause=0.01} (for instance) might smooth animation.}
	\item{...}{Optional arguments, including: \code{add}, a logical value indicating whether or not to add to the current plot (applies only to \code{show="p"}); \code{color}, change the color of the plotted line (works nicely with \code{add}, for obvious reasons); and \code{equil}, a logical value indicating whether or not to show the equilibrium value of \emph{p} using vertical (or horizontal) lines on the graph (defaults to \code{FALSE}).}
}
\description{
	This function performs numerical analysis of a simple biallelic selection model.
}
\value{
	The function creates one of several possible plots, depending on the value of \code{show}.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
  \code{\link{freqdep}}
}
\examples{
	selection(w=c(1.0,0.8,0.8),time=500)
	selection(w=c(1.0,1.0,0.0),show="surface")
}
\keyword{evolutionary theory}
\keyword{population genetics}
