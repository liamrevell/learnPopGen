\name{mutation.selection}
\alias{mutation.selection}
\title{Gene frequencies over time under mutation-selection balance}
\usage{
	mutation.selection(p0=1.0, w=c(1,1), u=0.001, time=100, show="q", pause=0, 
		ylim=c(0,1))
}
\arguments{
	\item{p0}{Starting frequency for the \emph{A} allele.}
	\item{w}{Fitnesses of the heterozygote (\emph{Aa}) and homozygote deleterious (\emph{aa}) genotypes. The fitness of genotype \emph{AA} is assumed to be 1.0.}
	\item{u}{Rate at which \emph{A} alleles are converted to \emph{a} alleles by mutation.}
	\item{time}{Number of generations to run the analysis.}
	\item{show}{Two options for plotting. \code{"q"} shows the frequency of a through time; \code{"fitness"} plots the mean fitness over time. The default is \code{show="q"}.}
	\item{pause}{Pause between generations. \code{pause=0.01} (for instance) might smooth animation.}
	\item{ylim}{Limits on the \emph{y}-axis for plotting.}
}
\description{
	This function performs numerical analysis of mutation-selection balance with mutation from \emph{A} to \emph{a} and selection against (either or both of) \emph{Aa} and \emph{aa}.
}
\value{
	The function creates one of three possible plots, depending on the value of \code{show}.
	
	The function also invisibly returns the frequency of the A allele through time and the mean population fitness as an object of class \code{"mutation.selection"} that can be printed or re-plotted with associated \code{print} and \code{plot} methods, respectively. The \code{plot} method also permits user control over various attributes of the appearance of the plot, such as the color of the plotted lines (\code{color}), the line widths (\code{lwd}), the limits of the \emph{y}-axis (\code{ylim}), and the type of line (e.g., \code{"l"} vs. \code{"s"}, via the argument \code{type}).
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{genetic.drift}}, \code{\link{selection}}
}
\examples{
	mutation.selection(w=c(1,0),time=100,ylim=c(0,0.1))
}
\keyword{population genetics}
\keyword{evolutionary theory}
