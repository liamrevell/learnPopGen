\name{coalescent.plot}
\alias{coalescent.plot}
\alias{print.coalescent.plot}
\alias{plot.coalescent.plot}
\title{Creates a (usually animated) simulation of gene coalescence within a population}
\usage{
	coalescent.plot(n=10, ngen=20, colors=NULL, ...)
	\method{print}{coalescent.plot}(x, ...)
	\method{plot}{coalescent.plot}(x, ...)
}
\arguments{
	\item{n}{number of haploid individuals or gene copies.}
	\item{ngen}{number of generations.}
	\item{colors}{colors to use for plotting individuals and lines. By default, the function tries to use a contrasting color scheme such that adjacent allele copies are dissimilar (to facilitate visualization of the coalescent process.)}
	\item{x}{object of class \code{"coalescent.plot"} for \code{print} and \code{plot} methods.}
	\item{...}{optional arguments. For \code{coalescent.plot} optional arguments include: \code{sleep}, the time to pause between generations (set to 0.2s by default); \code{lwd}, the line width for parent-offspring lines in the coalescent geneology; and \code{col.order}, if \code{colors=NULL}, whether to use 'sequential' (\code{col.order="sequential"}, the default) or 'alternating' (\code{col.order="alternating"}) colors for adjacent alleles.}
}
\description{
	Coalescence or coalescent theory is a model for genetic drift within a population in which we envision gene copies merging or "coalescing" into ancestors in the past. This function generates a(n) (optionally animated) visualization of this process of coalescence within a population.
}
\value{
	Creates a plot or animation.
	
	Invisibly returns an object of class \code{"coalescent.plot"} containing the alleles (coded numerically) and the parent-offspring relationships from the coalescent simulation. This object can be printed or re-plotted using \code{print} and \code{plot} methods. (See examples.)
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{drift.selection}}, \code{\link{genetic.drift}}
}
\examples{
	coalescent.plot()
	\dontrun{
		coalescent.plot(n=20,ngen=30,col.order="alternating")
		object<-coalescent.plot()
		print(object)
		plot(object)
	}
}
\keyword{population genetics}
