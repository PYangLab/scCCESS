% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scCCESS.R
\name{criteria_compare}
\alias{criteria_compare}
\title{Clustering concordance by ARI, NMI, FM or Jaccard index}
\usage{
criteria_compare(x, y, criteria)
}
\arguments{
\item{x}{clusters A}

\item{y}{clusters B}

\item{criteria}{measures include ARI, NMI}
}
\value{
the similarity value of two clusters.
}
\description{
Clustering concordance by ARI, NMI, FM or Jaccard index
}
\examples{

# criteria_compare(x=c(1,2,2,2,1,2),y=c(1,2,2,3,1,2),"ARI")

}
\keyword{internal}
