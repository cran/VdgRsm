\name{cpv}
\alias{cpv}
\title{
Cuboidal Prediction Variance
}
\description{
Create a Variance Dispesion Graph for a response surface design in a cuboidal region.
}
\usage{
cpv(design.matrix, add.pts = TRUE)
}
\arguments{
  \item{design.matrix}{
  A data frame containing design points in coded or uncoded units. There should be one column for each factor in the design, and one row for each run in the design. The maximum number of factors is 6.
}
  \item{add.pts}{
 Generate scaled prediction variances of random design points in the VDG. By default \code{add.pts = TRUE}.}
}

\value{
  \code{cpv} is called to generate a Variance Dispersion Graph when the number of factors k =  2, 3, or 4 and to generate side-by-side boxplots for k = 5 and 6. In the former case, a table of the minimum, maximum, and average of scaled prediction variances is also produced. 
}
\examples{
CCD1<- gen.CCD(n.vars = 4, n.center = 2, alpha = 1)
cpv(CCD1)
}

\keyword{Cuboidal Designs}

