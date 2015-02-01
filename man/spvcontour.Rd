\name{spvcontour}
\alias{spvcontour}
\title{
  Contour Plot of Scaled Prediction Variances 
}
\description{
  Create a contour plot of scaled prediction variances 
}
\usage{
  spvcontour(design.matrix, shape, max.radius = sqrt(2), length = 100, 
  nlevels = 10, title = "Contour of SPVs")
}

\arguments{
  \item{design.matrix}{
    A data frame of design points. There should be one column for each factor in the design, and one row for each run in the design. Only design with 2 factors is allowed.}
  \item{shape}{
    The shape can be "circle" or "square" which represent a shape of design space.
  }
  \item{max.radius}{
    The radius of a circle.
  }
  \item{length}{
    Argument from interp from library \code{akima}.}
  \item{nlevels}{
    Argument from interp from library \code{akima}.}
    \item{title}{
    The title of a contour plot.}
}
\value{
  \code{spvcontour} is called to generate a contour plot of scaled prediction variances for response surface designs.
}
\examples{ 
  library(akima)
  CCD1<- gen.CCD(n.vars = 2, n.center = 2, alpha = 1)
  spvcontour(CCD1, shape = "square")
  CCD2<-gen.CCD(n.var = 2, alpha = sqrt(2), n.center = 3)
  spvcontour(CCD2, shape = "circle")
  spvcontour(CCD2, shape = "circle", length = 200)
  spvcontour(CCD2, shape = "circle", length = 200, nlevels = 20)
}
\keyword{Contour Plot}
