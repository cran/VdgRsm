\name{spv}
\alias{spv}
\title{
Spherical Prediction Variance
}
\description{
Create Variance Dispesion Graphs (VDGs) for response surface designs in spherical regions.
}
\usage{
spv(design.matrix, design.matrix.2 = NULL, design.matrix.3 = NULL, 
    des.names = c("Design 1","Design 2","Design 3"),  
    scale = TRUE, add.pts = TRUE)
}

\arguments{
\item{design.matrix, design.matrix.2, design.matrix.3}{
  Data frames of design points to be compared in coded or uncoded units. There should be one column for each factor in the design, and one row for each run in the design. The maximum number of factors is 7.}
\item{des.names}{
  A vector of descriptive names for designs in character strings. 
}
\item{scale}{
  Design points are scaled by a factor equal to the square root of the number of factors divided by the maximum of radii across the set of design points. This factor makes two or more designs comparable by scaling the maximum design point radius to be the square root of the number of factors.
}
\item{add.pts}{
  Generate scaled prediction variances of random design points in the VDG. By default \code{add.pts = TRUE}}.
}
\value{
  \code{spv} is called to generate the Variance Dispersion Graph(s) and a table of the minimum, maximum, and average of scaled prediction variances for each design.
}
\examples{
CCD1<- gen.CCD(n.vars = 3, n.center = 2, alpha = 1)
CCD2<- gen.CCD(n.vars = 3, n.center = 2, alpha = sqrt(3))
spv(CCD1, CCD2, des.names = c("CCD 1","CCD 2"))
}
\keyword{Spherical Designs}

