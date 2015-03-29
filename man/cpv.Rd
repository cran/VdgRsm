\name{cpv}
\alias{cpv}
\title{
Cuboidal Prediction Variance
}
\description{
Create a variance dispersion graph for a response surface design in a cuboidal region.
}
\usage{
cpv(design.matrix, design.matrix.2 = NULL, des.names = c("Design 1","Design 2"),
    add.pts = TRUE)
}
\arguments{
\item{design.matrix, design.matrix.2}{
  Data frames of design points to be compared in coded or uncoded units. There should be one column for each factor in the design, and one row for each run in the design. The maximum number of factors is 6. If the number of factor is more than 4, only one design is allowed.}
  \item{add.pts}{
 Generate scaled prediction variances of random design points in the VDG. By default \code{add.pts = TRUE}.}
 \item{des.names}{
  A vector of descriptive names for designs in character strings. 
}
}

\value{
  \code{cpv} is called to generate a variance sispersion graph when the number of factors k =  2, 3, or 4 and to generate side-by-side boxplots for k = 5 and 6. In the former case, a table of the minimum, maximum, and average of scaled prediction variances is also produced. 
}
\examples{
CCD1<- gen.CCD(n.vars = 3, n.center = 2, alpha = 1)
CCD2<- gen.CCD(n.vars = 3, n.center = 5, alpha = 1)
cpv(CCD1, CCD2, des.names = c("CCD with nc=2", "CCD with nc=5"), add.pts = FALSE)
}

\keyword{Cuboidal Designs}

