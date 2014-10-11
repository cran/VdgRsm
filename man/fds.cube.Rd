\name{fds.cube}
\alias{fds.cube}
\title{
   The Fraction of Design Space (FDS) plots for cuboidal-region designs
}
\description{
  Create Fraction of Design Space plots for response surface designs in cuboidal regions.
}
\usage{
fds.cube(design.matrix, design.matrix.2 = NULL, design.matrix.3 = NULL, 
         des.names = c("Design 1","Design 2","Design 3"))
}

\arguments{
\item{design.matrix, design.matrix.2, design.matrix.3}{
  Data frames of design points to be compared in coded or uncoded units. There should be one column for each factor in the design, and one row for each run in the design. The maximum number of factors is 6.}
\item{des.names}{
  A vector of descriptive names for designs in character strings. 
}}

\value{
  \code{fds.sphere} is called to generate Fraction of Design Space plots for cuboidal-region designs.
}

\examples{
CCD1<- gen.CCD(n.vars = 4, n.center = 2, alpha = 1)
CCD2<- gen.CCD(n.vars = 4, n.center = 5, alpha = 1)
fds.cube(CCD1, CCD2)
}


