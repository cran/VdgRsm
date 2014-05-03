\name{fds.sphere}
\alias{fds.sphere}
\title{
  The Fraction of Design Space (FDS) plots for spherical-region designs
}
\description{
  Create Fraction of Design Space plots for response surface designs in spherical regions.
}
\usage{
fds.sphere(design.matrix, design.matrix.2 = NULL, design.matrix.3 = NULL, 
           des.names = c("Design 1","Design 2","Design 3"), scale = TRUE)
}
\arguments{
\item{design.matrix, design.matrix.2, design.matrix.3}{
 Data frames of design points to be compared in coded or uncoded units. There should be one column for each factor in the design, and one row for each run in the design. The maximum number of factors is 7.
}
\item{des.names}{
  A vector of descriptive names for designs in character strings. 
}
\item{scale}{
   Design points are scaled by a factor equal to the square root of the number of factors divided by the maximum of radii across the set of design points. This factor makes two or more designs comparable by scaling the maximum design point radius to be the square root of the number of factors.
}}
\value{
  \code{fds.sphere} is called to generate Fraction of Design Space plots for spherical-region designs.
}
\examples{
### Generate 3 CCDs for k = 5 with different alpha's
### library(DoE.wrapper) 
CCD5.1 <- ccd.design(nfactors=5,  ncenter = 2, alpha = 1.5, randomize = FALSE)
CCD5.1 <- CCD5.1[c(-17,-18),-1]
CCD5.2 <- ccd.design(nfactors=5,  ncenter = 2, alpha = 2, randomize = FALSE)
CCD5.2 <- CCD5.2[c(-17,-18),-1]
CCD5.3 <- ccd.design(nfactors=5,  ncenter = 2, alpha = sqrt(5), randomize = FALSE)
CCD5.3 <- CCD5.3[c(-17,-18),-1]
### create a FDS plot
fds.sphere(CCD5.1,CCD5.2,CCD5.3,
    des.names = c("CCD, alpha = 1.5","CCD, alpha = 2","CCD, alpha = sqrt(5)"))
}
