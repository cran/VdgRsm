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
    scale = TRUE, add.pts = TRUE, big.RAM = FALSE)
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
\item{big.RAM}{
  When the number of factors equals to 6 or 7, generated random points will not be generated unless \code{big.RAM = TRUE}. However, to save a memory by default \code{big.RAM = FALSE}.}
}
\value{
  \code{spv} is called to generate the Variance Dispersion Graph(s) and a table of the minimum, maximum, and average of scaled prediction variances for each design.
}
\examples{
### Generate 3 CCDs for k = 3 with different alpha's for the axial points
### library(DoE.wrapper) 
CCD3.1 <- ccd.design(nfactors=3,  ncenter = 2, alpha = 1.5, randomize = FALSE)
CCD3.1 <- CCD3.1[c(-9,-10),-1]
CCD3.2 <- ccd.design(nfactors=3,  ncenter = 2, alpha = 1.0, randomize = FALSE)
CCD3.2 <- CCD3.2[c(-9,-10),-1]
CCD3.3 <- ccd.design(nfactors=3,  ncenter = 2, alpha = sqrt(3), randomize = FALSE)
CCD3.3 <- CCD3.3[c(-9,-10),-1]
### create VDGs
spv(CCD3.1,CCD3.2,CCD3.3,
    des.names = c("CCD, alpha = 1.5","CCD, alpha = 1","CCD, alpha = sqrt(3)"))
}
\keyword{Spherical Designs}

