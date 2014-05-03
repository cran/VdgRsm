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
### Generate 3 FCDs for k = 3 with different numbers of center runs
### library(DoE.wrapper)
FCD3.1 <- ccd.design(nfactors=3,  ncenter = 1, alpha = 1, randomize = FALSE)
FCD3.1 <- FCD3.1[c(-9),-1]
FCD3.2 <- ccd.design(nfactors=3,  ncenter = 3, alpha = 1, randomize = FALSE)
FCD3.2 <- FCD3.2[c(-9,-10,-11),-1]
FCD3.3 <- ccd.design(nfactors=3,  ncenter = 5, alpha = 1, randomize = FALSE)
FCD3.3 <- FCD3.3[c(-9,-10,-11,-12,-13),-1]
### create a FDS plot
#  fds.cube(FCD3.1,FCD3.2,FCD3.3,
#           des.names = c("FCD 1","FCD 2","FCD 3"))
}


