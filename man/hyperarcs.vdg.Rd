\name{hyperarcs.vdg}
\alias{hyperarcs.vdg}
\title{
Scaled prediction variance in nested cubes
}
\description{
Create a graph of scaled prediction variances for points in nested cubes (hyperarcs)
}
\usage{
hyperarcs.vdg(design.matrix, design.matrix.2 = NULL, design.matrix.3 = NULL, 
           des.names = c("Design 1","Design 2","Design 3"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{design.matrix, design.matrix.2, design.matrix.3}{
  Data frames of design points to be compared in coded or uncoded units. There should be one column for each factor in the design, and one row for each run in the design. The minimum and maximum number of factors are 3 and 6, respectively.}
\item{des.names}{
  A vector of descriptive names for designs in character strings. 
}
}
\value{
  \code{hyperarcs.vdg} is called to generate a plot of scaled prediction variances on hyperarcs.
}

\examples{
# library(DoE.wrapper)
FCD3.1 <- ccd.design(nfactors=3,  ncenter = 1, alpha = 1, randomize = FALSE)
FCD3.1 <- FCD3.1[c(-9),-1]
# hyperarcs.vdg(FCD3.1, des.names = c("FCD 1"))
}

