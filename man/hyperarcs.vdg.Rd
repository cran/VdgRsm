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
CCD1<- gen.CCD(n.vars = 3, n.center = 2, alpha = 1)
CCD2<- gen.CCD(n.vars = 3, n.center = 5, alpha = 1)
hyperarcs.vdg(CCD1, CCD2)
}

