\name{gen.Roquemore}
\alias{gen.Roquemore}
\title{
  Roquemore's Hybrid Designs
}
\description{
  Generate Roquemore (1976) hybrid designs for k = 3, 4, and 6. For k = 3, R310, R311A, and R311B will be produced, for k = 4, R416A, R416B, and R416C will be generated, and for k = 6 R628A and R628B will be given. 
}
\usage{
  gen.Roquemore(k, n.center = 0)
}
\arguments{
  \item{k}{
    The number of factors or independent variables, k = 3, 4, and 6.
  }
      \item{n.center}{
    The number of center points
    }
}
\value{
  \code{gen.Roquemore} will retrieve the hybrid design points stored and the output is a list containing relevant Roquemore's designs given a k value. 
}
\examples{
 gen.Roquemore(3)
 gen.Roquemore(4, n.center = 2)
 gen.Roquemore(6, n.center = 1)
}
\references{
SAS 9.1 ADX Interface for Design of Experiments. Cary, NC: SAS Institute Inc.
}
\keyword{Roquemore's Hybrid Designs}
