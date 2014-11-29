\name{gen.BBD}
\alias{gen.BBD}
\title{
   Box-Behnken Designs (BBDs)
}
\description{
  Generate Box-Behnken Designs for k = 3 to 7
}
\usage{
  gen.BBD(k, n.center = 1)
}
\arguments{
  \item{k}{
    The number of factors or independent variables, k = 3 to 7.
  }
  \item{n.center}{
    The number of center points
    }
}
\value{
  \code{gen.BBD} is called to generate the data frame of the design matrix of the BBD.
}
\examples{
  gen.BBD(3)
  gen.BBD(4, n.center = 3)
  gen.BBD(7, n.center = 5)
}
\references{
SAS 9.1 ADX Interface for Design of Experiments. Cary, NC: SAS Institute Inc.
}
\keyword{Box-Behnken Designs}
