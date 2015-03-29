\name{gen.HSCD}
\alias{gen.HSCD}
\title{
  Hartley's Small Composite Designs (HSCDs)
}
\description{
  Generate Hartley's small composite designs for k = 2 to 7
}
\usage{
  gen.HSCD(k, alpha ="rotatable", n.center = 0)
}
\arguments{
  \item{k}{
    The number of factors or independent variables, k = 2, 3, 4, 5, 6, and 7.
  }
  \item{alpha}{
    Axial distance. User may specify "rotatable", "face-center", or other numeric numbers. See examples.
    }
      \item{n.center}{
    The number of center points
    }
}
\value{
  \code{gen.HSCD} is called to generate the data frame of the design matrix of the HSCD
}
\examples{
  gen.HSCD(3)
  gen.HSCD(4, alpha ="face-center")
  gen.HSCD(7, alpha = 2, n.center = 4)
}
\references{
SAS 9.1 ADX Interface for Design of Experiments. Cary, NC: SAS Institute Inc.
}
\keyword{Hartley's Small Composite Designs}
