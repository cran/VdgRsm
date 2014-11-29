\name{gen.PBCD}
\alias{gen.PBCD}
\title{
  Plackett-Burman Composite Designs (Draper and Lin's Method)
}
\description{
  Generate Plackett-Burman composite designs proposed by Draper and Lin (1990) for k = 3 to 7
}
\usage{
  gen.PBCD(k, alpha ="rotatable", n.center = 0)
}
\arguments{
  \item{k}{
    The number of factors or independent variables, k = 3 to 7.
  }
   \item{alpha}{
    Axial distance. User may specify "rotatanble", "face-center", or other numeric numbers. See examples.
    }
      \item{n.center}{
    The number of center points
    }
}
\value{
  \code{gen.PBCD} is called to generate the data frame of the design matrix of the PBCD
}
\examples{
 gen.PBCD(3)
 gen.PBCD(4, alpha = 1)
 gen.PBCD(5, alpha = "face-center", n.center = 3)
 gen.PBCD(6, alpha = 2, n.center = 5)
}
\references{
SAS 9.1 ADX Interface for Design of Experiments. Cary, NC: SAS Institute Inc.
}
\keyword{Plackett-Burman Composite Designs}
