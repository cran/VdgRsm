\name{gen.USD}
\alias{gen.USD}
\title{
   Doehlert's Uniform Shell Designs (USDs)
}
\description{
  Generate uniform shell designs for k = 2 to 6
}
\usage{
  gen.USD(k, alpha = 1)
}
\arguments{
  \item{k}{
    The number of factors or independent variables, k = 2 to 6.
  }
  \item{alpha}{
    A scaling factor. See examples.
    }
}
\value{
  \code{gen.USD} is called to generate the data frame of the design matrix of the USD.
}
\examples{
  gen.USD(2)
  gen.USD(3, alpha = sqrt(3))
  gen.USD(6)
  gen.USD(6, alpha = sqrt(6))
}
\references{
Doehlert, D. H (1970), \emph{Uniform Shell Designs}, Journal of the Royal Statistical Society, 19(3):231-239 
}
\keyword{Doehlert's Uniform Shell Designs}
