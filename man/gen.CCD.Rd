\name{gen.CCD}
\alias{gen.CCD}
\title{
  Central Composite Designs
}
\description{
  Generate central composite designs 
}
\usage{
  gen.CCD(n.vars, n.center, alpha, varNames)
}
\arguments{
  \item{n.vars}{
    The number of factors or independent variables
  }
  \item{n.center}{
    The number of center points
    }
      \item{alpha}{
    The axial distance
    }
          \item{varNames}{
    The variable names. If it is not provided, the default names are X1,X2,...,Xk
    }
}
\value{
  \code{gen.CCD} is called to generate the data frame of the design matrix of the CCD
}
\examples{
  CCD1<- gen.CCD(n.vars = 3, n.center = 2, alpha = 1)
  CCD2<- gen.CCD(n.vars = 3, n.center = 2, alpha = 1, varNames = c("T1","T2","T3"))
}

\keyword{Central Composite Designs}
