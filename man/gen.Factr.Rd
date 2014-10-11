\name{gen.Factr}
\alias{gen.Factr}
\title{
  Factorial Designs
}
\description{
  Generate factorial designs 
}
\usage{
  gen.Factr(n.vars, n.levels, varNames = NULL, scale = TRUE)
}
\arguments{
  \item{n.vars}{
    The number of factors or independent variables
  }
  \item{n.levels}{
    The number of levels of the factor
    }
      \item{varNames}{
    The names of factors. If it is not provided, the default names are X1, X2,..., Xk.
    }
          \item{scale}{
    If it is \code{scale = TRUE}, the level values will be scaled to -1 to 1.
    }
}
\value{
  \code{gen.Factr} is called to generate the data frame of the design matrix of the factorial design.
}
\examples{
CCD1<- gen.Factr(n.vars = 3, n.levels = 5)
CCD2<- gen.Factr(n.vars = 3, n.levels = 5, varNames = c("T1","T2","T3"), scale = FALSE)
}
\keyword{Factorial Designs}

