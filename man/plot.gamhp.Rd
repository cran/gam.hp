% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.gamhp.r
\name{plot.gamhp}
\alias{plot.gamhp}
\title{Plot for a \code{\link{gam.hp}} object}
\usage{
\method{plot}{gamhp}(x, plot.perc = FALSE, ...)
}
\arguments{
\item{x}{A \code{\link{gam.hp}} object.}

\item{plot.perc}{Logical;if TRUE, the bar plot (based on ggplot2 package) of the percentage to individual effects of variables towards total explained variation, the default is FALSE to show plot with original individual effects.}

\item{...}{unused}
}
\value{
a ggplot object
}
\description{
Plot for a \code{\link{gam.hp}} object
}
\examples{
library(mgcv)
mod1 <- gam(Sepal.Length ~ s(Petal.Length) + s(Petal.Width) + Sepal.Width,data = iris)
plot(gam.hp(mod1))
}
\author{
{Jiangshan Lai} \email{lai@njfu.edu.cn}
}
