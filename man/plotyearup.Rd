\name{plotyearup}
\alias{plotyearup}
\alias{plot.yearup}
\title{Plot Yearly Observability for a Single Target}
\description{
Produces a two-panel plot showing (i) the number of hours per night the
target is above 30 and 60 degrees altitude across a full year, and (ii)
the moon–target angular separation and moon phase over the same period.
}
\usage{
plotyearup(obsdata)

\method{plot}{yearup}(x, ...)
}
\arguments{
  \item{obsdata}{An object of class \code{yearup} as returned by
    \code{\link{yearup}}.}
  \item{x}{An object of class \code{yearup} (for the S3 \code{plot}
    method).}
  \item{\dots}{Additional arguments passed to \code{plotyearup} (for the
    S3 \code{plot} method; currently unused).}
}
\value{
Called for its side effect (a plot is drawn).  Returns \code{NULL}
invisibly.
}
\seealso{
\code{\link{yearup}}, \code{\link{plotdayup}}
}
\examples{
yr <- yearup(RA = "12:30:16", Dec = "-30:13:15",
             Date = c(2023, 1),
             Lon = 115.8178, Lat = -31.97333,
             UTCdiff = 8)

# S3 method
plot(yr)

# Direct function call
plotyearup(yr)
}
