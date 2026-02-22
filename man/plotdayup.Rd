\name{plotdayup}
\alias{plotdayup}
\alias{plot.dayup}
\title{Plot Daily Observability for a Single Target}
\description{
Produces a two-panel plot showing (i) altitude or air mass against
local time for the target, Moon and Sun, and (ii) azimuth against local
time.  Sunrise and sunset times are marked.
}
\usage{
plotdayup(obsdata, ytype = 'Alt', moonphase = TRUE, Name = '')

\method{plot}{dayup}(x, ...)
}
\arguments{
  \item{obsdata}{An object of class \code{dayup} as returned by
    \code{\link{dayup}}.}
  \item{x}{An object of class \code{dayup} (for the S3 \code{plot}
    method).}
  \item{ytype}{Either \code{'Alt'} (plot altitude in degrees, the
    default) or \code{'AM'} (plot air mass).}
  \item{moonphase}{Logical.  If \code{TRUE} (the default) the Moon
    track is coloured according to its illuminated fraction: bright
    white near full moon, dark grey near new moon.}
  \item{Name}{Character string label placed at the peak of the target
    track.}
  \item{\dots}{Additional arguments passed to \code{plotdayup} (for the
    S3 \code{plot} method).}
}
\value{
Called for its side effect (a plot is drawn).  Returns \code{NULL}
invisibly.
}
\seealso{
\code{\link{dayup}}, \code{\link{plotdayupmulti}}, \code{\link{plotyearup}}
}
\examples{
obs <- dayup(RA = "12:30:16", Dec = "-30:13:15",
             Date = c(2023, 6, 1), Time = c(12, 0, 0),
             UTCdiff = 8)

# Using the S3 method
plot(obs)

# Direct call with altitude axis and target label
plotdayup(obs, ytype = 'Alt', Name = "My Target")

# Air-mass version
plotdayup(obs, ytype = 'AM', Name = "My Target")
}
