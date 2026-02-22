\name{plotdayupmulti}
\alias{plotdayupmulti}
\alias{plot.dayupmulti}
\title{Plot Daily Observability for Multiple Targets}
\description{
Produces a two-panel plot showing (i) altitude or air mass against
local time for all targets (in blue), the Moon and the Sun, and (ii)
azimuth against local time.  Each target's peak is labelled with its
name.
}
\usage{
plotdayupmulti(obsdata, ytype = 'Alt', moonphase = TRUE, Names)

\method{plot}{dayupmulti}(x, ...)
}
\arguments{
  \item{obsdata}{An object of class \code{dayupmulti} as returned by
    \code{\link{dayupmulti}}.}
  \item{x}{An object of class \code{dayupmulti} (for the S3 \code{plot}
    method).}
  \item{ytype}{Either \code{'Alt'} (altitude in degrees, the default)
    or \code{'AM'} (air mass).}
  \item{moonphase}{Logical.  If \code{TRUE} (the default) the Moon
    track is coloured according to its phase.}
  \item{Names}{Character vector of target labels (one per target).  If
    absent, uses the names stored in \code{obsdata$names}.}
  \item{\dots}{Additional arguments passed to \code{plotdayupmulti} (for
    the S3 \code{plot} method).}
}
\value{
Called for its side effect (a plot is drawn).  Returns \code{NULL}
invisibly.
}
\seealso{
\code{\link{dayupmulti}}, \code{\link{plotdayup}}, \code{\link{plotyearup}}
}
\examples{
targets_mat <- cbind(
  Names = c("FieldA", "FieldB", "FieldC"),
  RA    = c("02:00:00", "09:30:00", "22:00:00"),
  Dec   = c("-05:00:00", "-10:00:00", "-30:00:00")
)
obs_multi <- dayupmulti(Targets = targets_mat,
                        Date = c(2023, 7, 1), Time = c(18, 0, 0),
                        Lon = 115.8178, Lat = -31.97333,
                        UTCdiff = 8)

# S3 method
plot(obs_multi)

# Direct call
plotdayupmulti(obs_multi, ytype = 'Alt')
}
