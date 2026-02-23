\name{nowup}
\alias{nowup}
\title{Find Targets Currently Observable}
\description{
Queries the built-in \code{\link{targets}} catalogue and returns those
sources that fall within the requested altitude, azimuth, and type
constraints at the specified (or current) date and time.  Results are
sorted by decreasing altitude.
}
\usage{
nowup(Ncut = 20, Azlim = c(0, 360), Altlim = c(60, 90),
      Date = 'get', Time = 'get',
      Lon = 115.8178, Lat = -31.97333, Loc = 'user',
      UTCdiff = 'get', Altitude = 10, Pressure = 1000,
      Temp = 20, select = c('G', 'S', 'C'))
}
\arguments{
  \item{Ncut}{Maximum number of targets to return (sorted by decreasing
    altitude).}
  \item{Azlim}{Length-2 numeric vector giving the azimuth window in
    degrees \code{c(min, max)}.  If \code{Azlim[1] > Azlim[2]} the
    window wraps around North (e.g. \code{c(300, 60)}).}
  \item{Altlim}{Length-2 numeric vector giving the altitude window in
    degrees \code{c(min, max)}.}
  \item{Date}{Either \code{'get'} (use current system date) or a
    length-3 numeric vector \code{c(year, month, day)}.}
  \item{Time}{Either \code{'get'} (use current system time) or a
    length-3 numeric vector \code{c(hour, min, sec)} in local time.}
  \item{Lon}{Observer longitude in decimal degrees (positive east).
    Ignored when \code{Loc} is not \code{'user'}.}
  \item{Lat}{Observer latitude in decimal degrees (positive north).
    Ignored when \code{Loc} is not \code{'user'}.}
  \item{Loc}{Name of a known observatory.  See \code{\link{gettelescope}}
    for valid names.  Use \code{'user'} (the default) to supply
    coordinates via \code{Lon} and \code{Lat}.}
  \item{UTCdiff}{UTC offset of local time in hours.  Use \code{'get'}
    (default) to determine automatically, or \code{'guess'} to estimate
    from longitude.}
  \item{Altitude}{Observatory altitude above sea level in metres.}
  \item{Pressure}{Atmospheric pressure in millibars.}
  \item{Temp}{Atmospheric temperature in degrees Celsius.}
  \item{select}{Character vector of target types to include.  Standard
    codes are \code{'G'} (galaxy), \code{'S'} (star), \code{'C'}
    (cluster).  Must match the \code{Type} column of \code{\link{targets}}.}
}
\value{
A data frame (a subset of the \code{targets} catalogue) with extra
columns:
\describe{
  \item{HA}{Hour angle in hours.}
  \item{Alt}{Altitude in degrees.}
  \item{Az}{Azimuth in degrees.}
}
Rows are sorted by decreasing \code{Alt}, and at most \code{Ncut} rows
are returned.
}
\seealso{
\code{\link{targets}}, \code{\link{nowwhere}}, \code{\link{dayup}}
}
\examples{
# Find galaxy and star targets above 30 degrees altitude
# at a fixed date/time from Perth, Australia
visible <- nowup(Ncut = 10, Altlim = c(30, 90),
                 Date = c(2023, 6, 15), Time = c(22, 0, 0),
                 Lon = 115.8178, Lat = -31.97333,
                 UTCdiff = 8, select = c('G', 'S'))
print(visible[, c("Name", "RA", "Dec", "Alt", "Az")])

# Use a named telescope location
visible2 <- nowup(Ncut = 5, Altlim = c(40, 90),
                  Date = c(2023, 9, 1), Time = c(23, 0, 0),
                  Loc = "AAT", UTCdiff = 10)
print(visible2[, c("Name", "Alt", "Az")])
}
