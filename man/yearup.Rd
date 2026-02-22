\name{yearup}
\alias{yearup}
\title{Compute Yearly Observability for a Single Target}
\description{
Computes how many hours per night a target is above 30 or 60 degrees of
altitude (during astronomical night) over a full year, sampled every
5 days.  Also records moon phase and moon–target angular separation.
Returns an object of class \code{yearup} suitable for plotting with
\code{\link{plot.yearup}} / \code{\link{plotyearup}}.
}
\usage{
yearup(RA = "12:30:16", Dec = "-30:13:15", Target = 'user',
       Date = 'get', Lon = 115.8178, Lat = -31.97333,
       Loc = 'user', UTCdiff = 'get', Altitude = 10,
       Pressure = 1000, Temp = 20, select = c('G', 'S'))
}
\arguments{
  \item{RA}{Character string giving the Right Ascension of the target in
    \code{"HH:MM:SS"} format.  Ignored when \code{Target} is not
    \code{'user'}.}
  \item{Dec}{Character string giving the Declination of the target in
    \code{"sDD:MM:SS"} format.  Ignored when \code{Target} is not
    \code{'user'}.}
  \item{Target}{Name of a known target or solar-system body (e.g.
    \code{'mars'}, \code{'moon'}, \code{'sun'}, or any name in the
    \code{targets} dataset).  Use \code{'user'} (the default) to supply
    coordinates via \code{RA} and \code{Dec}.}
  \item{Date}{Either \code{'get'} (use the current system date) or a
    length-2 numeric vector \code{c(year, month)} giving the start of the
    year-long window.}
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
  \item{select}{Target type filter passed to the internal catalogue
    look-up (not used when \code{Target = 'user'}).}
}
\value{
A list of class \code{yearup} with components:
\describe{
  \item{obs}{Data frame with one row per 5-day sample containing columns
    \code{jd} (Julian Date), \code{LTPOSIX} (POSIXct date), \code{up30}
    (number of 10-min intervals above 30 deg at night), \code{up60}
    (number of 10-min intervals above 60 deg at night), \code{moonsep}
    (moon–target angular separation in degrees), and \code{moonphase}
    (fractional lunar illumination).}
  \item{at}{POSIXct sequence of monthly tick positions for the x-axis.}
  \item{startyear}{Integer: the starting calendar year.}
  \item{yearlines}{POSIXct vector of first-of-month positions spanning
    the two-year window, for drawing month boundaries.}
}
Note: \code{up30} and \code{up60} are in units of 10-minute intervals;
divide by 6 to convert to hours.
}
\seealso{
\code{\link{plotyearup}}, \code{\link{dayup}}, \code{\link{nowup}}
}
\examples{
# Yearly observability for a fixed target from Perth, Australia
yr <- yearup(RA = "12:30:16", Dec = "-30:13:15",
             Date = c(2023, 1),
             Lon = 115.8178, Lat = -31.97333,
             UTCdiff = 8)

# Median observable hours per night above 30 deg
cat("Median hours above 30 deg:",
    round(median(yr$obs$up30) / 6, 1), "\n")

# Median observable hours per night above 60 deg
cat("Median hours above 60 deg:",
    round(median(yr$obs$up60) / 6, 1), "\n")

# From a named observatory
yr2 <- yearup(RA = "05:34:32", Dec = "+22:00:52",
              Date = c(2023, 6), Loc = "AAT", UTCdiff = 10)
cat("Max hours above 30 deg from AAT:",
    round(max(yr2$obs$up30) / 6, 1), "\n")
}
