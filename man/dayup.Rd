\name{dayup}
\alias{dayup}
\title{Compute Daily Observability for a Single Target}
\description{
Computes the altitude, azimuth, hour angle and airmass of a target (and the Moon
and Sun) over the course of a single day, from a given observer location.  Returns
an object of class \code{dayup} that can be passed to \code{\link{plot.dayup}} /
\code{\link{plotdayup}}.
}
\usage{
dayup(RA = "12:30:16", Dec = "-30:13:15", Target = 'user',
      Date = 'get', Time = c(12, 0, 0),
      Lon = 115.8178, Lat = -31.97333, Loc = 'user',
      UTCdiff = 'get', Altitude = 10, Pressure = 1000,
      Temp = 20, step = 0.1)
}
\arguments{
  \item{RA}{Character string giving the Right Ascension of the target in
    \code{"HH:MM:SS"} format, e.g. \code{"12:30:16"}.  Ignored when
    \code{Target} is not \code{'user'}.}
  \item{Dec}{Character string giving the Declination of the target in
    \code{"sDD:MM:SS"} format, e.g. \code{"-30:13:15"}.  Ignored when
    \code{Target} is not \code{'user'}.}
  \item{Target}{Name of a known target or solar-system body (e.g.
    \code{'mars'}, \code{'moon'}, \code{'sun'}, or any name found in the
    \code{targets} dataset).  Use \code{'user'} (the default) to supply
    coordinates via \code{RA} and \code{Dec}.}
  \item{Date}{Either \code{'get'} (use the current system date) or a
    length-3 numeric vector \code{c(year, month, day)}.}
  \item{Time}{Either \code{'get'} (use the current system time) or a
    length-3 numeric vector \code{c(hour, min, sec)} in local time.}
  \item{Lon}{Observer longitude in decimal degrees (positive east).
    Ignored when \code{Loc} is not \code{'user'}.}
  \item{Lat}{Observer latitude in decimal degrees (positive north).
    Ignored when \code{Loc} is not \code{'user'}.}
  \item{Loc}{Name of a known observatory (e.g. \code{'AAT'}).  See
    \code{\link{gettelescope}} for valid names.  Use \code{'user'} (the
    default) to supply coordinates via \code{Lon} and \code{Lat}.}
  \item{UTCdiff}{Offset of local time from UTC in hours.  Use
    \code{'get'} (default) to determine automatically from the system
    clock, or \code{'guess'} to estimate from longitude.}
  \item{Altitude}{Observatory altitude above sea level in metres.}
  \item{Pressure}{Atmospheric pressure in millibars (used for
    refraction corrections).}
  \item{Temp}{Atmospheric temperature in degrees Celsius (used for
    refraction corrections).}
  \item{step}{Time step in hours between successive calculations over
    the 24-hour period.}
}
\value{
A list of class \code{dayup} with the following components:
\describe{
  \item{obs}{Data frame with one row per time step containing columns:
    \code{JD}, \code{LST}, \code{LT.*} (local time components), \code{LTPOSIX}
    (POSIXct local time), \code{Alt}, \code{Az}, \code{HA}, \code{AirMass},
    \code{AltMoon}, \code{AzMoon}, \code{HAMoon}, \code{AirMassMoon},
    \code{AltSun}, \code{AzSun}, \code{HASun}, \code{AirMassSun}.}
  \item{UTCdiff}{UTC offset used (hours).}
  \item{RA}{Right Ascension string used.}
  \item{Dec}{Declination string used.}
  \item{RAdeg}{Right Ascension in decimal degrees.}
  \item{Decdeg}{Declination in decimal degrees.}
  \item{Lon}{Observer longitude used (decimal degrees).}
  \item{Lat}{Observer latitude used (decimal degrees).}
  \item{Altitude}{Observatory altitude used (metres).}
  \item{Pressure}{Pressure used (millibars).}
  \item{Temp}{Temperature used (degrees Celsius).}
  \item{at}{POSIXct sequence of hourly tick positions for the x-axis.}
  \item{dayline}{POSIXct giving midnight of the observation date.}
  \item{rise}{Two-element POSIXct vector of sunrise times.}
  \item{set}{Two-element POSIXct vector of sunset times.}
  \item{moonsep}{Angular separation between target and Moon (degrees).}
  \item{moonphase}{Fractional lunar illumination (0 = new, 1 = full).}
  \item{sunsep}{Angular separation between target and Sun (degrees).}
}
}
\seealso{
\code{\link{dayupmulti}}, \code{\link{plotdayup}}, \code{\link{yearup}},
\code{\link{nowwhere}}
}
\examples{
# Compute observability for a user-supplied RA/Dec on a fixed date
# from the default location (Perth, Western Australia)
obs <- dayup(RA = "12:30:16", Dec = "-30:13:15",
             Date = c(2023, 6, 1), Time = c(12, 0, 0),
             UTCdiff = 8)
cat("Moon phase:", round(obs$moonphase, 2), "\n")
cat("Moon separation:", round(obs$moonsep, 1), "deg\n")
cat("Sun separation:", round(obs$sunsep, 1), "deg\n")
cat("Max altitude:", round(max(obs$obs$Alt), 1), "deg\n")

# Use a named observatory location
obs2 <- dayup(RA = "05:34:32", Dec = "+22:00:52",
              Date = c(2023, 1, 15), Time = c(20, 0, 0),
              Loc = "AAT", UTCdiff = 11)
cat("Target max altitude from AAT:", round(max(obs2$obs$Alt), 1), "deg\n")
}
