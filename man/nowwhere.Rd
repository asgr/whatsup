\name{nowwhere}
\alias{nowwhere}
\title{Current Sky Position of One or More Targets}
\description{
Computes the current (or specified) altitude, azimuth, hour angle and
airmass for one or more targets, as well as the Sun and Moon, from a
given observer location.
}
\usage{
nowwhere(RA = "12:30:16", Dec = "-30:13:15", Target = 'user',
         Date = 'get', Time = 'get',
         Lon = 115.8178, Lat = -31.97333, Loc = 'user',
         UTCdiff = 'get', Altitude = 10, Pressure = 1000,
         Temp = 20, Name = 'Auto')
}
\arguments{
  \item{RA}{Character string or vector of Right Ascensions in
    \code{"HH:MM:SS"} format.  Ignored when \code{Target} is not
    \code{'user'}.}
  \item{Dec}{Character string or vector of Declinations in
    \code{"sDD:MM:SS"} format.  Ignored when \code{Target} is not
    \code{'user'}.}
  \item{Target}{Name of a known target or solar-system body (e.g.
    \code{'mars'}, \code{'moon'}, \code{'sun'}, or any name in the
    \code{targets} dataset).  Use \code{'user'} (the default) to supply
    coordinates via \code{RA} and \code{Dec}.}
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
  \item{Name}{Character vector of target names used for labelling the
    output table.  Use \code{'Auto'} (the default) to generate names
    automatically as \code{'Target-1'}, \code{'Target-2'}, etc.}
}
\value{
A named list with components:
\describe{
  \item{obs}{Data frame with one row per target containing columns
    \code{Name}, \code{RA}, \code{Dec}, \code{Alt}, \code{AirMass},
    \code{Az}, \code{HA}, \code{MoonSep}, and \code{Sunsep}.}
  \item{moonsun}{Data frame with two rows (Sun, Moon) containing columns
    \code{Name}, \code{RA}, \code{Dec}, \code{Alt}, \code{AirMass},
    \code{Az}, and \code{HA}.}
  \item{JD}{Julian Date of the observation.}
  \item{RJD}{Reduced Julian Date (JD - 2400000).}
  \item{LST}{Local Sidereal Time in hours.}
  \item{LT.year, LT.mon, LT.mday, LT.hour}{Components of the local
    date and time.}
  \item{UTCdiff}{UTC offset used (hours).}
  \item{Lon}{Observer longitude used (decimal degrees).}
  \item{Lat}{Observer latitude used (decimal degrees).}
  \item{Altitude}{Observatory altitude used (metres).}
  \item{Pressure}{Pressure used (millibars).}
  \item{Temp}{Temperature used (degrees Celsius).}
}
}
\seealso{
\code{\link{dayup}}, \code{\link{nowup}}, \code{\link{gettarget}}
}
\examples{
# Where is a user-supplied target right now (fixed date for reproducibility)?
pos <- nowwhere(RA = "05:34:32", Dec = "+22:00:52",
                Date = c(2023, 1, 15), Time = c(22, 0, 0),
                Lon = 115.8178, Lat = -31.97333,
                UTCdiff = 8, Name = "CrabNebula")
print(pos$obs)
cat("Moon altitude:", round(pos$moonsun$Alt[2], 1), "deg\n")
cat("Sun altitude:", round(pos$moonsun$Alt[1], 1), "deg\n")

# Query multiple targets at once
pos2 <- nowwhere(
  RA  = c("12:30:16", "02:42:40"),
  Dec = c("-30:13:15", "-00:00:47"),
  Date = c(2023, 6, 1), Time = c(20, 0, 0),
  UTCdiff = 8,
  Name = c("TargetA", "TargetB")
)
print(pos2$obs[, c("Name", "Alt", "Az", "AirMass")])
}
