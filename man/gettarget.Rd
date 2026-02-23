\name{gettarget}
\alias{gettarget}
\title{Look Up a Target's Sky Coordinates}
\description{
Returns the Right Ascension and Declination of a named target.  The
function first checks the built-in \code{\link{targets}} catalogue,
then handles planets and solar-system bodies internally, and finally
falls back to the online Sesame name resolver (via
\code{\link{nameresolve}}) for unrecognised names.
}
\usage{
gettarget(name, JD = 2440000)
}
\arguments{
  \item{name}{Character string giving the target name.  Valid values
    include any entry in the \code{\link{targets}} dataset, the planet
    names \code{'mercury'}, \code{'venus'}, \code{'earth'},
    \code{'mars'}, \code{'jupiter'}, \code{'saturn'}, \code{'uranus'},
    \code{'neptune'}, \code{'pluto'}, or the special names \code{'moon'}
    and \code{'sun'}.  Case insensitive.  For unrecognised names the
    Sesame name resolver is queried (requires internet access and the
    \pkg{XML} and \pkg{httr} packages).}
  \item{JD}{Julian Date used when evaluating the position of moving
    targets (planets, Moon, Sun).}
}
\value{
A data frame with columns:
\describe{
  \item{Name}{Target name.}
  \item{RA}{Right Ascension in \code{"HH:MM:SS"} format.}
  \item{Dec}{Declination in \code{"sDD:MM:SS"} format.}
  \item{RAdeg}{Right Ascension in decimal degrees.}
  \item{Decdeg}{Declination in decimal degrees.}
  \item{Type}{Source type code (e.g. \code{'G'} for galaxy, \code{'S'}
    for star, \code{'P'} for planet, \code{'M'} for Moon, \code{'R'} for
    resolved name).}
}
}
\seealso{
\code{\link{targets}}, \code{\link{getplanet}}, \code{\link{nameresolve}},
\code{\link{dayup}}
}
\examples{
# Look up a catalogue source
gettarget("G09")

# Look up a planet (position depends on the Julian Date supplied)
gettarget("mars", JD = 2460000)

# Look up the Moon
gettarget("moon", JD = 2460000)

# Look up the Sun
gettarget("sun", JD = 2460000)
}
