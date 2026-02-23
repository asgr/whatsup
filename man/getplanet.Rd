\name{getplanet}
\alias{getplanet}
\title{Get the Current Sky Position of a Solar-System Planet}
\description{
Computes the J2000 Right Ascension and Declination of a named planet at
a given Julian Date using heliocentric orbital elements from the
\pkg{astrolibR} package.
}
\usage{
getplanet(name = 'mars', JD = 2440000)
}
\arguments{
  \item{name}{Character string giving the planet name (case
    insensitive).  Must be one of \code{'mercury'}, \code{'venus'},
    \code{'earth'}, \code{'mars'}, \code{'jupiter'}, \code{'saturn'},
    \code{'uranus'}, \code{'neptune'}, or \code{'pluto'}.}
  \item{JD}{Julian Date at which to evaluate the planetary position.}
}
\value{
A data frame with one row and columns:
\describe{
  \item{Name}{The planet name.}
  \item{RA}{Right Ascension in \code{"HH:MM:SS"} format.}
  \item{Dec}{Declination in \code{"sDD:MM:SS"} format.}
  \item{RAdeg}{Right Ascension in decimal degrees.}
  \item{Decdeg}{Declination in decimal degrees.}
  \item{Type}{Always \code{'P'}.}
}
}
\seealso{
\code{\link{gettarget}}, \code{\link{dayup}}
}
\examples{
# Position of Mars at a fixed Julian Date
getplanet("mars", JD = 2460000)

# Position of Jupiter
getplanet("jupiter", JD = 2460000)

# All outer planets
planets <- c("mars", "jupiter", "saturn", "uranus", "neptune")
do.call(rbind, lapply(planets, getplanet, JD = 2460000))
}
