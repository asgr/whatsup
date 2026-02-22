\name{whatsupeuler}
\alias{whatsupeuler}
\title{Euler Angle Coordinate Transformation}
\description{
Performs one of six standard astronomical coordinate transformations
using Euler angle rotation matrices.  This is used internally by
\code{\link{getplanet}} to convert heliocentric ecliptic coordinates to
equatorial (RA/Dec) coordinates, but is available for direct use.
}
\usage{
whatsupeuler(ai, bi, select, fk4 = FALSE, radian = FALSE)
}
\arguments{
  \item{ai}{Numeric vector of input longitudes (or Right Ascensions) in
    degrees (or radians if \code{radian = TRUE}).}
  \item{bi}{Numeric vector of input latitudes (or Declinations) in
    degrees (or radians if \code{radian = TRUE}).}
  \item{select}{Integer from 1 to 6 selecting the transformation:
    \enumerate{
      \item Galactic -> FK4 B1950 equatorial
      \item FK4 B1950 equatorial -> Galactic
      \item FK4 B1950 equatorial -> FK4 B1950 ecliptic
      \item FK4 B1950 ecliptic -> FK4 B1950 equatorial
      \item FK4 B1950 equatorial -> supergalactic
      \item Supergalactic -> FK4 B1950 equatorial
    }
    When \code{fk4 = FALSE} (the default), the J2000 FK5 Euler constants
    are used instead.}
  \item{fk4}{Logical.  If \code{TRUE}, use the FK4 B1950 Euler constants;
    otherwise use FK5 J2000 constants (the default).}
  \item{radian}{Logical.  If \code{TRUE}, input and output angles are in
    radians; otherwise degrees (the default).}
}
\value{
A named list with two elements:
\describe{
  \item{ao}{Numeric vector of output longitudes (or Right Ascensions)
    in the same angular units as the input.}
  \item{bo}{Numeric vector of output latitudes (or Declinations) in
    the same angular units as the input.}
}
}
\seealso{
\code{\link{getplanet}}
}
\examples{
# Convert J2000 ecliptic longitude/latitude to equatorial RA/Dec
# select = 4: ecliptic -> equatorial (FK4 sign convention requires
#              using select=4 with fk4=FALSE for J2000 ecliptic)
ecl_lon <- 120   # degrees
ecl_lat <-  30   # degrees
result <- whatsupeuler(ecl_lon, ecl_lat, select = 4)
cat("RA:", round(result$ao, 3), "deg\n")
cat("Dec:", round(result$bo, 3), "deg\n")

# Round-trip: equatorial -> galactic -> equatorial
ra  <- 187.5
dec <- 12.0
gal <- whatsupeuler(ra, dec, select = 2)  # eq -> galactic (FK4=F -> J2000)
back <- whatsupeuler(gal$ao, gal$bo, select = 1)
cat("Original RA:", ra, "  Recovered RA:", round(back$ao, 6), "\n")
}
