\name{airmass}
\alias{airmass}
\title{Compute Airmass from Altitude}
\description{
Computes the airmass at a given altitude angle using the Rozenberg
approximation (1966).  Returns \code{NA} for targets below the horizon.
}
\usage{
airmass(Alt = 30)
}
\arguments{
  \item{Alt}{Numeric vector of altitude angles in degrees (can contain
    negative values for targets below the horizon).}
}
\value{
Numeric vector of airmass values.  Values for \code{Alt < 0} are
returned as \code{NA}.  At the zenith (\code{Alt = 90}) the airmass is
1; at \code{Alt = 30} it is approximately 2.
}
\seealso{
\code{\link{dayup}}, \code{\link{nowwhere}}
}
\examples{
# Airmass at zenith
airmass(90)

# Airmass at typical low-elevation angles
airmass(c(10, 20, 30, 45, 60, 90))

# Below the horizon returns NA
airmass(-5)

# Convert an altitude track to airmass
alt_deg <- seq(0, 90, by = 10)
data.frame(Alt = alt_deg, AirMass = airmass(alt_deg))
}
