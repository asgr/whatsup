\name{telescopes}
\alias{telescopes}
\docType{data}
\title{Observatory and Telescope Location Data}
\description{
A data frame of observatory locations and telescope properties used by
\code{\link{gettelescope}}, \code{\link{dayup}}, \code{\link{nowup}},
\code{\link{nowwhere}}, and \code{\link{yearup}} when specifying a named
location via the \code{Loc} argument.
}
\usage{data(telescopes)}
\format{
A data frame with one row per facility and six variables:
\describe{
  \item{Name}{Character. Telescope or facility name (used as the key for
    look-up; matching is case insensitive).}
  \item{Lat}{Numeric. Observer latitude in decimal degrees (positive north).}
  \item{Lon}{Numeric. Observer longitude in decimal degrees (positive east).}
  \item{Height}{Numeric. Altitude above sea level in metres.}
  \item{Diam}{Numeric. Primary mirror (or dish) diameter in metres.}
  \item{N}{Integer. Number of individual telescopes / antennas at the
    facility.}
}
}
\seealso{
\code{\link{gettelescope}}, \code{\link{targets}}
}
\examples{
data(telescopes)
head(telescopes)

# Find all facilities in the Southern Hemisphere
southern <- telescopes[telescopes$Lat < 0, ]
cat("Number of southern facilities:", nrow(southern), "\n")

# Look up a specific telescope
gettelescope("AAT")
}
