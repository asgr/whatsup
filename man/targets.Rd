\name{targets}
\alias{targets}
\docType{data}
\title{Astronomical Target Catalogue}
\description{
A data frame of astronomical targets (galaxies, star clusters, and
notable stars) used by \code{\link{nowup}}, \code{\link{gettarget}},
and \code{\link{yearup}} for catalogue look-ups.  Includes Messier
objects, NGC/IC objects, and other commonly observed sources.
}
\usage{data(targets)}
\format{
A data frame with one row per target and six variables:
\describe{
  \item{Name}{Character. Object identifier (e.g. \code{"M31"},
    \code{"NGC_4889"}).}
  \item{RA}{Character. Right Ascension in \code{"HH:MM:SS"} format
    (J2000).}
  \item{Dec}{Character. Declination in \code{"sDD:MM:SS"} format
    (J2000).}
  \item{RAdeg}{Numeric. Right Ascension in decimal degrees (J2000).}
  \item{Decdeg}{Numeric. Declination in decimal degrees (J2000).}
  \item{Type}{Character. Object type code: \code{'G'} (galaxy),
    \code{'S'} (star), or \code{'C'} (cluster / nebula).}
}
}
\seealso{
\code{\link{gettarget}}, \code{\link{nowup}}, \code{\link{telescopes}}
}
\examples{
data(targets)

# Number of targets by type
table(targets$Type)

# List all galaxy entries
galaxies <- targets[targets$Type == "G", ]
cat("Number of galaxies:", nrow(galaxies), "\n")
head(galaxies[, c("Name", "RA", "Dec")])

# Find a specific target
targets[targets$Name == "M31", ]
}
