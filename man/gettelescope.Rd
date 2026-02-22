\name{gettelescope}
\alias{gettelescope}
\title{Look Up Observatory / Telescope Coordinates}
\description{
Returns the latitude, longitude, altitude, primary mirror diameter and
number of telescopes for one or more named observatories from the
built-in \code{\link{telescopes}} dataset.
}
\usage{
gettelescope(name)
}
\arguments{
  \item{name}{Character string (or vector of strings) giving the
    telescope or observatory name(s) to look up.  Matching is case
    insensitive.  See \code{\link{telescopes}} for the full list of
    valid names.}
}
\value{
A data frame (one or more rows) from the \code{\link{telescopes}}
dataset with columns \code{Name}, \code{Lat}, \code{Lon}, \code{Height},
\code{Diam}, and \code{N}.
}
\seealso{
\code{\link{telescopes}}, \code{\link{dayup}}, \code{\link{nowup}}
}
\examples{
# Look up the Anglo-Australian Telescope
gettelescope("AAT")

# Look up multiple telescopes
gettelescope(c("AAT", "ALMA", "VLT"))

# Case is insensitive
gettelescope("aat")
}
