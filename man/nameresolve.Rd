\name{nameresolve}
\alias{nameresolve}
\title{Resolve an Astronomical Object Name to Sky Coordinates}
\description{
Queries the CDS Sesame name resolver service to obtain the J2000 Right
Ascension and Declination of one or more astronomical objects by name.
Requires internet access and the \pkg{XML} and \pkg{httr} packages.
}
\usage{
nameresolve(name = "M31")
}
\arguments{
  \item{name}{Character string or character vector of object name(s) to
    resolve (e.g. \code{"M31"}, \code{"NGC_4889"}, \code{"Andromeda_Galaxy"}).
    Spaces in names are automatically replaced with underscores before
    the query is submitted.}
}
\value{
A data frame with one row per name and columns:
\describe{
  \item{Name}{The requested name.}
  \item{RA}{Right Ascension in \code{"HH:MM:SS"} format, or \code{NA}
    if the name could not be resolved.}
  \item{Dec}{Declination in \code{"sDD:MM:SS"} format, or \code{NA} if
    unresolved.}
  \item{RAdeg}{Right Ascension in decimal degrees, or \code{NA}.}
  \item{Decdeg}{Declination in decimal degrees, or \code{NA}.}
  \item{Type}{Always \code{'R'} (resolved).}
}
}
\note{
This function requires internet access to contact
\url{http://cdsweb.u-strasbg.fr/cgi-bin/nph-sesame}.
The \pkg{XML} and \pkg{httr} packages must be installed.  If the
service is unavailable or the name is ambiguous/not found, \code{NA}
values are returned for the coordinate columns.
}
\seealso{
\code{\link{gettarget}}, \code{\link{dayup}}
}
\examples{
\dontrun{
# Resolve a single Messier object
nameresolve("M31")

# Resolve multiple objects
nameresolve(c("M31", "NGC_4889", "Crab_Nebula"))
}
}
