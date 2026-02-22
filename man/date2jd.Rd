\name{date2jd}
\alias{date2jd}
\title{Convert a Calendar Date to Julian Date}
\description{
Converts a Gregorian calendar date (and optional time) to the
corresponding Julian Date (JD).  Vectorised: all arguments may be
vectors of the same length.
}
\usage{
date2jd(year = 1968, mon = 5, mday = 23, hour = 12)
}
\arguments{
  \item{year}{Integer (or numeric vector) giving the four-digit year.}
  \item{mon}{Integer (or numeric vector) giving the month (1 = January,
    12 = December).}
  \item{mday}{Integer (or numeric vector) giving the day of the month.}
  \item{hour}{Numeric (or numeric vector) giving the time as a decimal
    hour (0–24).  Defaults to noon (12).}
}
\value{
Numeric vector of Julian Dates.
}
\note{
The Julian Date starts at noon UT on 1 January 4713 BC (Julian calendar).
A Julian Date of 2440000.5 corresponds to 1968-05-23 00:00:00 UT.
}
\seealso{
\code{\link{jd2date}}, \code{\link{dayup}}
}
\examples{
# Julian Date for the Unix epoch (1970-01-01 00:00:00 UTC)
date2jd(1970, 1, 1, hour = 0)

# Multiple dates at once
date2jd(year = c(2000, 2010, 2020),
        mon  = c(1, 6, 12),
        mday = c(1, 15, 31),
        hour = 12)

# Round-trip check
jd <- date2jd(2023, 6, 15, hour = 18.5)
back <- jd2date(jd)
cat("Year:", back$year, "Month:", back$mon, "Day:", back$mday,
    "Hour:", round(back$hour, 4), "\n")
}
