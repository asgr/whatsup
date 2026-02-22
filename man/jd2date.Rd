\name{jd2date}
\alias{jd2date}
\title{Convert a Julian Date to a Calendar Date}
\description{
Converts one or more Julian Dates to Gregorian calendar date components
(year, month, day, and decimal hour).  Vectorised.
}
\usage{
jd2date(JD = 2440000)
}
\arguments{
  \item{JD}{Numeric vector of Julian Dates to convert.}
}
\value{
A named list with four numeric vector components:
\describe{
  \item{year}{Four-digit year.}
  \item{mon}{Month (1 = January, 12 = December).}
  \item{mday}{Day of the month.}
  \item{hour}{Decimal hour (0–24) of the day in UT.}
}
}
\seealso{
\code{\link{date2jd}}, \code{\link{dayup}}
}
\examples{
# Convert a single Julian Date
jd2date(2440000)

# Convert multiple Julian Dates
jd2date(c(2451545.0, 2455927.5, 2460000.0))

# Round-trip check: date -> JD -> date
jd <- date2jd(2023, 3, 21, hour = 9.5)
back <- jd2date(jd)
cat("Year:", back$year, "Month:", back$mon, "Day:", back$mday,
    "Hour:", round(back$hour, 4), "\n")
}
