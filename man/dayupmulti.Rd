\name{dayupmulti}
\alias{dayupmulti}
\title{Compute Daily Observability for Multiple Targets}
\description{
Runs \code{\link{dayup}} for each row in a targets matrix and bundles the
results into a single object of class \code{dayupmulti} for combined plotting
with \code{\link{plot.dayupmulti}} / \code{\link{plotdayupmulti}}.
}
\usage{
dayupmulti(Targets = cbind(
             Names = c('G02', 'G09', 'G12', 'G15', 'G23'),
             RA    = c("02:00:00", "09:00:00", "12:00:00", "14:30:00", "23:00:00"),
             Dec   = c("-05:00:00", "01:00:00", "-01:00:00", "01:00:00", "-32:30:00")),
           ...)
}
\arguments{
  \item{Targets}{A matrix or data frame with at least three columns named (case
    insensitively) \code{Names}, \code{RA}, and \code{Dec}.  Each row describes
    one target: \code{RA} in \code{"HH:MM:SS"} format and \code{Dec} in
    \code{"sDD:MM:SS"} format.}
  \item{\dots}{Additional arguments passed to \code{\link{dayup}} (e.g.
    \code{Date}, \code{Time}, \code{Lon}, \code{Lat}, \code{Loc},
    \code{UTCdiff}).}
}
\value{
A list of class \code{dayupmulti}.  The first \emph{N} elements (one per
target) are each a \code{dayup} object as returned by \code{\link{dayup}}.
The final element, \code{names}, is a character vector of target names
taken from the \code{Names} column of \code{Targets}.
}
\seealso{
\code{\link{dayup}}, \code{\link{plotdayupmulti}}, \code{\link{yearup}}
}
\examples{
# Define a small set of targets and compute their observability
targets_mat <- cbind(
  Names = c("FieldA", "FieldB", "FieldC"),
  RA    = c("02:00:00", "09:30:00", "22:00:00"),
  Dec   = c("-05:00:00", "-10:00:00", "-30:00:00")
)

obs_multi <- dayupmulti(Targets = targets_mat,
                        Date = c(2023, 7, 1),
                        Time = c(18, 0, 0),
                        Lon = 115.8178, Lat = -31.97333,
                        UTCdiff = 8)

# Access per-target results
cat("Max altitude of FieldA:",
    round(max(obs_multi[[1]]$obs$Alt), 1), "deg\n")
cat("Max altitude of FieldB:",
    round(max(obs_multi[[2]]$obs$Alt), 1), "deg\n")
}
