# whatsup (R package)

<!-- badges: start -->
![R-CMD-check](https://github.com/asgr/whatsup/workflows/R-CMD-check/badge.svg)
<!-- badges: end -->

## Synopsis

Package containing tools for producing observability plots for arbitrary locations, dates and times. Useful for planning astronomical observations by showing which targets are visible from a given site on a given night.

## Installation

### Getting whatsup

The **whatsup** package can be installed directly from GitHub:

```R
install.packages('remotes')
remotes::install_github("asgr/whatsup")
library(whatsup)
```

A few Mac people seem to have issues with the above due to the backend used to download files. A work around seems to be to either use devtools (which I do not use as the default since it has a few more dependencies, and is tricky to install on HPCs):

```R
install.packages('devtools')
devtools::install_github("asgr/whatsup")
library(whatsup)
```

Or try the following:

```R
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
remotes::install_github("asgr/whatsup")
```

I also have these options set by default in my .Rprofile, which seems to help with some of the remote install issues some people face:

```R
options(download.file.method = "libcurl")
options(repos="http://cran.rstudio.com/")
options(rpubs.upload.method = "internal")
```

If all of these do not work then the nuclear option is to download (or clone) the GitHub repo, cd to where the tar.gz file is and run in the **console** (or **Terminal** on Mac):

```console
R CMD install whatsup_X.Y.Z.tar.gz
```

where X, Y and Z should be set as appropriate for the version downloaded (check the name of the file basically).

If none of the above works then you should consider burning your computer in sacrifice to the IO Gods. Then buy a newer *better* computer, and try all the above steps again.

Failing all of the above, please email me for help (or perhaps raise an Issue here, if it really does not seem like a local issue).

#### Package Dependencies

The above should also install the required packages. If you have trouble with this you can try installing the required packages manually first and then retry the installation for **whatsup**:

```R
install.packages('remotes')
remotes::install_github("asgr/astrolibR")
remotes::install_github("asgr/moonsun")
remotes::install_github("asgr/magicaxis")
remotes::install_github("asgr/celestial")
remotes::install_github("asgr/whatsup")
```

Assuming you have installed all of the packages that you need/want, you should now be able to load **whatsup** within **R** with the usual:

```R
library(whatsup)
```

## Code Example

The standard usage is to call `whatsup` to produce an observability plot for a target from a given location on a given date:

```R
library(whatsup)

# Show observability of the Centaurus A galaxy from Perth tonight
whatsup(RA="13:25:27", Dec="-43:01:09", Target="Cen A")

# Use a named telescope location
whatsup(RA="13:25:27", Dec="-43:01:09", Target="Cen A", Loc="AAT")

# Check what is up right now
whatsup(Target="moon")
```

## Motivation

A package for planning astronomical observations by generating observability plots that show altitude and related information for targets over the course of a night, from any location on Earth. It accounts for the position of the Moon and Sun, and supports named telescope sites as well as arbitrary coordinates.

## Contributors

Aaron Robotham

## License

LGPL-3
