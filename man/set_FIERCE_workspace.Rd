% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FIERCE.R
\name{set_FIERCE_workspace}
\alias{set_FIERCE_workspace}
\title{set_FIERCE_workspace}
\usage{
set_FIERCE_workspace()
}
\value{
The function loads in the current R workspace all the necessary helper functions
}
\description{
Sets the current R workspace for the FIERCE pipeline, with all the helper functions that are necessary to interface R with python. This function must be always launched before using FIERCE
}
\examples{
Simply launch the function with no arguments

set_FIERCE_workspace()

}
