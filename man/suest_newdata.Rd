\name{suest_newdata}
\alias{suest_newdata}
\title{Stack the Estimation Samples from a SUEST Object}
\usage{
suest_newdata(object)
}
\arguments{
\item{object}{A \code{suest_model} returned by \code{suest()}.}
}
\value{
A data frame containing the component model frames and internal model-routing
columns.
}
\description{
Creates model-specific new data using each component model's own estimation
sample. This is useful when effects should be averaged separately within
different samples.
}
