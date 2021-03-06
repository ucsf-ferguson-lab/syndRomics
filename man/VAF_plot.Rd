% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/user_functions.R
\name{VAF_plot}
\alias{VAF_plot}
\title{Variance accounted for plot}
\usage{
VAF_plot(
  pca,
  pca_data,
  ndim = 1:5,
  resample_ci = NULL,
  style = "line",
  colors = c("steelblue", "orange")
)
}
\arguments{
\item{pca}{Object of class \emph{prcomp}, \emph{princals}, or \emph{data.frame}. If object is a \emph{prcomp} or \emph{princals} object, \emph{pca_data} is required, and
the loadings will be extracted. If object is a data.frame object, the dataframe needs to be formatted as: first column named \emph{Variables} and all other columns corresponding to a PC.
One row per variable. The values are the loadings.}

\item{pca_data}{Data passed to the \emph{prcomp} or \emph{princals} function.}

\item{ndim}{Numeric. Number of PCs to plot}

\item{resample_ci}{dataframe. Dataframe containing the columns "original", "mean","ci_low" and "ci_high" containing the center, the lower bound and the upper bound of the confident intervals to plot. Each row contain the values of 1 PC in order (PC1 first row, PC2 second row, etc).
This can be obtained by the permut_pc_test "results" element with statistic="VAF" or computed independently.}

\item{style}{Character. There are two styles of VAF plots "line" by default or "reduced".}

\item{colors}{Character vector of length 2. Vector with the character name or
hexadecimal number (e.g. "#FF0000") of two colors, the lower color and the higher color
for the gradient used in the plot. Hexadecimal number can be obtained using \emph{rgb} for example.}
}
\value{
Returns a \emph{ggplot2} object.
}
\description{
Plot of VAF from a PCA solution given the first ndim PCs.
}
\examples{
data(mtcars)
pca_mtcars<-prcomp(mtcars, center = TRUE, scale. = TRUE)

VAF_plot(pca = pca_mtcars, pca_data = mtcars, ndim = 1:7)

}
\author{
Abel Torres Espin
}
