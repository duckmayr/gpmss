% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CovSEiso.R
\name{CovSEiso}
\alias{CovSEiso}
\title{Isotropic Squared Exponential Covariance Function}
\description{
The isometric squared exponential covariance function, also  called the
Gaussian kernel, the squared exponential kernel, or the radial basis
function.
}
\details{
The isometric squared exponential covariance function has two
hyperparameters, \eqn{\sigma_f}, the scale factor, and \eqn{\ell}{ell},
the characteristic length scale. The scale factor governs (on average)
how far from the mean function values can be, while the length scale governs
how quickly the function can change; or, in other words, as the function
output covariance is given as a function of distance in the covariate space,
the length scale governs what "closeness" means.

The covariance between \eqn{f(\mathbf{x}_i)}{f(x_i)}
and \eqn{f(\mathbf{x}_j)}{f(x_j)}
is given by

\deqn{
    k \left( \mathbf{x}_i, \mathbf{x}_j \right)
    =
    \sigma_f^2
    \exp \left(
        \left( \mathbf{x}_i - \mathbf{x}_j \right)^T
        M \left( \mathbf{x}_i - \mathbf{x}_j \right)
    \right),
}{%
    k ( x_i, x_j ) = \sigma_f^2 exp [ ( x_i - x_j )^T M ( x_i - x_j ) ],
}

where \eqn{M} is a matrix whose diagonal entries are
\eqn{1 / \ell^2}{1 / ell^2}.
}
\section{Warning}{

Note that the hyperparameters should be stored on the log scale;
that is, you should supply the log of the scale factor and the log of the
length scale (in that order).
}

\section{Super class}{
\code{\link[gpmss:CovarianceFunction]{gpmss::CovarianceFunction}} -> \code{CovSEiso}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{name}}{A character vector of length one giving the covariance
function's name; "isotropic squared exponential"}

\item{\code{hypers}}{A numeric vector giving the covariance function's
hyperparameters; a vector of length two giving the log of the
scale factor and the log of the length scale, in that order}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-cov}{\code{CovSEiso$cov()}}
\item \href{#method-parameter_derivative}{\code{CovSEiso$parameter_derivative()}}
\item \href{#method-input_derivative}{\code{CovSEiso$input_derivative()}}
\item \href{#method-new}{\code{CovSEiso$new()}}
\item \href{#method-clone}{\code{CovSEiso$clone()}}
}
}
\if{html}{
\out{<details open ><summary>Inherited methods</summary>}
\itemize{
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-cov"></a>}}
\if{latex}{\out{\hypertarget{method-cov}{}}}
\subsection{Method \code{cov()}}{
Compute function covariance
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CovSEiso$cov(X, Z = X, hypers = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{X}}{The first set of input values (should be a numeric matrix)}

\item{\code{Z}}{The second set of input values (should be a numeric matrix);
The default is Z = X.}

\item{\code{hypers}}{A numeric vector giving hyperparameters for the
covariance function. If NULL (the default), the hypers data
member is used.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-parameter_derivative"></a>}}
\if{latex}{\out{\hypertarget{method-parameter_derivative}{}}}
\subsection{Method \code{parameter_derivative()}}{
Compute partial derivatives of covariance function with respect to
its hyperparameters
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CovSEiso$parameter_derivative(X, Z = X, hypers = NULL, param = 1, K = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{X}}{The first set of input values (should be a numeric matrix)}

\item{\code{Z}}{The second set of input values (should be a numeric matrix);
The default is Z = X.}

\item{\code{hypers}}{A numeric vector giving hyperparameters for the
covariance function. If NULL (the default), the hypers data
member is used.}

\item{\code{param}}{An integer vector of length one; which element of
\code{hypers} should the derivative be taken with respect to?
If 1 (the default), the derivative is taken with respect to the
(log of the) scale factor; if 2, it is taken with respect to the
(log of the) length scale.}

\item{\code{K}}{An optional provision of the pre-computed kernel;
this is useful if parameter_derivative() will be called
repeatedly (for the different hypers) without the kernel
itself changing}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-input_derivative"></a>}}
\if{latex}{\out{\hypertarget{method-input_derivative}{}}}
\subsection{Method \code{input_derivative()}}{
Compute partial derivatives of covariance function with respect to
its inputs
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CovSEiso$input_derivative(
  X,
  Z = X,
  hypers = NULL,
  dimension = 1,
  order = 1,
  K = NULL
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{X}}{The first set of input values (should be a numeric matrix)}

\item{\code{Z}}{The second set of input values (should be a numeric matrix);
The default is Z = X.}

\item{\code{hypers}}{A numeric vector giving hyperparameters for the
covariance function. If NULL (the default), the hypers data
member is used.}

\item{\code{dimension}}{an integer vector of length one giving the dimension
of X with respect to which the derivative is being taken; the
default is 1}

\item{\code{order}}{An integer vector of length one indicating whether the
first partial derivative (order = 1) is desired, or the cross
partial (order = 2); the default is 1}

\item{\code{K}}{An optional provision of the pre-computed kernel;
this is useful if parameter_derivative() will be called
repeatedly (for the different hypers) without the kernel
itself changing}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new CovSEiso object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CovSEiso$new(hypers = c(0, 0))}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{hypers}}{A numeric vector giving hyperparameters for the
covariance function; a vector of length two giving the log of the
scale factor and the log of the length scale, in that order}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CovSEiso$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
