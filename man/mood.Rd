% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mood.R
\name{mood}
\alias{mood}
\title{Setting up the parameters of a factorial experiment to search for multi-objective optimal completely randomised design.}
\usage{
mood(
  K,
  Levels,
  Nruns,
  criterion.choice = c("GD", "GL", "GDP", "GLP", "MSE.D", "MSE.L", "MSE.P"),
  kappa = list(),
  control = list(),
  prob = list(),
  model_terms = list(primary.model = "first_order")
)
}
\arguments{
\item{K}{Number of factors.}

\item{Levels}{Either (a) a common number of levels for each factor or (b) a list of length K of the vectors containing levels of each factor.}

\item{Nruns}{Number of runs of the experiment.}

\item{criterion.choice}{Compound criterion to be used for the optimal design search or evaluation.
Possible values are:
\itemize{
\item \code{GL}, \code{GD} for Generalised D- and L-optimality (Goos et al., 2005)
\item \code{GDP} and \code{GLP} for Generalised DP- and LP-optimality (Trinca and Gilmour, 2012)
\item \code{MSE.D}, \code{MSE.L} and \code{MSE.P} for compound criteria with MSE-component: determinant-based, trace-based and determinant-based but with point estimates for the MSE(D)-component
}}

\item{kappa}{List specifying the weights used in the compound criterion. Each named entry must be between 0 and 1.
\itemize{
\item \code{kappa.Ds} Weight of the Ds-criterion (default = 1 if \code{criterion.choice = GD})
\item \code{kappa.DP} Weight of the DP-criterion (default = 1 if \code{criterion.choice = GDP})
\item \code{kappa.Ls} Weight of the Ls-criterion (default = 1 if \code{criterion.choice = GL})
\item \code{kappa.LP} Weight of the LP-criterion (default = 1 if \code{criterion.choice = GLP})
\item \code{kappa.LoF} Weight of the lack-of-fit criterion
\item \code{kappa.bias} Weight of the bias criterion
\item \code{kappa.mse} Weight of the MSE criterion (default = 1 if 'criterion.choice = MSE.*')
}}

\item{control}{Named list specifying control parameters for the design search.
\itemize{
\item \code{Nstarts} The number of randomly generated start designs of the search algorithm (default = \code{10}).
\item \code{Cubic} Indicator of whether the experimental region is cubic (\code{TRUE}, default) or spherical (\code{FALSE}).
\item \code{tau2} The variance scaling parameter for the prior distribution of the potential terms (default = 1).
\item \code{Biter} Number of samples for evaluating the MSE determinant-based component criterion (default = \code{50}).
\item \code{MC} Indicator of whether to apply a multiple comparison (Bonferroni) correction for trace-based criteria (\code{TRUE}, default) or not (\code{FALSE}).
\item \code{orth} Indicator of whether to orthonormalise the potential and primary terms (\verb{TRUE', default) or not (}FALSE').
}}

\item{prob}{Named list specifying confidence levels for DP- (\code{prob.DP}), LP- (\code{prob.LP}) and Lack-of-fit (\code{prob.LoF}) criteria (pre-Bonferroni correction).
All default to 0.95.}

\item{model_terms}{A list specifying the primary (fitted) and potential (biased) models with the following named elements (see Details).
\itemize{
\item \code{primary.model} The order of the fitted polynomial model. Alternatively polynomial terms can be directly specified through the \code{primary.terms} parameter.
\item \code{potential.model} The order of the potential/extra polynomial terms. Alternatively can be specified through the \code{potential.terms} parameter.
\item \code{primary.terms} Alternative specification of the primary model via character vector of the names of the primary terms.
\item \code{potential.terms} Alternative specification of the potential model via character vector of the names of the potential terms.
}}
}
\value{
List of parameters of the experiment, compound criterion of choice, and primary and potential model terms.
\itemize{
\item \code{K} Number of factors.
\item \code{Klev} Number of levels of each factor, if all factors have the same number of levels.
\item \code{Levels} List of length K of the vectors containing values of the factors.
\item \code{Nruns} Number of runs of the experiment.
\item \code{criterion.choice} Compound criterion to be used for the optimal design search or evaluation.
\item \code{Nstarts} The number of randomly generated start designs of the search algorithm.
\item \code{Biter} Number of samples for evaluating the MSE determinant-based component criterion.
\item \code{tau2} The variance scaling parameter for the prior distribution of the potential terms.
\item \code{tau}  The square root of \code{tau2}
\item \code{Cubic} Whether the experimental region is cubic (\code{TRUE}) or spherical (\code{FALSE}).
\item \code{MC} Indicator of the multiple comparison (Bonferroni) correction for trace-based criteria.
\item \verb{prob.DP, prob.LP, prob.LoF, prob.LoFL} Confidence levels for the DP-, LP-, lack of fit determinant- and trace-based criteria.
\item \verb{alpha.DP, alpha.LP, alpha.LoF, alpha.LoFL} Significance levels for the DP-, LP-, lack of fit determinant- and trace-based criteria.
\item \code{orth} Whether the candidate sets are orthonormalised (\code{TRUE}) or not (\code{FALSE}).
\item \code{Z0} Z0 matrix.
\item \code{W} Weight matrix for Ls criterion.
\item \code{primary.terms} Fitted (primary) model terms.
\item \code{potential.terms} Potential terms.
\item \code{P} The number of terms in the fitted model (including intercept).
\item \code{Q} The number of potential terms.
\item \verb{kappa.Ds, kappa.DP, kappa.Ls, kappa.LP,  kappa.LoF, kappa.bias, kappa.mse} Compound criterion weights.
\item \code{warning.msg} Warning messages.
}
}
\description{
Creates an object containing the parameters of the experiment, compound optimality criterion with the
weights and parameters of the search.
}
\details{
The function provides different ways of specifying the levels of the factors and the models. Although some default options are provided
for, e.g., \code{criterion.choice} and \verb{kappa.*} values, specification of these input parameters should be carefully chosen to reflect the aims of the experiment and available prior information.

Specifying the factors and levels

If all \code{K} factors have the same number of levels, \code{Levels} parameter is used to input that number.
Otherwise, \code{Levels} is set to be a list of vectors containing the values of the factors, e.g.
\code{list(1:3, 1:2, 1:4)} for 3 factors with equally spaced \eqn{3, 2} and \eqn{4} levels respectively.

Specifying the fitted model and the potential terms

There are two ways to describe the primary and potential sets of model terms via the \code{model_terms} list.
Named elements \code{primary.model} and \code{potential.model} can be used to specify the fitted model and the potential terms via a string form.
They are used to generate the sets of \code{primary.terms} and \code{potential.terms} which alternatively can be input directly.
Possible values of \code{primary.model} are:
\itemize{
\item \code{main_effects} -- main effects for all the factors (default for all criteria)
\item \code{first_order} -- main effects and linear interactions
\item \code{second_order} -- full second order polynomial
\item \code{third_order} -- full second order polynomial model and all interactions of degree 3
\item \code{cubic} -- third order polynomial model with cubic terms
}
The intercept is always included as a primary term.

Possible elements of the \code{potential.model} vector argument:
\itemize{
\item \code{linear_interactions} -- linear interactions among the factors (default for MSE criteria)
\item \code{quadratic_terms} -- quadratic terms for all the factors
\item \code{third_order_terms} --  all interactions of degree 3: linear-by-linear-by-linear and quadratic-by-linear terms
\item \code{cubic_terms} -- cubic terms for all the factors
\item \code{fourth_order_terms} -- all interactions of degree 4, similar to \code{third_order_terms}
}
\code{primary.terms} and \code{potential.terms} arguments are used to specify the fitted model and the potential terms explicitly -- up to the total power of 4.
\itemize{
\item Single factor powers,  are coded as a string starting with with "x" and followed by
the index of the factor and the power: \code{"x32"}.
For example, \eqn{x_3^2} is coded as \code{"x32"}; \code{"x22"} stands for \eqn{x_2^2}, and \code{"x4"} stands for the linear term \eqn{x_4}.
\item Interaction of factors' powers are coded by merging the individual factors'
powers, ordered by the factors' indexes. For example, \eqn{x_2^2\times x_1} is coded as \code{"x1x22"},
\eqn{x_3x_12x_4} -- as \code{"x12x3x4"}.
}
}
\examples{

example1 <- mood(K = 5, Levels = 3, Nruns = 40, criterion.choice = "GDP", 
kappa = list(kappa.Ds = 1./3, kappa.DP = 1./3, kappa.LoF = 1./3), 
control = list(Nstarts = 50, tau2 = 0.1),
model_terms = list(primary.model = "second_order", potential.terms = c("x12x2", "x22x3", "x32x4", "x42x5")))
example1

example2 <- mood(K = 3, Nruns = 12, Levels = list(1:3, 1:2, 1:2), criterion.choice = "MSE.L",
kappa = list(kappa.LP = 1./2, kappa.LoF = 1./4, kappa.mse = 1./4),
control = list(Nstarts = 50, tau2 = 1), 
model_terms = list(primary.terms = "first_order",
potential.terms = c("x12", "x12x2", "x12x3")))
example2

}
