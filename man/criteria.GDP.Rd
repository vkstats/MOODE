% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/criteria_GDP.R
\name{criteria.GDP}
\alias{criteria.GDP}
\title{Calculates the values of the Generalised DPs-criterion and its components}
\usage{
criteria.GDP(X1, X2, search.object, eps = 10^-23)
}
\arguments{
\item{X1}{The primary model matrix, with the first column containing the labels of treatments, and the second -- the intercept term.}

\item{X2}{The matrix of potential terms, with the first column containing the labels of treatments.}

\item{search.object}{Object of class \code{\link[=mood]{mood()}} specifying experiment parameters.}

\item{eps}{Computational tolerance, the default value is 10^-23}
}
\value{
A list of values: indicator of whether the evaluation was successful ("eval"), Ds-criterion value -- intercept excluded ("Ds"),
DPs-criterion value -- intercept excluded ("DPs"), Lack-of-fit(DP) criterion value ("LoF"), the bias component value ("bias"),
the number of pure error degrees of freedom ("df") and the value of the compound criterion ("compound").
}
\description{
This function evaluates the Generalised DPs-criterion for given primary and potential model matrices.
Components: Ds-, DPs-, LoF(DP)- and Bias(D)-optimality.
}
\examples{
# Experiment: one 5-level factor, primary model -- full quadratic, X^3 and X^4 potential terms.
ex.mood <- mood(K = 1, Levels = 5, Nruns = 8, criterion.choice = "GDP", 
               kappa = list(kappa.Ds = .25, kappa.LoF = .25, kappa.bias = .25, kappa.DP = .25), 
               model_terms = list(primary.model = "second_order", potential.terms = "x14"))
# Generating candidate sets: primary and full orthonormalised ones
K <- ex.mood$K
Levels <- ex.mood$Levels 
cand.not.orth <- candidate_set_full(candidate_trt_set(Levels, K), K)
cand.full.orth <- candidate_set_orth(cand.not.orth, ex.mood$primary.terms, ex.mood$potential.terms)
# Choosing a design
index <- c(rep(1,2),3,rep(4,2),rep(5,3))
X.primary <- cand.full.orth[index, c(1, match(ex.mood$primary.terms, colnames(cand.full.orth)))]
X.potential <- cand.full.orth[index, (c(1, match(ex.mood$potential.terms, colnames(cand.full.orth))))]
# Evaluating a compound GDP-criterion
criteria.GDP(X1 = X.primary, X2 = X.potential, ex.mood)

# Output: eval = 1, Ds = 0.6884783, DP = 4.4538023, LoF = 3.895182, bias = 1.03807, df = 4, compound = 2.465318

}
