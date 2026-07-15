\name{suest}
\alias{suest}
\alias{coef.suest_model}
\alias{nobs.suest_model}
\alias{print.suest_model}
\alias{vcov.suest_model}
\title{Combine Two Fitted Models with Seemingly Unrelated Estimation}
\usage{
suest(model1, model2, model_names = NULL)

\method{coef}{suest_model}(object, ...)

\method{vcov}{suest_model}(object, ...)

\method{nobs}{suest_model}(object, ...)

\method{print}{suest_model}(x, ...)
}
\arguments{
\item{model1, model2}{Two supported fitted model objects.}

\item{model_names}{Optional character vector containing two display names.}

\item{object, x}{A \code{suest_model} object.}

\item{...}{Additional arguments.}
}
\value{
An object of class \code{"suest_model"} containing the component models,
a joint coefficient vector, a joint model-robust covariance matrix, and
sample-alignment information.
}
\description{
Combines two separately fitted models and constructs a joint model-robust
covariance matrix from their observation-level score contributions.
}
\details{
Exactly two models are supported. Models may use identical, partially
overlapping, or disjoint samples. Supported model classes include linear,
binary logit and probit, Poisson, negative binomial, ordered logit and probit,
and multinomial logit models.

The returned object works with the \pkg{marginaleffects} package.
}
\examples{
dat <- mtcars
dat$am <- factor(dat$am)
model1 <- glm(am ~ wt, family = binomial(), data = dat)
model2 <- glm(am ~ wt + hp, family = binomial(), data = dat)
fit <- suest(model1, model2, model_names = c("Base", "Adjusted"))
fit
}
