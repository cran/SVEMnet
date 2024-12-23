% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/glmnet_with_cv.R
\name{glmnet_with_cv}
\alias{glmnet_with_cv}
\title{Fit a glmnet Model with Cross-Validation}
\usage{
glmnet_with_cv(
  formula,
  data,
  glmnet_alpha = c(0, 0.5, 1),
  standardize = TRUE,
  nfolds = 10,
  ...
)
}
\arguments{
\item{formula}{A formula specifying the model to be fitted.}

\item{data}{A data frame containing the variables in the model.}

\item{glmnet_alpha}{Elastic Net mixing parameter(s) (default is \code{c(0, 0.5, 1)}).
If multiple values are provided, \code{cv.glmnet} is run for each \code{alpha}, and the model with the lowest cross-validation
error is selected.}

\item{standardize}{Logical flag passed to \code{\link[glmnet]{glmnet}}. If \code{TRUE} (default), each variable is standardized
before model fitting.}

\item{nfolds}{Number of cross-validation folds (default is \code{10}).}

\item{...}{Additional arguments passed to \code{\link[glmnet]{cv.glmnet}}.}
}
\value{
A list containing:
\itemize{
  \item \code{parms}: Coefficients from the selected \code{cv.glmnet} model at \code{lambda.min}.
  \item \code{debias_fit}: A linear model of the form \code{y ~ y_pred} used for debiasing (if applicable).
  \item \code{glmnet_alpha}: The vector of \code{alpha} values considered.
  \item \code{best_alpha}: The selected \code{alpha} value that gave the best cross-validation result.
  \item \code{best_lambda}: The \code{lambda} value chosen by cross-validation at the selected \code{alpha}.
  \item \code{actual_y}: The response vector used in the model.
  \item \code{training_X}: The predictor matrix used in the model.
  \item \code{y_pred}: The fitted values from the final model (no debiasing).
  \item \code{y_pred_debiased}: Debiased fitted values if \code{debias_fit} is available.
  \item \code{formula}: The formula used for model fitting.
  \item \code{terms}: The terms object extracted from the model frame.
}
}
\description{
A wrapper function for \code{\link[glmnet]{cv.glmnet}} that takes input arguments in a manner similar to \code{\link{SVEMnet}}.
This function searches over multiple \code{alpha} values by running \code{cv.glmnet()} for each provided \code{alpha}, and then
selects the combination of \code{alpha} and \code{lambda} with the best cross-validation performance.
}
\details{
This function uses \code{\link[glmnet]{cv.glmnet}} to fit a generalized linear model with elastic net regularization,
performing k-fold cross-validation to select the regularization parameter \code{lambda}. If multiple \code{alpha} values are
provided, it selects the best-performing \code{alpha}-\code{lambda} pair based on the minimal cross-validation error.

After fitting, the function calculates a debiasing linear model (if possible). This is done by regressing the actual responses
on the fitted values obtained from the selected model. The resulting linear model is stored in \code{debias_fit}.
}
\examples{
set.seed(0)
n <- 50
X1 <- runif(n)
X2 <- runif(n)
y <- 1 + 2*X1 + 3*X2 + rnorm(n)
data <- data.frame(y, X1, X2)

model_cv <- glmnet_with_cv(y ~ X1 + X2, data = data, glmnet_alpha = c(0,0.5,1))
predictions <- predict_cv(model_cv, data)

}
\references{
Friedman, J., Hastie, T., & Tibshirani, R. (2010). Regularization Paths for Generalized Linear Models via Coordinate Descent.
\emph{Journal of Statistical Software}, 33(1), 1-22. \doi{10.18637/jss.v033.i01}
}
\seealso{
\code{\link[glmnet]{glmnet}}, \code{\link[glmnet]{cv.glmnet}}, \code{\link{SVEMnet}}
}
