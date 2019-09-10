#' @title generate LASSO-Cox model using input data
#'
#' @description Use input expression, clinical and gene list to generate LASSO-Cox model
#'
#' @param exp data.frame, expression matrix, log transformed, first column is patient id (PID), other column is gene
#'
#' @param survdata data.frame, include columnd os.time and os.status,first column must be identical as first column in exp
#'
#' @param feature vector, genes names, consitent with expression data frame colnames
#'
#' @return List: coefficient and lambda (if required)
#'
#' @examples lassocoef <- lasso_cox(x, f, y, lambda = T)
#'
#' @import glmnet
#'
#' @import tidyverse
#'
#' @import survival
#'
#' @import magrittr
#'
#' @importFrom rlang .data
#'
#' @export

lasso_cox <- function(exp, feature, surv, lambda = FALSE) {
  require(glmnet)
  require(tidyverse)
  require(survival)

  ## sanity check
  cln.exp <- colnames(exp)
  cln.surv <- colnames(surv)
  # browser()
  stopifnot(
    ## feature in x
    all(feature %in% cln.exp),
    ## x and y include same patients and same order
    all.equal(
      exp %>% pull(cln.exp[1]) %>% as.character(),
      surv %>% pull(cln.surv[1]) %>% as.character()
    )
  )
  survmat <- Surv(surv$os.time, surv$os.status)
  cvfit <- cv.glmnet(exp[, feature] %>% as.matrix(), survmat,
                     family = "cox")
  coef <- coef(cvfit, s = "lambda.min") %>% as.vector()
  lambda.min <- cvfit$lambda.min
  # browser()
  lassocoef <-
    data.frame(covariate = colnames(exp[, feature])[which(coef != 0)],
               beta = coef[which(coef != 0)])
  if (lambda) {
    return(list(coef=lassocoef,
                lambdamin=lambda.min,
                cvmin=cvfit$cvm[which.min(cvfit$lambda)]))
  } else
    return(lassocoef)
}
