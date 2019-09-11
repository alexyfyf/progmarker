#' @title generate LASSO-Cox model based on occurence in permutations
#'
#' @description Use input expression, clinical and gene list to generate LASSO-Cox model, permutate and obtain most often occured model
#'
#' @param x data.frame, expression matrix, log transformed, first column is patient id (PID), other column is gene
#'
#' @param y data.frame, include columnd os.time and os.status,first column must be identical as first column in exp
#'
#' @param f vector, genes names, consitent with expression data frame colnames
#'
#' @param n numeric, top n models will be output
#'
#' @param perm numeric, number of permutations
#'
#' @return List: coefficient and lambda (if required)
#'
#' @examples
#' \dontrun{
#' top <- topN(exp, feature, surv, perm=100, n=3)
#' }
#'
#' @import glmnet
#'
#' @import tidyverse
#'
#' @import dplyr
#'
#' @import survival
#'
#' @import magrittr
#'
#' @importFrom rlang .data
#'
#' @import future.apply
#'
#' @import future
#'
#' @export

## for reproducibility, not true RNG
topN <- function(x, f, y, perm, n=3) {
  # require(future.apply)
  plan(multiprocess)
  list <- future_lapply(1:perm, function(a) {
    lasso_cox(x, f, y, lambda = T)
  },
  ## here, because I want the sampling to be random for each permutation
  ## set future.seed = 1234
  future.seed = 1234)

  lbd <- sapply(list, function(x)
    x[[2]])
  lbd_c <- lbd %>% round(6) %>% as.character()
  order <- lbd_c %>% table() %>% sort(.data, decreasing = T)
  # return(order[1:3])

  lst <- lapply(list, function(x)
    x[[1]])
  res <- lapply(1:n, function(t) {
    idx <- which(lbd_c == names(order[t]))[1]
    lst[[idx]]
  })
  names(res) <-
    paste(deparse(substitute(x)),
          deparse(substitute(y)),
          deparse(substitute(f)),
          1:n,
          sep = "_")
  return(list(res=res,
              lbd=data.frame(lambda=names(order), freq=as.vector(order))))
}
