#' @title Calculate scores based on LASSO coefficient_ simplified
#'
#' @description Calculate scores based on LASSO coefficient using expression matrix
#'
#' @param exp matrix, expression matrix
#'
#' @param lassocoef data.frame, table of coefficient, result from lass_cox
#'
#' @param fd data.frame, feature data including id and gene symbol mapping
#'
#' @return matrix, rownames are sample ids, first coloumn is score
#'
#' @examples lassocoef <- lasso_cox(exp, feature, surv)
#' score <- calc_score("GSE7440", lassocoef, destdir = ".", log = TRUE)
#' score2 <- calc_scores_simply(score$exp, score$feature, lassocoef)
#'
#' @import GEOquery
#'
#' @import tidyverse
#'
#' @import edgeR
#'
#' @import magrittr
#'
#' @importFrom rlang .data
#'
#' @export

calc_scores_simply <- function(exp, fd, lassocoef) {
  idx <- fd$`Gene Symbol` %in% lassocoef$covariate
  # all(rownames(exp) == fd$ID)
  if (all(rownames(exp) == fd$ID)) {
    exp <-
      exp[idx,] %>% data.frame() %>% mutate(symbol = fd$`Gene Symbol`[idx]) %>%
      group_by(.data$symbol) %>% summarise_all(mean)
  } else
    exp <-
      exp[fd$ID[idx],] %>% data.frame() %>% mutate(symbol = fd$`Gene Symbol`[idx]) %>%
      group_by(.data$symbol) %>% summarise_all(mean)

  score <-
    t(exp[,-1] %>% mutate_all(list( ~ replace(., is.na(
      .
    ), 0)))) %*%
    lassocoef$beta[match(exp$symbol, lassocoef$covariate)]
}
