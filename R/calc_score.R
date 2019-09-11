#' @title Calculate scores based on LASSO coefficient
#'
#' @description Calculate scores based on LASSO coefficient using data from GEO
#'
#' @param GEO character, GEO accession number
#'
#' @param lassocoef data.frame, table of coefficient, result from lass_cox
#'
#' @param fd data.frame (optional), feature data including id and gene symbol mapping
#'
#' @param is.count logical, default is FALSE. is the GEO data count matrix (RNA-seq) or log transformed expression (microarray or RNA-seq)
#'
#' @param log logical, default is FALSE, whether or not need log transform. Depends on the normalization of expression matrix
#'
#' @param destdir where to store downloaded GEO data, default is tempdir()
#'
#' @return List of data, include scores for each sample, feature data, pheno data, expression matrix for all and subset genes
#'
#' @examples lassocoef <- lasso_cox(exp, feature, surv)
#' score <- calc_score("GSE7440", lassocoef, destdir = ".", log = TRUE)
#'
#' @import GEOquery
#'
#' @import tidyverse
#'
#' @import Biobase
#'
#' @import edgeR
#'
#' @import magrittr
#'
#' @importFrom rlang .data
#'
#' @export

calc_score <- function(GEO, lassocoef=lassocoef, fd=NULL, is.count=FALSE, log = FALSE,
                       destdir = tempdir()) {
  ## updated 2019-04-19
  # library(GEOquery)
  # library(edgeR)
  # browser()
  data1 <- getGEO(GEO, destdir = destdir)
  data1 <- data1[[1]]

  # library(tidyverse)
  pd <- pData(data1)
  if (is.null(fd)) {
    fd <- fData(data1)
  }
  exp.all <- exprs(data1)
  if (is.count) {
    exp.all <- cpm(exp.all, log=T)
  }
  if (log) {
    exp.all <- apply(exp.all,2,log2)
  }
  ## add batch correct?
  # if (combat) {
  #   exp.all
  # }
  idx <- fd$`Gene Symbol` %in% lassocoef$covariate
  # all(rownames(exp) == fd$ID)
  if (all(rownames(exp.all) == fd$ID)) {
    exp.sub <- exp.all[idx, ] %>% data.frame() %>% mutate(symbol = fd$`Gene Symbol`[idx]) %>%
      group_by(.data$symbol) %>% summarise_all(mean)
  } else
    exp.sub <- exp.all[fd$ID[idx], ] %>% data.frame() %>% mutate(symbol = fd$`Gene Symbol`[idx]) %>%
    group_by(.data$symbol) %>% summarise_all(mean)


  # exp <- exp %>% data.frame() %>% mutate(symbol = fd$`Gene Symbol`) %>%
  #   group_by(symbol) %>% summarise_all(mean)
  #
  # idx <- match(lassocoef$covariate, exp$symbol)

  score <- t(exp.sub[, -1] %>% mutate_all(list(~replace(., is.na(.), 0)))) %*%
    lassocoef$beta[match(exp.sub$symbol, lassocoef$covariate)]
  # score <- data.frame(score=score, sample=rownames(score), group=pd$source_name_ch1)

  # score <- score %>% separate(., group, c("leukemia","respond","radiation"), remove = F, sep = ", ")
  return(list(score=score, pheno=pd, exp=exp.all, feature=fd, exp.sub=exp.sub))
}
