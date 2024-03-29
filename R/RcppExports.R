# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

score_cpp <- function(ipar, resp, th, p, sigma, maxIter = 30L, conv = 0.001, D = 1.7, Fisher = TRUE) {
    .Call('_MAT_score_cpp', PACKAGE = 'MAT', ipar, resp, th, p, sigma, maxIter, conv, D, Fisher)
}

SCORE_cpp <- function(ipar, resp_full, p, sigma, maxIter = 30L, conv = 0.001, D = 1.7, Fisher = TRUE) {
    .Call('_MAT_SCORE_cpp', PACKAGE = 'MAT', ipar, resp_full, p, sigma, maxIter, conv, D, Fisher)
}

selectItem_cpp <- function(ipar, available, given, th, p, sigma, D = 1.7, method = "D", selectionType = "FISHER", c_weights = 0L, content_balancing = FALSE, topN = 1L) {
    .Call('_MAT_selectItem_cpp', PACKAGE = 'MAT', ipar, available, given, th, p, sigma, D, method, selectionType, c_weights, content_balancing, topN)
}

