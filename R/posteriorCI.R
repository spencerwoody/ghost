
posteriorCI <- function(draws, terms, alpha = 0.05) {

  est <- rowMeans(draws)
  lo <- apply(draws, 1, quantile, probs = alpha / 2)
  hi <- apply(draws, 1, quantile, probs = 1 - alpha / 2)

  data.frame(
    terms = terms,
    est = est,
    lo = lo,
    hi = hi
  )

}
