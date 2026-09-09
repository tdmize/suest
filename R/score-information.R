# Central differences of analytic score sums avoid cancellation from second
# differences of the log likelihood. Preserve curvature away from the optimum.
.suest_score_information <- function(parameters, score_sum) {
  step <- .Machine$double.eps^(1/3) * pmax(1, abs(parameters))
  H <- vapply(seq_along(parameters), function(j) {
    plus <- minus <- parameters
    plus[j] <- plus[j] + step[j]
    minus[j] <- minus[j] - step[j]
    (score_sum(plus) - score_sum(minus))/(2*step[j])
  }, numeric(length(parameters)))
  -(H + t(H))/2
}
