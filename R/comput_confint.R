#' Computation of confidence intervals of biases, and slopes and intercepts estimations by Bootstrap
#'
#' @param comput.bias.res Resulting object from `comput_bias()` function.
#' @param nbBootstrap Number of bootstrap re-samplings used to get the confidence intervals (default to 500).
#' @param alpha The probability to not belong to the confidence interval (default to 0.05).
#'
#' @returns A list containing the slopes, intercepts and biases confidence intervals.
#' @export
#'
#' @importFrom stats lm
#'
#' @examples
#' bias_result <- comput_bias(nbleaves = c(2^5, 2^6), freg.name = "sinus",
#'    xdim = 1, nbobs = 640, nbobs_test = 60, nfor = 10, var_estim = TRUE,
#'    mc.cores = 2)
#' comput_confint(bias_result, nbBootstrap = 250)
comput_confint <- function(
    comput.bias.res, nbBootstrap = 500, alpha = 0.05) {

  resamplings <- lapply(1:nbBootstrap, function(i) {
    indU <- sample(
      1:length(comput.bias.res$for.bias[[1]]),
      size = length(comput.bias.res$for.bias[[1]]), replace = TRUE)
    indX <- sample(1:nrow(comput.bias.res$for.bias[[1]][[1]]),
                   size = nrow(comput.bias.res$for.bias[[1]][[1]]),
                   replace = TRUE)
    bias <- t(sapply(comput.bias.res$for.bias, function(biask) {
      subBiask <- lapply(biask[indU], function(subB) {
        subB[indX, , drop = FALSE]})
      subSquaredErrorsk <- do.call("rbind", subBiask)
      return(colMeans(subSquaredErrorsk))
    }))
    log.bias <- log(bias, base = 2)
    log.nbleaves <- log(comput.bias.res$for.bias.res[, "nbleaves"], base = 2)
    slopes <- stats::lm(log.bias[1, ] ~ log.nbleaves)$coef[2]
    intercepts <- stats::lm(log.bias[1, ] ~ log.nbleaves)$coef[1]
    return(list(bchaps = bias, slopes = slopes, intercepts = intercepts))
  })

  resampledBchaps <- sapply(resamplings, function(resamp) { resamp$bchaps })
  bchapsCI <- t(apply(resampledBchaps, MARGIN = 1, FUN = function(bchapsk) {
    sortedBchapsk <- sort(bchapsk)
    bchapsBounds <- sortedBchapsk[c(floor(nbBootstrap * alpha / 2),
                                    ceiling(nbBootstrap * (1 - alpha / 2)))]
    return(bchapsBounds)
  }))

  resampledSlopes <- sapply(resamplings, function(resamp) resamp$slopes)
  slopesCI <- sort(resampledSlopes)[c(floor(nbBootstrap * alpha / 2),
                                      ceiling(nbBootstrap * (1 - alpha / 2)))]

  resampledIntercepts <- sapply(resamplings, function(resamp) resamp$intercepts)
  interceptsCI <- sort(resampledIntercepts)[c(floor(nbBootstrap * alpha / 2),
                                     ceiling(nbBootstrap * (1 - alpha / 2)))]
  return(list(slopesCI = slopesCI, interceptsCI = interceptsCI,
              biasesCI = bchapsCI))
}
