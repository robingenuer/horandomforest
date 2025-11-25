#' Plot of the evolution of the bias of forests and trees against the number of leaves of trees, with confidence intervals for each bias estimation and for slopes and intercepts of linear fits to curves
#'
#'
#' @param forest.comput.bias.res Resulting object from `comput_bias()` function for forests (with `tree = FALSE` in the call of the function)
#' @param forest.confint.res Resulting object from `comput_confint()` function applied to the content of `forest.comput.bias.res`
#' @param tree.comput.bias.res Resulting object from `comput_bias()` function for trees (with `tree = TRUE` in the call of the function)
#' @param tree.confint.res Resulting object from `comput_confint()` function applied to the content of `tree.comput.bias.res`
#'
#' @returns Nothing is returned, only a plot is made.
#' @export
#'
#' @import ggplot2
#' @import scales
#' @importFrom viridisLite viridis
#'
#' @examples
#' forest_bias_result <- comput_bias(nbleaves = c(2^5, 2^6), freg.name = "sinus",
#'    xdim = 1, nbobs = 640, nbobs_test = 60, nfor = 10, var_estim = TRUE,
#'    mc.cores = 2)
#' forest_confint_result <- comput_confint(forest_bias_result, nbBootstrap = 250)
#' tree_bias_result <- comput_bias(nbleaves = c(2^5, 2^6), freg.name = "sinus",
#'    xdim = 1, nbobs = 640, nbobs_test = 60, nfor = 10, tree = TRUE,
#'    var_estim = TRUE, mc.cores = 2)
#' tree_confint_result <- comput_confint(tree_bias_result, nbBootstrap = 250)
#' plot_bias(forest_bias_result, forest_confint_result,
#'     tree_bias_result, tree_confint_result)
plot_bias <- function(forest.comput.bias.res, forest.confint.res,
  tree.comput.bias.res, tree.confint.res) {

  for.bias <- data.frame(forest.comput.bias.res$for.bias.res)
  for.slopesCI <- forest.confint.res$slopesCI
  for.interceptsCI <- forest.confint.res$interceptsCI
  for.bchapsCI <- forest.confint.res$biasesCI

  tree.bias <- data.frame(tree.comput.bias.res$for.bias.res)
  tree.slopesCI <- tree.confint.res$slopesCI
  tree.interceptsCI <- tree.confint.res$interceptsCI
  tree.bchapsCI <- tree.confint.res$biasesCI

  if (!all.equal(tree.bias[, 1], for.bias[, 1])) {
    stop("number of leaves must be the same for tree and forest")
  } else {
    bias <- cbind(tree.bias[, -1], for.bias[, -1])
  }
  means <- seq(1, ncol(bias), 2)
  sds <- seq(2, ncol(bias), 2)
  log.bias <- log(bias, base = 2)
  log.nbleaves <- log(tree.bias$nbleaves, base = 2)

  pentes <- apply(log.bias[, means], MARGIN = 2,
                  FUN = function(y, x = log.nbleaves) {
                    return(lm(y ~ x, data = data.frame(x, y))$coef[2])})
  pentes <- round(pentes, 3)

  intercepts <- apply(log.bias[, means], MARGIN = 2,
                      FUN = function(y,x = log.nbleaves) {
                        return(lm(y ~ x, data = data.frame(x, y))$coef[1])})
  # consts <- round(2^intercepts, 3)
  intercepts <- round(intercepts, 3)

  slopesCI <- round(rbind(tree.slopesCI, for.slopesCI), 3)
  rownames(slopesCI) <- c("treeSquaredErr", "forestSquaredErr")

  interceptsCI <- round(rbind(tree.interceptsCI, for.interceptsCI), 3)
  rownames(interceptsCI) <- c("treeSquaredErr", "forestSquaredErr")

  tree.bias[, 3:4] <- tree.bchapsCI
  colnames(tree.bias)[3:4] <- c("Binf", "Bsup")
  tree.bias$Type <- "tree"
  for.bias[, 3:4] <- for.bchapsCI
  colnames(for.bias)[3:4] <- c("Binf", "Bsup")
  for.bias$Type <- "forest"
  bias <- rbind(tree.bias, for.bias)
  .x <- NULL # avoid a NOTE on package check

  pbiasBars <- ggplot2::ggplot(
    data = bias, ggplot2::aes
    (x = .data$nbleaves, y = .data$mean, color = .data$Type)) +
    ggplot2::geom_point() + ggplot2::geom_line() +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$Binf, ymax = .data$Bsup), width = 0.2) +
    ggplot2::scale_y_continuous(trans = scales::log2_trans(),
                       breaks = scales::trans_breaks("log2", function(x) 2^x),
                       labels = scales::trans_format(
                         "log2", scales::label_math(2^.x))) +
    ggplot2::scale_x_continuous(trans = scales::log2_trans(),
                       breaks = scales::trans_breaks("log2", function(x) 2^x),
                       labels = scales::trans_format(
                         "log2", scales::label_math(2^.x))) +
    ggplot2::scale_color_manual(
      values = viridisLite::viridis(3)[1:2],
      breaks = c("tree", "forest"),
      labels = c(paste0("tree", "\n",
                        "r=", pentes[1], "(",
                        paste(slopesCI[1, ], collapse = ","), ")", "\n",
                        "i=", intercepts[1], "(",
                        paste(interceptsCI[1, ], collapse = ","), ")"),
                 paste0("forest", "\n",
                        "r=", pentes[2], "(",
                        paste(slopesCI[2, ], collapse = ","), ")" , "\n",
                        "i=", intercepts[2], "(",
                        paste(interceptsCI[2, ], collapse = ","), ")"))) +
    ggplot2::xlab("Number of leaves of trees (log-scale)") +
    ggplot2::ylab("Bias (log-scale)") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank())

  print(pbiasBars)
}

