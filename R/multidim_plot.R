library(ggplot2)
library(scales)

# Slope computation
slope <- function(y, x) {
  data <- data.frame(y)
  return(lm(y ~ x, data = data)$coef[2])
}

# Intercept computation
intercept <- function(y, x) {
  data <- data.frame(y)
  return(lm(y ~ x, data = data)$coef[1])
}

bias_plot <- function(nbleaves, variante, freg.name, d, nbobs, beta, nfor_for,
                      nfor_tree, b, mtry, bootstrap_part = TRUE, BL = TRUE) {
  
  load(file = paste0(
    "outputs/", paste(
      freg.name, variante, "dim", d, "nbleaves", paste(nbleaves, collapse = "_"),
      "nbobs", nbobs, "beta", beta, "nfor", nfor_tree, "mtry", mtry, sep="_"),
    "_TREE", "_BIAS.Rdata"))
  tree.bias <- data.frame(for.bias.res)
  tree.slopesCI <- slopesCI
  tree.interceptsCI <- interceptsCI
  tree.bchapsCI <- bchapsCI
  load(file = paste0(
    "outputs/", paste(
      freg.name, variante, "dim", d, "nbleaves", paste(nbleaves, collapse = "_"),
      "nbobs", nbobs, "beta", beta, "nfor", nfor_for, "mtry", mtry, sep="_"),
    "_FOREST", "_BIAS.Rdata"))
  for.bias <- data.frame(for.bias.res)
  for.slopesCI <- slopesCI
  for.interceptsCI <- interceptsCI
  for.bchapsCI <- bchapsCI
    
  if (!BL) {
    tree.bias <- tree.bias[, -(4:5)]
    for.bias <- for.bias[, -(4:5)]
  }
  if (!all.equal(tree.bias[, 1], for.bias[, 1])) {
    stop("number of leaves must be the same for tree and forest")
  } else {
    bias <- cbind(tree.bias[, -1], for.bias[, -1])
  }
  means <- seq(1, ncol(bias), 2)
  sds <- seq(2, ncol(bias), 2)
  log.bias <- log(bias, base = 2)
  if (variante %in% c("toy", "purfprod") & d > 1) {
    tree.bias$nbleaves <- (tree.bias$nbleaves + 1)^d
    for.bias$nbleaves <- (for.bias$nbleaves + 1)^d
  }
  log.nbleaves <- log(tree.bias$nbleaves, base = 2)
  
  pentes <- apply(log.bias[, means], MARGIN=2, FUN=slope, x = log.nbleaves)
  pentes <- round(pentes, 3)
  
  intercepts <- apply(log.bias[, means], MARGIN=2, FUN=intercept,
                      x = log.nbleaves)
  consts <- round(2^intercepts, 3)
  intercepts <- round(intercepts, 3)
  
  slopesCI <- round(rbind(tree.slopesCI, for.slopesCI), 3)
  rownames(slopesCI) <- c("treeSquaredErr", "treeSquaredErrBL",
                          "forestSquaredErr", "forestSquaredErrBL")
  
  interceptsCI <- round(rbind(tree.interceptsCI, for.interceptsCI), 3)
  rownames(interceptsCI) <- c("treeSquaredErr", "treeSquaredErrBL",
                              "forestSquaredErr", "forestSquaredErrBL")
  if (!BL) {
    slopesCI <- slopesCI[-c(2, 4), ]
    interceptsCI <- interceptsCI[-c(2, 4), ]
  }
  
  if (BL) {
    tree.biasBF <- tree.bias[, 1:3]
    tree.biasBF$Type <- "tree"
    tree.biasBL <- tree.bias[, c(1, 4:5)]
    tree.biasBL$Type <- "tree BL"
    colnames(tree.biasBL)[2:3] <- c("mean", "sd")
    tree.bias <- rbind(tree.biasBF, tree.biasBL)
    
    for.biasBF <- for.bias[, 1:3]
    for.biasBF$Type <- "forest"
    for.biasBL <- for.bias[, c(1, 4:5)]
    for.biasBL$Type <- "forest BL"
    colnames(for.biasBL)[2:3] <- c("mean", "sd")
    for.bias <- rbind(for.biasBF, for.biasBL)
    bias <- rbind(tree.bias, for.bias)
    
    pbiasBars <- ggplot(data = bias, aes(x = nbleaves, y = mean, color = Type)) +
      geom_point() + geom_line() +
      geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +
      scale_y_continuous(trans = log2_trans(),
                         breaks = trans_breaks("log2", function(x) 2^x),
                         labels = trans_format("log2", math_format(2^.x))) +
      scale_x_continuous(trans = log2_trans(),
                         breaks = trans_breaks("log2", function(x) 2^x),
                         labels = trans_format("log2", math_format(2^.x))) +
      scale_color_manual(
        values = viridis::viridis(4)[1:4],
        breaks = c("tree", "tree BL", "forest", "forest BL"),
        labels = c(paste0("tree", "\n",
                          "r=", pentes[1], "(",
                          paste(slopesCI[1, ], collapse = ","), ")"),
                   paste0("tree BL", "\n",
                          "r=", pentes[2], "(",
                          paste(slopesCI[2, ], collapse = ","), ")"),
                   paste0("forest", "\n",
                          "r=", pentes[3], "(",
                          paste(slopesCI[3, ], collapse = ","), ")"),
                   paste0("forest BL", "\n",
                          "r=", pentes[4],"(",
                          paste(slopesCI[4, ], collapse = ","), ")"))) +
      xlab("Number of leaves of trees (log-scale)") + ylab("Bias (log-scale)") +
      theme_bw() +
      theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
  } else {
    tree.bias[, 3:4] <- tree.bchapsCI
    colnames(tree.bias)[3:4] <- c("Binf", "Bsup")
    tree.bias$Type <- "tree"
    for.bias[, 3:4] <- for.bchapsCI
    colnames(for.bias)[3:4] <- c("Binf", "Bsup")
    for.bias$Type <- "forest"
    bias <- rbind(tree.bias, for.bias)
    
    pbiasBars <- ggplot(data = bias, aes(x = nbleaves, y = mean, color = Type)) +
      geom_point() + geom_line() +
      geom_errorbar(aes(ymin = Binf, ymax = Bsup), width = 0.2) +
      scale_y_continuous(trans = log2_trans(),
                         breaks = trans_breaks("log2", function(x) 2^x),
                         labels = trans_format("log2", math_format(2^.x))) +
      scale_x_continuous(trans = log2_trans(),
                         breaks = trans_breaks("log2", function(x) 2^x),
                         labels = trans_format("log2", math_format(2^.x))) +
      scale_color_manual(
        values = viridis::viridis(3)[1:2],
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
      xlab("Number of leaves of trees (log-scale)") + ylab("Bias (log-scale)") +
      theme_bw() +
      theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
  }
  print(pbiasBars)
  ggsave(filename = paste0(
    "fig/", paste(
      freg.name, variante, "dim", d, "nbleaves", paste(nbleaves, collapse = "_"),
      "nbobs", nbobs, "betaRound", round(beta), "nfor_for", nfor_for,
      "nfor_tree", nfor_tree, "border", b, "mtry", mtry, sep="_"), "_BIAS.pdf"),
    height=3, width=5)
}

