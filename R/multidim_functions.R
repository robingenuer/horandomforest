# regression function
freg <- function(x, freg.name) {
  if (freg.name == "sum") {
    y <- rowSums(x)
  }
  
  if (freg.name=="fried1") {
    if (length(x[1, ]) < 5) {
      stop("x must have at least 5 coordinates")
    } else {     
      y <- 10 * sin(pi * x[, 1] * x[, 2]) + 20 * (x[, 3] - 0.5)^2 +
        10 * x[, 4] + 5 * x[, 5]
      y <- y / 10
    }
  }
  
  if (freg.name=="sinus") {
     y <-  sin(2 * pi * x)
  }
  
  if (freg.name=="abs") {
    y <- abs(x - 0.5)
  }
  
  return(y)
}

FREG <- NULL

FREG <- function(x, freg.name) {
  if (freg.name=="sinus") {
    y <- - cos(2*pi*x) / (2*pi)
  }  
  return(y)
}

# Function computing the mean value of the regression function in an
# hyper-rectangle h for Fried1 regression function:
beta_fried1 <- function(h) {
  cosc <- function(x) {
    1 / (pi * x) * (cos(pi * h[1, 1] * x) - cos(pi * h[1, 2] * x))
  }
  res <- 10 / ((h[2, 2] - h[2, 1]) * (h[1, 2] - h[1, 1])) *
    integrate(cosc, lower = h[2, 1], upper = h[2, 2])$value +
    20/3 * ((h[3, 2] - 1/2)^2 + (h[3, 2] - 1/2) * (h[3, 1] - 1/2) +
              (h[3, 1] - 1/2)^2) +
    5 * (h[4, 1] + h[4, 2]) + 5/2 * (h[5, 1] + h[5, 2])
  return(res / 10)
}



# Data simulation
simul_data <- function(n, d, freg.name, sigma = 1/4) {
  x <- matrix(runif(d * n), ncol = d)
  eps <- rnorm(n, mean = 0, sd = sigma)
  y <- freg(x, freg.name) + eps
  return(list(x = x, y = y))
}


# Toy model: partition generation (dimension 1)
# k is the parameter of toy model (induces k+1 leaves)
# toy_partition <- function(k) {
#   part <- c(0, seq(from = 1/k, to = 1, by = 1/k) - runif(1) / k, 1)
#   return(part)
# }

# Finds in which interval along each dimension observation z is in (for toy)
intfind <- function(z, d, part) {
  int <- sapply(1:d, FUN = function(j) {
    max(which(part[, j] < z[j]))
    })
  return(int)
}

# Descent of an observation in a tree
descent <- function(z, forest, ind_tree, variante, freg.name,
                    tree_type = "final") {
  if (variante == "toy") {
    d <- ifelse(freg.name == "fried1",
                ncol(forest[[ind_tree]]$beta_margi$beta_margi_others),
                ncol(forest[[ind_tree]]$beta_margi))
    part <- forest[[ind_tree]]$part
    node <- intfind(z = z, d = d, part = part)
  } else {
    node <- 1
    
    if (tree_type == "final") {
      status <- -1
    }
    if (tree_type == "building") {
      status <- 0
    }
    
    while (forest$nodestatus[node, ind_tree] != status) {
      if (z[forest$bestvar[node, ind_tree]] <=
          forest$xbestsplit[node, ind_tree]) {
        node <- forest$leftDaughter[node, ind_tree]
      } else {
        node <- forest$rightDaughter[node, ind_tree]
      }
    }
  }
  return(node)
}


# Bounds computation associated to a node of a tree (of a forest)
interval_bounds <- function(
    forest, ind_tree, node, dim.int = forest$bestvar[node, ind_tree]) {
  
  bounds <- matrix(NA, nrow = length(dim.int), ncol = 2)
  ancetres.gauche <- NULL
  ancetres.droit <- NULL
  root <- node
  
  if (node == 1) {
    bounds[, 1] <- 0
    bounds[, 2] <- 1
  } else {
    
    while (root != 1) {
      ancetres.gauche <- c(ancetres.gauche,
                           which(forest$leftDaughter[, ind_tree] == root))
      ancetres.droit <- c(ancetres.droit,
                          which(forest$rightDaughter[, ind_tree] == root))
      root <- min(ancetres.gauche[length(ancetres.gauche)],
                  ancetres.droit[length(ancetres.droit)])
    }
    
    dim_nb <- 1
    for (dim in dim.int) {
      ancetres.gauche.dim <-
        ancetres.gauche[which(forest$bestvar[ancetres.gauche, ind_tree] == dim)]
      ancetres.droit.dim <-
        ancetres.droit[which(forest$bestvar[ancetres.droit, ind_tree] == dim)]
      
      bsup <- 1
      binf <- 0
      if (length(ancetres.gauche.dim) >0) {
        bsup <- forest$xbestsplit[max(ancetres.gauche.dim), ind_tree]
      }
      if (length(ancetres.droit.dim) > 0) {        
        binf <- forest$xbestsplit[max(ancetres.droit.dim), ind_tree]
      }
      bounds[dim_nb, ] <- c(binf, bsup)
      dim_nb <- dim_nb + 1
    }
  }
  return(bounds)
}


# Forest structure generation
forest_structure <- function(
    k, d, variante, freg.name, ntree, nbobs = NULL, mtry = max(1, floor(d/3)),
    bootstrap_part = TRUE, mc.cores = 1, seed = NULL) {
  
  if (variante == "toy") {
    forest <- replicate(
      ntree,
      tree_tilde(k = k, d = d, variante = variante, freg.name = freg.name),
      simplify=FALSE)
  }

  if (variante == "bb") {
    if (!is.null(seed)) set.seed(seed)
    
    z <- simul_data(n = nbobs, d = d, freg.name = freg.name)
    x <- z$x
    y <- z$y
    
    if (bootstrap_part == TRUE) {
      if (d == 1) {
        rf <- randomForest::randomForest(
          y ~ x, data = z, ntree = ntree, maxnodes = k, nodesize = 1, mtry = mtry)
      } else {
        rf <- randomForest::randomForest(
          x = x, y = y, ntree = ntree, maxnodes = k, nodesize = 1, mtry = mtry)
      }
    } else {
      rf <- randomForest::randomForest(
        x = x, y = y, ntree = ntree, maxnodes = k, nodesize = 1, mtry = mtry,
        replace = FALSE, sampsize = nrow(x))
    }
    forest <- rf$forest
  }
  return(forest)
}


# Hyper-rectangles coordinates computation
hyper_rec_base <- function(ind_tree, forest, hrec, d) {
  hrec.base <- hrec[, ind_tree, , , drop = FALSE]
  nrnodes <- forest$nrnodes
  tree <- lapply(forest[-(1:2)], FUN = function(x)
    subset(x, select = ind_tree, subset = rep(TRUE, nrnodes)))
  
  for (i in 2: nrnodes) {
    hrec.base[i, , , ] <-
      interval_bounds(
        forest = forest, ind_tree = ind_tree, node = i, dim.int = 1:d)
  }
  
  return(hrec.base)
}

hyper_rec <- function(forest, d, mc.cores = 1) {
  ntree <- forest$ntree
  nrnodes <- forest$nrnodes
  hrec <- array(data = NA, dim=c(nrnodes, ntree, d, 2))
  hrec[1, , , 1] <- 0
  hrec[1, , , 2] <- 1    
  hrec_bases_list <- parallel::mclapply(
    1:ntree, hyper_rec_base, forest = forest, hrec = hrec, d = d,
    mc.cores = mc.cores)
  for (j in 1:ntree) {
    hrec[, j, , ] <- hrec_bases_list[[j]]
  }
  
  return(hrec)
}


# Ideal forest built
forest_tilde <- function(
    k, d, variante, freg.name, ntree, nbobs = NULL, mtry = max(1, floor(d/3)),
    bootstrap_part = TRUE, seed = NULL, mc.cores = 1) {
  
    forest <- forest_structure(
      k = k, d = d, variante = variante, freg.name = freg.name, ntree = ntree,
      nbobs = nbobs, mtry = mtry, bootstrap_part = bootstrap_part,
      mc.cores = mc.cores, seed = seed)
    
    if (variante == "bb") {
      hrec <- hyper_rec(forest = forest, d = d, mc.cores = mc.cores)
      
      if (freg.name == "sum") {
        forest$nodepred <- apply(X = hrec, MARGIN = 1:2, FUN = sum) / 2
      }
      
      if (freg.name == "fried1") {
        forest$nodepred <- apply(X = hrec, MARGIN = 1:2, FUN = beta_fried1)
      }
      
      if (freg.name == "sinus") {
        forest$nodepred <- matrix(apply(X = hrec, MARGIN = 1:3, FUN = function(u) {
          diff(FREG(u, freg.name = freg.name)) / diff(u)
        }), ncol = ntree)
      }
      
      if (freg.name == "abs") {
        forest$nodepred <- matrix(apply(X = hrec, MARGIN = 1:3, FUN = function(u) {
          if (u[1] < 0.5 & u[2] > 0.5) {
            out <- (0.5 + u[2]^2 - u[2] + u[1]^2 - u[1]) / (2*u[2] - u[1])
          } else {
            out <- abs((u[1] + u[2] - 1) / 2)
          }
          return(out)
        }), ncol = ntree)
      }
    }
  
  return(forest)
}


# Predict one observation with a forest
f_estim_base <- function(ind_obs, x, forest, variante, freg.name) {
  pred.trees <- rep(NA, forest$ntree)
  for (ind_tree in 1:forest$ntree) {
    node <- descent(z = x[ind_obs, ], forest = forest, ind_tree = ind_tree,
                    variante = variante, freg.name = freg.name)
    pred.trees[ind_tree] <- forest$nodepred[node, ind_tree]
  }
  return(mean(pred.trees))
}


# Predict all observations of x with a forest
f_estim <- function(x, forest, variante, freg.name, mc.cores=1) {
  if (variante == "toy") {
    all_estims <- simplify2array(parallel::mclapply(
      forest, FUN = function(tree) {
        t_estim(x = x, tree = tree, variante = variante, freg.name = freg.name)
        }, mc.cores = mc.cores))
    estims <- rowMeans(all_estims)
  } else {
    n <- length(x[, 1])
    estims <- unlist(parallel::mclapply(
      1:n, FUN = f_estim_base, x = x, forest = forest, variante = variante,
      freg.name = freg.name, mc.cores=mc.cores))
  }
  return(estims)
}


## TREE
# Ideal tree computation
tree_tilde <- function(
    k, d, variante, freg.name, nbobs = NULL, mtry = max(1, floor(d/3)),
    bootstrap_part = TRUE, seed = NULL) {
  
  if (variante == "toy") {
    k_prod <- k
    if (floor(k_prod) != k_prod) {
      warning("The number of leaves for toy model must be of the
    form (k+1)^d with d the dimension and k an integer")
    }
    part <- matrix(NA, nrow = k_prod + 2, ncol = d)
    beta_margi <- matrix(NA, nrow = k_prod + 1, ncol = d)
    
    for (j in 1:d) {
      part[, j] <- toy_partition(k_prod)
      if (freg.name == "sum") {
        beta_margi[, j] <- (part[1:(k_prod + 1), j] + part[2:(k_prod + 2), j]) / 2
      }
    }
    
    if (freg.name == "fried1") {
      beta_1st_plan <- matrix(NA, nrow = k_prod + 1, ncol = k_prod + 1)
      cosc <- function(x) {
        1/(pi*x) * (cos(pi*part[leaf_x1, 1] * x) -
                      cos(pi*(part[leaf_x1, 1] + 1) * x))
      }
      for (leaf_x1 in 1:(k_prod + 1)) {
        for (leaf_x2 in 1:(k_prod + 1)) {
          beta_1st_plan[leaf_x1, leaf_x2] <- 
            (part[leaf_x2 + 1, 2] - part[leaf_x2, 2]) *
            (part[leaf_x1 + 1, 1] - part[leaf_x1, 1]) *
            integrate(cosc, lower = part[leaf_x2, 2],
                      upper = part[leaf_x2 + 1, 2])$value 
        }
      }
      beta_margi[, 3] <- 20/3 *
        ((part[2:(k_prod + 2), 3] - 1/2)^2 +
           (part[2:(k_prod + 2), 3] - 1/2) * (part[1:(k_prod + 1), 3] - 1/2) +
           (part[1:(k_prod + 1), 3] - 1/2)^2)
      beta_margi[, 4] <- 5 * (part[1:(k_prod + 1), 4] + part[2:(k_prod + 2), 4])
      beta_margi[, 5] <- 5/2 * (part[1:(k_prod + 1), 5] + part[2:(k_prod + 2), 5])
      
      beta_margi <- list(beta_1st_plan = 1/10 * beta_1st_plan,
                         beta_margi_others = 1/10 * beta_margi)
    }
    
    tree <- list(part = part, beta_margi = beta_margi)
  } else {
    tree <- forest_tilde(
      k = k, d = d, variante = variante, freg.name = freg.name, ntree = 1,
      nbobs = nbobs, mtry = mtry, bootstrap_part = bootstrap_part,
      seed = seed, mc.cores = 1)
  }
  
  return(tree)
}


# Predict all observations of x with one tree
t_estim <- function(x, tree, variante, freg.name) {
  if (variante == "toy") {
    part <- tree$part
    
    if (freg.name == "sum") {
      beta_margi <- tree$beta_margi
    } else {
      beta_margi <- tree$beta_margi$beta_margi_others
    }

    d <- ncol(beta_margi)
    y_margi <- matrix(NA, nrow = nrow(x), ncol = d)
    ints <- apply(x, MARGIN = 1, FUN = intfind, d = d, part = part)
    y_margi <- mapply(FUN = function(i, j) {beta_margi[i, j]}, ints, 1:d)
    y_margi <- matrix(y_margi, ncol = d, byrow = TRUE)
    
    if (freg.name == "fried1") {
      beta_1st_plan <- tree$beta_margi$beta_1st_plan
      y_1stplan <- mapply(FUN = function(i, j) {
        beta_1st_plan[i, j]}, ints[1, ], ints[2, ])
    }
    
    if (freg.name == "sum") {
      estims <- rowSums(y_margi)
    } else {
      estims <- rowSums(y_margi[, -(1:2)]) + y_1stplan
    }
  } else {
    estims <- f_estim(
      x = x, forest = tree, variante = variante, freg.name = freg.name)
  }
  
  return(estims)
}


# Forest bias computation
bias_for <- function(
    nbleaves, variante, freg.name, d, nbobs = NULL, beta = 1, nfor = 1,
    r = 500, b = 1, mtry = max(1, floor(d/3)), tree = FALSE,
    bootstrap_part = TRUE, var_estim = FALSE, seeds = FALSE, mc.cores = 1,
    mc.cores2 = 1) {
  
  debut <- Sys.time()
  m <- 1
  for.bias <- vector(mode = "list", length = length(nbleaves))
  for.bias.res <- matrix(NA, nrow = length(nbleaves), ncol = 4)
  x.test <- matrix(runif(d*r), ncol = d)
  s <- freg(x.test, freg.name)
  if (seeds) {seeds <- dget("outputs/seeds")}
  
  for (k in nbleaves) {
    debutk <- Sys.time()
    if (variante == "toy") {
      x.test.bl <- matrix(runif(d*r, min = b/k, max = 1 - b/k), ncol = d)
    } else {
      x.test.bl <- matrix(runif(d*r, min = b*log(k^(1/d)) / (k^(1/d)),
                                max = 1 - b*log(k^(1/d)) / (k^(1/d))), ncol = d)
    }
    
    ntree <- ifelse(tree, 1, ifelse(
      variante == "toy", floor(((k+1)^d)^beta), max(100, floor(k^beta))))
    
    for.bias[[m]] <- parallel::mclapply(1:nfor, function(j) {
      forest <- forest_tilde(k = k, d = d, variante = variante,
                             freg.name = freg.name, ntree = ntree,
                             nbobs = nbobs, mtry = mtry,
                             bootstrap_part = bootstrap_part,
                             # seed = ifelse(seeds, seed[m, j], NULL),
                             mc.cores = mc.cores)
      s.tilde <- f_estim(x = x.test, forest = forest, variante = variante,
                         freg.name = freg.name, mc.cores = mc.cores)
      # mse <- mean((s - s.tilde)^2)
      squaredErrors <- (s - s.tilde)^2
      
      s.bl <- freg(x.test.bl, freg.name)
      s.tilde.bl <- f_estim(x = x.test.bl, forest = forest, variante = variante,
                            freg.name = freg.name, mc.cores = mc.cores)
      # mse.bl <- mean((s.bl - s.tilde.bl)^2)
      squaredErrors.bl <- (s.bl - s.tilde.bl)^2
      
      # return((mse, mse.bl))
      return(cbind(squaredErrors, squaredErrors.bl))
    }, mc.cores = mc.cores2)
    
    allSquaredErrors <- do.call("rbind", for.bias[[m]])
    for.bias.res[m, c(1, 3)] <- colMeans(allSquaredErrors)
    
    if (var_estim) {
      allCrossedSquaredErrors <- tcrossprod(allSquaredErrors[, 1])
      for (i in 1:nfor) {
        allCrossedSquaredErrors[1:r + (i-1)*r, 1:r + (i-1)*r] <- NA
      }
      allCrossedSquaredErrors[lower.tri(allCrossedSquaredErrors)] <- NA
      
      for.bias.res[m, 2] <-
        sqrt(for.bias.res[m, 1]^2 - mean(allCrossedSquaredErrors, na.rm = TRUE))
      
      allCrossedSquaredErrorsBL <- tcrossprod(allSquaredErrors[, 2])
      for (i in 1:nfor) {
        allCrossedSquaredErrorsBL[1:r + (i-1)*r, 1:r + (i-1)*r] <- NA
      }
      
      for.bias.res[m, 4] <-
        sqrt(for.bias.res[m, 3]^2 - mean(allCrossedSquaredErrorsBL, na.rm = TRUE))
    }
    
    # indPairs <- seq(1, nrow(allSquaredErrors), by = r+1)
    # indCrossed <- expand.grid(indPairs, indPairs)
    # allIndCrossed <- indCrossed[-indPairs, ]
    # 
    # allIndCrossed <- NULL
    # for (l in 1:10) {
    #   indPairs <- seq(l, nrow(allSquaredErrors), by = r+1)
    #   indCrossed <- expand.grid(indPairs, indPairs)
    #   allIndCrossed <- rbind(allIndCrossed, indCrossed[-indPairs,])
    # }
    
    # allCrossedSquaredErrors <- t(mapply(function(i, j) {
    #   allSquaredErrors[i,] * allSquaredErrors[j, ]},
    #   allIndCrossed$Var1, allIndCrossed$Var2))
    # 
    # for.bias.res[m, c(2, 4)] <-
    #   for.bias.res[m, c(1, 3)]^2 - colMeans(allCrossedSquaredErrors)
    
    # for.bias.res[m, c(2, 4)] <- apply(for.bias[[m]], MARGIN = 1, sd)
    if (tree) {
      print(paste(variante, "bias tree k =", k, "d =", d, "mtry = ", mtry,
                  "bootstrap_part = ", bootstrap_part, "done"))
    } else {
      print(paste(variante, "bias forest k =", k, "d =", d, "mtry = ", mtry,
                  "bootstrap_part = ", bootstrap_part, "done"))
    }
    print(Sys.time() - debutk)
    m <- m+1
  }
  
  names(for.bias) <- nbleaves
  for.bias.res <- cbind(nbleaves, for.bias.res)
  colnames(for.bias.res) <- c("nbleaves", "mean", "sd", "BLmean", "BLsd")
  total_time <- Sys.time() - debut
  time_unit <- attr(total_time, "units")
  print(paste0("Total time = ", total_time))
  save(
    for.bias.res, for.bias, freg.name, variante, d, nbleaves, nbobs, beta, nfor,
    r, b, mtry, tree, bootstrap_part, total_time, mc.cores2,
    file = paste0("outputs/",
                  paste(freg.name, variante, "dim", d,
                        "nbleaves", paste(nbleaves, collapse = "_"),
                        "nbobs", nbobs, "beta", beta, "nfor", nfor, "mtry",
                        mtry, sep="_"),
                  ifelse(tree, "_TREE", "_FOREST"),
                  "_BIAS.Rdata"))
  return(for.bias)
}

# Slope computation
slope <- function(y, x) {
  data <- data.frame(x, y)
  return(lm(y ~ x, data = data)$coef[2])
}

# Intercept computation
intercept <- function(y, x) {
  data <- data.frame(x, y)
  return(lm(y ~ x, data = data)$coef[1])
}


CIslope <- function(
    for.bias, nbleaves, variante, freg.name, d, nbobs = NULL, beta = 1,
    tree = FALSE, nfor = 1, nbSamples = 500, alpha = 0.05) {
  
  resamplings <- lapply(1:nbSamples, function(i) {
    indU <- sample(
      1:length(for.bias[[1]]), size = length(for.bias[[1]]), replace = TRUE)
    indX <- sample(1:nrow(for.bias[[1]][[1]]), size = nrow(for.bias[[1]][[1]]),
                   replace = TRUE)
    bias <- t(sapply(for.bias, function(biask) {
      subBiask <- lapply(biask[indU], function(subB) { subB[indX, ] })
      subSquaredErrorsk <- do.call("rbind", subBiask)
      return(colMeans(subSquaredErrorsk))
    }))
    log.bias <- log(bias, base = 2)
    log.nbleaves <- log(nbleaves, base = 2)
    slopes <- apply(log.bias, MARGIN = 2, FUN = slope, x = log.nbleaves)
    intercepts <- apply(log.bias, MARGIN = 2, FUN = intercept, x = log.nbleaves)
    return(list(bchaps = bias, slopes = slopes, intercepts = intercepts))
  })
  
  resampledBchaps <- sapply(resamplings, function(resamp) { resamp$bchaps[, 1] })
  bchapsCI <- t(apply(resampledBchaps, MARGIN = 1, FUN = function(bchapsk) {
    sortedBchapsk <- sort(bchapsk)
    bchapsBounds <- sortedBchapsk[c(floor(nbSamples * alpha / 2),
                                    ceiling(nbSamples * (1 - alpha / 2)))]
    return(bchapsBounds)
  }))
  
  resampledBchapsBL <- sapply(resamplings, function(resamp) { resamp$bchaps[, 2] })
  bchapsBLCI <- t(apply(resampledBchapsBL, MARGIN = 1, FUN = function(bchapsk) {
    sortedBchapsk <- sort(bchapsk)
    bchapsBounds <- sortedBchapsk[c(floor(nbSamples * alpha / 2),
                                    ceiling(nbSamples * (1 - alpha / 2)))]
    return(bchapsBounds)
  }))
  
  resampledSlopes <- sapply(resamplings, function(resamp) { resamp$slopes })
  slopesCI <- t(apply(resampledSlopes, MARGIN = 1, FUN = function(slopes) {
    sortedSlopes <- sort(slopes)
    bounds <- sortedSlopes[c(floor(nbSamples * alpha / 2),
                             ceiling(nbSamples * (1 - alpha / 2)))]
    return(bounds)
  }))
  
  resampledIntercepts <- sapply(resamplings, function(resamp) { resamp$intercepts })
  interceptsCI <- t(apply(resampledIntercepts, MARGIN = 1, FUN = function(intercepts) {
    sortedIntercepts <- sort(intercepts)
    bounds <- sortedIntercepts[c(floor(nbSamples * alpha / 2),
                                 ceiling(nbSamples * (1 - alpha / 2)))]
    return(bounds)
  }))
  load(file = paste0("outputs/",
                     paste(freg.name, variante, "dim", d,
                           "nbleaves", paste(nbleaves, collapse = "_"),
                           "nbobs", nbobs, "beta", beta, "nfor", nfor, "mtry",
                           mtry, sep="_"),
                     ifelse(tree, "_TREE", "_FOREST"),
                     "_BIAS.Rdata"))
  save(for.bias.res, for.bias, freg.name, variante, d, nbleaves, nbobs, beta,
       nfor, r, b, mtry, tree, bootstrap_part, total_time, mc.cores2, slopesCI,
       interceptsCI, bchapsCI, bchapsBLCI,
       file = paste0("outputs/",
                     paste(freg.name, variante, "dim", d,
                           "nbleaves", paste(nbleaves, collapse = "_"),
                           "nbobs", nbobs, "beta", beta, "nfor", nfor, "mtry",
                           mtry, sep="_"),
                     ifelse(tree, "_TREE", "_FOREST"),
                     "_BIAS.Rdata"))
  return(list(slopesCI = slopesCI, interceptsCI = interceptsCI,
              bchapsCI = bchapsCI, bchapsBLCI = bchapsBLCI))
}

