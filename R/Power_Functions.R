################################## Functions for power analysis #############################################

################ Simplifying NB distributions:
power.groups.NB <- function (n1, n2, pars1, pars2, 
                             distribution = "NegBin", zero.inflation = FALSE,
                             test = c("permute", "KW"), 
                             alternative = "two", alpha = 0.05, 
                             nsims = 1000, nreps = 999) 
{
  if (zero.inflation == FALSE) {
    # Negbin parameters:
    m1 = pars1[1]
    m2 = pars2[1]
    size1 = pars1[2]
    size2 = pars2[2]
    # Starting count:
    count = 0
    if (test == "KW") {
      count_min = 0
      count_max = 0
    }
    for (j in 1:nsims) {
      y1 = stats::rnbinom(n = n1, size = size1, mu = m1)
      y2 = stats::rnbinom(n = n2, size = size2, mu = m2)
      if (test == "permute") {
        p = emon::permute.groups(y1, y2, alternative = alternative, 
                                 nreps = nreps)$p.val
      }
      if (test == "KW") {
        TestDat <- data.frame("Y" = c(y1, y2),
                              "Group" = as.factor(c(rep("1", length(y1)), rep("2", length(y2)))))
        p = coin::pvalue(coin::kruskal_test(Y ~ Group, data = TestDat, alpha = alpha,
                                            distribution = coin::approximate(nresample = nreps)))
        p_min <- attributes(p)$conf.int[1]
        p_max <- attributes(p)$conf.int[2]
      }
      if (p < alpha) 
        count = count + 1
      if(test  == "KW") {
        if (p_min < alpha) {
          count_min = count_min +1
        }
        if (p_max < alpha) {
          count_max = count_max +1
        }
      }
    }
    power = count/nsims
    if(test  == "KW") {
      power_min <- count_min/nsims
      power_max <- count_max/nsims
    }
  }
  if (zero.inflation == TRUE) {
    # Proportitions of 0s:
    P0_1 = pars1[1]
    P0_2 = pars2[1]
    # Negbin parameters:
    m1 = pars1[2]
    m2 = pars2[2]
    size1 = pars1[3]
    size2 = pars2[3]
    # Starting count:
    count = 0
    if (test == "KW") {
      count_min = 0
      count_max = 0
    }
    for (j in 1:nsims) {
      y1 = ifelse(stats::rbinom(n1, size = 1, prob = P0_1) > 0, 0, 
                  stats::rnbinom(n = n1, size = size1, mu = m1))
      y2 = ifelse(stats::rbinom(n2, size = 1, prob = P0_2) > 0, 0, 
                  stats::rnbinom(n = n2, size = size2, mu = m2))
      if (test == "permute") {
        p = emon::permute.groups(y1, y2, alternative = alternative, 
                                 nreps = nreps)$p.val
      }
      if (test == "KW") {
        TestDat <- data.frame("Y" = c(y1, y2),
                              "Group" = as.factor(c(rep("1", length(y1)), rep("2", length(y2)))))
        p = coin::pvalue(coin::kruskal_test(Y ~ Group, data = TestDat, alpha = alpha,
                                            distribution = coin::approximate(nresample = nreps)))
        p_min <- attributes(p)$conf.int[1]
        p_max <- attributes(p)$conf.int[2]
      }
      if (p < alpha) 
        count = count + 1
      if(test  == "KW") {
        if (p_min < alpha) {
          count_min = count_min +1
        }
        if (p_max < alpha) {
          count_max = count_max +1
        }
      }
    }
    power = count/nsims
    if(test  == "KW") {
      power_min <- count_min/nsims
      power_max <- count_max/nsims
    }
  }
  if(test == "KW") {
    power_results <- list()
    power_results[["power"]] <- power
    power_results[["power_alphamin"]] <- power_min
    power_results[["power_alphamax"]] <- power_max
    return(power_results)
  }
  return(power)
}

################ Implementing beta distribution and ZI beta:
power.groups.beta <- function (n1, n2, pars1, pars2, distribution = "beta", zero.inflation = FALSE,
                             test = c("permute", "KW"), alternative = "two", alpha = 0.05, 
                             nsims = 1000, nreps = 999) 
{
  if (zero.inflation == FALSE) { 
    # Shape of beta distribution:
    shape1_1 = pars1[2]
    shape1_2  = pars2[2]
    shape2_1 = pars1[3]
    shape2_2  = pars2[3]
    # Startng count:
    count = 0
    if (test == "KW") {
      count_min = 0
      count_max = 0
    }
    for (j in 1:nsims) {
      # Group 1:
      y1 = stats::rbeta(n1, shape1 = shape1_1, shape2 = shape2_1)
      # Group 2:
      y2 = stats::rbeta(n2, shape1 = shape1_2, shape2 = shape2_2)
      if (test == "permute") {
        p = emon::permute.groups(y1, y2, alternative = alternative, 
                                 nreps = nreps)$p.val
      }
      if (test == "KW") {
        TestDat <- data.frame("Y" = c(y1, y2),
                              "Group" = as.factor(c(rep("1", length(y1)), rep("2", length(y2)))))
        p = coin::pvalue(coin::kruskal_test(Y ~ Group, data = TestDat, alpha = alpha,
                                            distribution = coin::approximate(nresample = nreps)))
        p_min <- attributes(p)$conf.int[1]
        p_max <- attributes(p)$conf.int[2]
      }
      if (p < alpha) 
        count = count + 1
      if(test  == "KW") {
        if (p_min < alpha) {
          count_min = count_min +1
        }
        if (p_max < alpha) {
          count_max = count_max +1
        }
      }
    }
    power = count/nsims
    if(test  == "KW") {
      power_min <- count_min/nsims
      power_max <- count_max/nsims
    }
  }
  if (zero.inflation == TRUE) { 
    # Proportitions of 0s:
    P0_1 = pars1[1]
    P0_2 = pars2[1]
    # Shape of beta distribution:
    shape1_1 = pars1[2]
    shape1_2  = pars2[2]
    shape2_1 = pars1[3]
    shape2_2  = pars2[3]
    # Startng count:
    count = 0
    if (test == "KW") {
      count_min = 0
      count_max = 0
    }
    for (j in 1:nsims) {
      # Group 1:
      y1 = ifelse(stats::rbinom(n1, size = 1, prob = P0_1) > 0, 0, 
                  stats::rbeta(n1, shape1 = shape1_1, shape2 = shape2_1))
      # Group 2:
      y2 = ifelse(stats::rbinom(n2, size = 1, prob = P0_2) > 0, 0, 
                  stats::rbeta(n2, shape1 = shape1_2, shape2 = shape2_2))
      if (test == "permute") {
        p = emon::permute.groups(y1, y2, alternative = alternative, 
                                 nreps = nreps)$p.val
      }
      if (test == "KW") {
        TestDat <- data.frame("Y" = c(y1, y2),
                              "Group" = as.factor(c(rep("1", length(y1)), rep("2", length(y2)))))
        p = coin::pvalue(coin::kruskal_test(Y ~ Group, data = TestDat, alpha = alpha,
                                            distribution = coin::approximate(nresample = nreps)))
        p_min <- attributes(p)$conf.int[1]
        p_max <- attributes(p)$conf.int[2]
      }
      if (p < alpha) 
        count = count + 1
      if(test  == "KW") {
        if (p_min < alpha) {
          count_min = count_min +1
        }
        if (p_max < alpha) {
          count_max = count_max +1
        }
      }
    }
    power = count/nsims
    if(test  == "KW") {
      power_min <- count_min/nsims
      power_max <- count_max/nsims
    }
  }
  if(test == "KW") {
    power_results <- list()
    power_results[["power"]] <- power
    power_results[["power_alphamin"]] <- power_min
    power_results[["power_alphamax"]] <- power_max
    return(power_results)
  }
  return(power)
}


################# Implementing multimodal beta & NB:

power.groups.bimode <- function (n1, n2,
                                 pars1 = list(mod1 = NULL, mod2 = NULL, prop = NULL), 
                                 pars2 = list(mod1 = NULL, mod2 = NULL, prop = NULL), 
                                 distribution = c("NegBin", "beta"), 
                                 test = c("permute", "KW"), alternative = "two", alpha = 0.05, 
                                 nsims = 1000, nreps = 999) 
{
  if (distribution == "NegBin") {
    ## Negbin parameters:
    # First group:
    n1_mod1 <- round(n1*pars1$prop[1])
    n1_mod2 <- round(n1*pars1$prop[2])
    Fmod1_mu <- pars1$mod1[1]
    Fmod2_mu <- pars1$mod2[1]
    Fmod1_size <- pars1$mod1[2]
    Fmod2_size <- pars1$mod2[2]
    # Second group:
    n2_mod1 <- round(n2*pars2$prop[1])
    n2_mod2 <- round(n2*pars2$prop[2])
    Smod1_mu <- pars2$mod1[1]
    Smod2_mu <- pars2$mod2[1]
    Smod1_size <- pars2$mod1[2]
    Smod2_size <- pars2$mod2[2]
    # Starting count:
    count = 0
    if (test == "KW") {
      count_min = 0
      count_max = 0
    }
    for (j in 1:nsims) {
      # First group:
      y1_mod1 <- stats::rnbinom(n = n1_mod1, size = Fmod1_size, mu = Fmod1_mu)
      y1_mod2 <- stats::rnbinom(n = n1_mod2, size = Fmod2_size, mu = Fmod2_mu)
      y1 <- c(y1_mod1, y1_mod2)
      # Second group:
      y2_mod1 <- stats::rnbinom(n = n2_mod1, size = Smod1_size, mu = Smod1_mu)
      y2_mod2 <- stats::rnbinom(n = n2_mod2, size = Smod2_size, mu = Smod2_mu)
      y2 <- c(y2_mod1, y2_mod2)
      if (test == "permute") {
        p = emon::permute.groups(y1, y2, alternative = alternative, 
                                 nreps = nreps)$p.val
      }
      if (test == "KW") {
        TestDat <- data.frame("Y" = c(y1, y2),
                              "Group" = as.factor(c(rep("1", length(y1)), rep("2", length(y2)))))
        p = coin::pvalue(coin::kruskal_test(Y ~ Group, data = TestDat, alpha = alpha,
                                            distribution = coin::approximate(nresample = nreps)))
        p_min <- attributes(p)$conf.int[1]
        p_max <- attributes(p)$conf.int[2]
      }
      if (p < alpha) 
        count = count + 1
      if(test  == "KW") {
        if (p_min < alpha) {
          count_min = count_min +1
        }
        if (p_max < alpha) {
          count_max = count_max +1
        }
      }
    }
    power = count/nsims
    if(test  == "KW") {
      power_min <- count_min/nsims
      power_max <- count_max/nsims
    }
  }
  if (distribution == "beta") { 
    # Shape of beta distribution:
    # First group:
    n1_mod1 <- round(n1*pars1$prop[1])
    n1_mod2 <- round(n1*pars1$prop[2])
    Fmod1_shape1 <- pars1$mod1[1]
    Fmod2_shape1 <- pars1$mod2[1]
    Fmod1_shape2 <- pars1$mod1[2]
    Fmod2_shape2 <- pars1$mod2[2]
    # Second group:
    n2_mod1 <- round(n2*pars2$prop[1])
    n2_mod2 <- round(n2*pars2$prop[2])
    Smod1_shape1 <- pars2$mod1[1]
    Smod2_shape1 <- pars2$mod2[1]
    Smod1_shape2 <- pars2$mod1[2]
    Smod2_shape2 <- pars2$mod2[2]
    # Startng count:
    count = 0
    if (test == "KW") {
      count_min = 0
      count_max = 0
    }
    for (j in 1:nsims) {
      # Group 1:
      y1_mod1  <- stats::rbeta(n1_mod1, shape1 = Fmod1_shape1, shape2 = Fmod1_shape2)
      y1_mod2  <- stats::rbeta(n1_mod2, shape1 = Fmod2_shape1, shape2 = Fmod2_shape2)
      y1 <- c(y1_mod1, y1_mod2)
      # Group 2:
      y2_mod1  <- stats::rbeta(n2_mod1, shape1 = Smod1_shape1, shape2 = Smod1_shape2)
      y2_mod2  <- stats::rbeta(n2_mod2, shape1 = Smod2_shape1, shape2 = Smod2_shape2)
      y2 <- c(y2_mod1, y2_mod2)
      if (test == "permute") {
        p = emon::permute.groups(y1, y2, alternative = alternative, 
                                 nreps = nreps)$p.val
      }
      if (test == "KW") {
        TestDat <- data.frame("Y" = c(y1, y2),
                              "Group" = as.factor(c(rep("1", length(y1)), rep("2", length(y2)))))
        p = coin::pvalue(coin::kruskal_test(Y ~ Group, data = TestDat, alpha = alpha,
                                            distribution = coin::approximate(nresample = nreps)))
        p_min <- attributes(p)$conf.int[1]
        p_max <- attributes(p)$conf.int[2]
      }
      if (p < alpha) 
        count = count + 1
      if(test  == "KW") {
        if (p_min < alpha) {
          count_min = count_min +1
        }
        if (p_max < alpha) {
          count_max = count_max +1
        }
      }
    }
    power = count/nsims
    if(test  == "KW") {
      power_min <- count_min/nsims
      power_max <- count_max/nsims
    }
  }
  if(test == "KW") {
    power_results <- list()
    power_results[["power"]] <- power
    power_results[["power_alphamin"]] <- power_min
    power_results[["power_alphamax"]] <- power_max
    return(power_results)
  }
  return(power)
}



