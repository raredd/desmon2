### from desmon if package not available
# bin1samp, binci, simon, twocon, twostg
### 


bin1samp <- function(p0, pa, alpha = 0.1, beta = 0.1, n.min = 20L) {
  ## desmon::bin1samp
  if (p0 == pa)
    stop('p0 == pa')
  b <- 1
  x <- round(p0*n.min)
  n <- n.min - 1L
  
  if (pa > p0) {
    while (b > beta) {
      n <- n + 1
      ## determine cutoff: reject if X > x
      l <- x:n
      s <- 1 - pbinom(l, n, p0)
      sub <- s <= alpha
      x <- l[sub][1L]
      size <- s[sub][1L]
      b <- pbinom(x, n, pa)
    }
  } else if (pa < p0) {
    while (b > beta) {
      n <- n + 1L
      ## determine cutoff: reject if X < x
      l <- x:0
      s <- pbinom(l, n, p0)
      sub <- s <= alpha
      x <- l[sub][1L]
      size <- s[sub][1L]
      b <- pbinom(x, n, pa, lower.tail = FALSE)
      x <- x + 1L
    }
  }
  
  structure(
    c(n = n, r = x, p0 = p0, pa = pa, size = size, type2 = b),
    class = 'bin1samp'
  )
}

binci <- function(r, n, conf = 0.95) {
  ## desmon::binci
  if (r < 0 | r > n)
    stop('invalid value for r')
  if (conf <= 0 | conf >= 1)
    stop('must have 0 < conf < 1')
  
  ff <- function(p, r, n, alpha) {
    pbinom(r - 1L, n, p, lower.tail = FALSE) - alpha
  }
  ff2 <- function(p, r, n, alpha) {
    pbinom(r, n, p) - alpha
  }
  
  alpha <- (1 - conf) / 2
  pl <- if (r <= 0)
    0 else uniroot(ff, c(1e-8, 1 - 1e-8), r = r, n = n, alpha = alpha)$root
  pu <- if (r >= n)
    1 else uniroot(ff2, c(1e-8, 1 - 1e-8), r = r, n = n, alpha = alpha)$root
  
  c(pl, pu)
}

simon <- function(p0, pa, n1max = 0, ntmax = 1e5, alpha = 0.1, beta = 0.1,
                  del = 1L, minimax = FALSE) {
  ## desmon::simon
  # prob calculations for 2-stage phase II design
  # optimal simon designs (minimize E(n|H0))
  # p0, pa, null & alt response probabilities
  # n1max=max # subject in first stage
  if (alpha > 0.5 | alpha <= 0 | 1 - beta <= alpha | beta <= 0 | p0 <= 0 |
      p0 >= pa | pa >= 1 | n1max > ntmax)
    stop('invalid arguments')
  
  # determine min sample size first stage
  n1min <- max(ceiling(log(beta) / log(1 - pa)), 2)
  if (n1min > ntmax)
    stop('no valid designs')
  
  # optimal one-sample design
  u1 <- bin1samp(p0, pa, alpha, beta)
  
  #include even if n>n1max
  z <- matrix(c(u1[1:2], 0, u1[2L], 1, u1[5:6], u1[1L]), nrow = 1L)
  if (n1min < u1[1L]) {
    n1max <- if (n1max > 0)
      min(n1max, u1[1L] - 1) else u1[1L] - 1
  } else if (n1min == u1[1L]) {
    n1max <- if (n1max > 0)
      min(n1max, u1[1L]) else u1[1L]
  } else {
    stop('no valid designs')
  }
  n1max <- min(n1max, ntmax)
  # determine min total sample size (use randomized decision rule on 
  # boundary for single sample)
  b <- 0
  n <- u1[1L]
  while (b <= beta) {
    n <- n - 1
    u <- c(1, pbinom(0:(u1[2]+1), n, p0, lower.tail = FALSE))
    index <- min((1:(u1[2L] + 3))[u <= alpha])
    pi <- (alpha - u[index]) / (u[index - 1] - u[index])
    if (index > 2) {
      u <- pbinom((index-3):(index-2), n, pa, lower.tail = FALSE)
    } else {
      u <- c(1, pbinom(0, n, pa, lower.tail = FALSE))
    }
    b <- 1 - u[2L] - pi * (u[1L] - u[2L])
  }
  if (n > ntmax)
    stop('no valid designs')
  e0 <- u1[1L]
  for (i in n1min:n1max) { # cases 1st stage
    # stage I
    # feasible stopping rules
    x1 <- 0:i
    w1 <- dbinom(x1, i, p0)
    w2 <- dbinom(x1, i, pa)
    sub <- cumsum(w2) <= beta
    y1 <- x1[sub] # feasible 1st stage stopping rules
    for (r1 in y1) {
      j <- n - i # additional cases at 2nd stage
      if (j > 0) {
        q <- 0
        pi0 <- 1 - sum(w1[1:(r1 + 1)])
        while (q == 0) {
          x2 <- 0:j
          w3 <- dbinom(x2, j, p0)
          w4 <- dbinom(x2, j, pa)
          # b0,ba=marginal dist of # responses H0, Ha
          u3 <- c(outer(x1[-(1:(r1 + 1))], x2, `+`))
          u4 <- c(outer(w1[-(1:(r1 + 1))], w3))
          b0 <- cumsum(c(w1[1:(r1 + 1)], tapply(u4, u3, sum)))
          sub <- b0 < 1 - alpha
          r2 <- sum(sub)
          u4 <- c(outer(w2[-(1:(r1 + 1))], w4))
          ba <- cumsum(c(w2[1:(r1 + 1)], tapply(u4, u3, sum)))
          e1 <- i + pi0 * j
          if (ba[r2 + 1] <= beta) {
            q <- 1
            if (minimax) {
              if (i + j < ntmax) {
                ntmax <- i + j
                e0 <- e1 # need to reset to min E.H0 at min # cases
              } else if (i + j == ntmax)
                e0 <- min(e0, e1)
            } else {
              e0 <- min(e0, e1)
            }
            z <- rbind(z, c(i, r1, j, r2, b0[r1 + 1L], 1 - b0[r2 + 1], ba[r2 + 1], e1))
          } else {
            j <- j + 1
            if (minimax) {
              if (i + j > ntmax)
                q <- 1 
            } else {
              if (e1 > e0 + del | i + j > ntmax)
                q <- 1
            }
          }
        }
      }
    }
  }
  dimnames(z) <- list(NULL, c('n1', 'r1', 'n2', 'r2', 'Pstop1.H0',
                              'size', 'type2', 'E.tot.n.H0'))
  z <- z[z[, 1L] + z[, 3L] <= ntmax, , drop = FALSE]
  z <- z[z[, 8L] <= e0 + del, , drop = FALSE]
  
  list(
    designs = z[order(z[, 8L]), , drop = FALSE],
    call = match.call(),
    description = c(
      'n1, n2 = cases 1st stage and additional # in 2nd',
      'r1, r2 = max # responses 1st stage and total to declare trt inactive')
  )
}

twocon <- function(n1, n2, r1, r, conf = 0.95, dp = 1) {
  ## desmon::twocon
  # confidence interval for two-stage phase II
  # n1=# entered at first stage, n2=additional # entered at second stage
  # r1=max number responses stage 1 trt 1 and still stop
  # r= total number responses observed
  if (n1 < 1 | n2 < 1 | r1 < 0 | r1 > n1 | r < 0 |
      r > n2 + n1 | conf <= 0 | conf >= 1)
    stop('invalid arguments')
  alpha <- (1 - conf) / 2
  x1 <- 0:n1
  x2 <- 0:n2
  u1 <- c(outer(x1[-(1:(r1 + 1))], x2, `+`))
  
  dbin2 <- function(p1, x1, x2, u1, n1, n2, r1) {
    w1 <- dbinom(x1, n1, p1)
    w3 <- dbinom(x2, n2, p1)
    u2 <- c(outer(w1[-(1:(r1 + 1))], w3))
    # ith component is P(R=i-1)
    c(w1[1:(r1 + 1)], tapply(u2, u1, sum))
  }
  ff <- function(p, x1, x2, u1, n1, n2, r1, mm, dbin2, mle) {
    sum(dbin2(p, x1, x2, u1, n1, n2, r1) * mm) - mle
  }
  ff2 <- function(p, x1, x2, u1, n1, n2, r1, s, dbin2, alpha) {
    sum(dbin2(p, x1, x2, u1, n1, n2, r1)[s]) - alpha
  }
  
  n <- n1 + n2
  mle <- if (r > r1)
    r / n else r / n1
  # bias corrected estimator
  mm <- c((0:r1) / n1, (r1 + 1):n / n)
  pm <- if (r <= 0)
    0 else if (r >= n)
      1 else uniroot(ff, c(1e-8, 1 - 1e-8), x1 = x1, x2 = x2, u1 = u1,
                     n1 = n1, n2 = n2, r1 = r1, mm = mm, dbin2 = dbin2,
                     mle = mle)$root
  # unbiased estimator
  if (r1 >= r) { 
    ube <- r / n1
  } else {
    aa <- dhyper((r1 + 1):r, n1, n2, r)
    ube <- sum(((r1 + 1):r) * aa) / (n1 * sum(aa))
  }
  mm <- if (r > r1)
    c((n / n1) ^ dp * (0:r1), (r1 + 1):n)
  else c(0:r1,(n1 / n) ^ dp * ((r1 + 1):n))
  s1 <- mm >= r # for P(R>=r)
  s2 <- mm <= r
  pl <- if (r <= 0)
    0 else uniroot(ff2, c(1e-8, 1 - 1e-8), x1 = x1, x2 = x2, u1 = u1,
                   n1 = n1, n2 = n2, r1 = r1, s = s1, dbin2 = dbin2,
                   alpha = alpha)$root
  pu <- if (r >= n)
    1 else uniroot(ff2, c(1e-8, 1 - 1e-8), x1 = x1, x2 = x2, u1 = u1,
                   n1 = n1, n2 = n2, r1 = r1, s = s2, dbin2 = dbin2,
                   alpha = alpha)$root
  
  c(lower = pl, upper = pu, bcmle = pm, mle = mle, unbiased = ube)
}

twostg <- function(n1, n2, p1, r1, r2) {
  ## desmon::twostg
  # prob calculations for 2-stage phase II design
  # n1=# entered at first stage, n2=additional # entered at second stage
  # p1, response probability
  # r1=max number responses stage 1 trt 1 and still stop
  # r2=max total number responses on trt 1 to be inactive
  if (n1 < 1 | n2 < 1 | r1 < 0 | r2 < 0 | p1 <= 0 | r1 > n1 |
      r2 > n2 + n1 | p1 >= 1)
    stop('invalid arguments')
  x1 <- 0:n1
  x2 <- 0:n2
  w1 <- dbinom(x1, n1, p1)
  w3 <- dbinom(x2, n2, p1)
  # b1=marginal dist of # responses in group 1
  u1 <- c(outer(x1[-(1:(r1 + 1))], x2, `+`))
  u2 <- c(outer(w1[-(1:(r1 + 1))], w3))
  b1 <- c(w1[1:(r1 + 1)], tapply(u2, u1, sum))
  
  u <- list(
    inputs = c(n1 = n1, n2 = n2, p1 = p1, r1 = r1, r2 = r2),
    prob.inactive = c(total = sum(b1[1:(r2 + 1L)]), sum(b1[1:(r1 + 1)]))
  )
  
  structure(u, class = 'twostg')
}
