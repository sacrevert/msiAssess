# Deterministic small example to mirror the JAGS recursions
recur_true  <- function(estimate, spgrowth, FY) {
  T <- length(estimate); stopifnot(length(spgrowth) == T - 1)
  spindex <- rep(NA_real_, T)
  # anchor (Y1perfect = TRUE)
  spindex[FY] <- estimate[FY]
  # forward
  if (FY < T) for (t in (FY+1):T) spindex[t] <- estimate[FY] + sum(spgrowth[FY:(t-1)])
  # backward
  if (FY > 1) for (t in (FY-1):1) spindex[t] <- spindex[t+1] - spgrowth[t]
  spindex
}

recur_false <- function(estimate, spgrowth, FY) {
  T <- length(estimate); stopifnot(length(spgrowth) == T - 1)
  spindex <- rep(NA_real_, T)
  # backward *includes* t = FY (uses spindex[FY+1])
  # forward defines spindex[FY+1], but notice the boundary case FY = T
  if (FY < T) {
    # define spindex[FY+1]..spindex[T]
    for (t in (FY+1):T) spindex[t] <- estimate[FY] + sum(spgrowth[FY:(t-1)])
    # define spindex[FY] down to 1
    for (t in FY:1) spindex[t] <- spindex[t+1] - spgrowth[t]
  } else {
    # FY == T: forward loop is empty; backward will need spindex[T+1] (undefined)
    stop("Boundary bug: needs spindex[FY+1] but FY == nyears.")
  }
  spindex
}

## 1) Interior anchor: Y1perfect TRUE vs FALSE give identical results
T  <- 6
FY <- 3
estimate <- c(NA, NA, 2, NA, NA, NA)
spgrowth <- rep(0.1, T - 1)

a <- recur_true(estimate, spgrowth, FY)
b <- recur_false(estimate, spgrowth, FY)
all.equal(a, b)  # TRUE

## 2) Boundary anchor: FY == nyears triggers the bug in the FALSE variant
FY <- T
try(recur_false(estimate = c(rep(NA, T-1), 2), spgrowth = rep(0.1, T-1), FY = FY))
# Error: "Boundary bug: needs spindex[FY+1] but FY == nyears."



