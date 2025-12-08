## ================================================================
## Berazneva et al. (2019, AJAE)
## ================================================================

rm(list = ls())

## ----------------------------
## 0) Horizon & discounting
## ----------------------------
T_h   <- 100        
delta <- 0.10
beta  <- 1/(1 + delta)
disc  <- beta^(0:(T_h-1))

## ----------------------------
## 1) Prices & costs (USD)
## ----------------------------
p_maize <- 331      # $/Mg maize
nP      <- 4434     # $/Mg N
qR      <- 58       # $/Mg residues (opportunity cost if not retained)
m_fix   <- 375      # $/ha other/overhead

## ----------------------------
## 2) Biophysical parameters (paper-consistent)
## ----------------------------
k_rg <- 1.50        # residue:grain ratio
F_C  <- 0.43        # carbon fraction of residues
Dmin <- 0.11        # SOC mineralization rate
Bexp <- 0.45        # exponent in SOC formation

## Reported baseline steady state (δ = 10%)
f_ss <- 0.13
a_ss <- 0.54
c_ss <- 25.63
y_ss_target <- 3.91  # paper steady-state yield (≈)

## ----------------------------
## 3) Yield function (quadratic; Table 2 coefficients)
## y = b0 + bc*c + bcc*c^2 + bf*f + bff*f^2 + bcf*c*f
## ----------------------------
b0   <- -1.461
bc   <-  0.113
bcc  <- -0.000413
bf   <-  27.04
bff  <- -41.30
bcf  <- -0.218

yield_fun_raw <- function(c, f) {
  b0 + bc*c + bcc*c^2 + bf*f + bff*f^2 + bcf*c*f
}

# Adjust intercept to include average fixed effects: y(c_ss, f_ss) ≈ 3.91 ## add missing fixed effect so that yield at (c_ss,f_ss) matches paper’s steady-state yield (3.91 Mg/ha)
b0 <- b0 + (y_ss_target - max(1e-9, yield_fun_raw(c_ss, f_ss)))

yield_fun <- function(c, f) {
  pmax(1e-9, b0 + bc*c + bcc*c^2 + bf*f + bff*f^2 + bcf*c*f)
}

## ----------------------------
## 4) Calibrate A from steady state (Eq. 16): A*(aFk*y)^B = D*c
# Calibrate A so steady-state inflow equals outflow:
#   at steady state: A*(a*F*k*y)^B = D*c  =>  A = D*c / (a*F*k*y)^B.
# Paper does not report A directly, so we solve for A using reported (c*,a*,f*).
## ----------------------------
calibrate_A <- function(c_star, f_star, a_star, D, F, k, B) {
  y_star <- yield_fun(c_star, f_star)
  inflow <- (a_star * F * k * y_star)^B
  if (inflow <= 0) stop("Non-positive inflow term; check parameters/yield.")
  D * c_star / inflow
}
Ain <- calibrate_A(c_ss, f_ss, a_ss, Dmin, F_C, k_rg, Bexp)
cat("Calibrated A =", round(Ain, 4),
    " (y_ss ~", round(yield_fun(c_ss, f_ss), 3), "Mg/ha)\n")

## ----------------------------
## 5) Dynamics & profits
## ----------------------------
c_next <- function(c, f, a, D = Dmin, A = Ain, F = F_C, k = k_rg, B = Bexp) {
  inflow <- (pmax(a,0) * F * k * yield_fun(c, f))^B      # c_{t+1} = c_t - D*c_t + A*(a * F * k * y(c,f))^B
  pmax(c - D * c + A * inflow, 0)                        # inflow from retained residue; pmax prevents negative SOC
}

profit_fun <- function(c, f, a) {
  y <- yield_fun(c, f)
  (p_maize + qR * k_rg * (1 - a)) * y - nP * f - m_fix   # π(c,f,a) = (p_maize + qR * k_rg * (1-a)) * y(c,f)  - nP*f - m_fix
}

terminal_penalty <- function(cT, c_ss_target = c_ss, rho_disc = delta, tau = 1.0) {
  - tau * ((1 + rho_disc) / rho_disc) * (c_ss_target - cT)^2
}

## ----------------------------
## 6) Simulator for a path {f_t, a_t}
## ----------------------------
simulate_path <- function(f_vec, a_vec, c0, A = Ain) {
  c <- numeric(T_h + 1); c[1] <- c0  #creates a vector to store soil carbon over all periods (including time 0).
  y <- numeric(T_h) 
  pr <- numeric(T_h) #storage for yield and profit each year.
  for (t in 1:T_h) {
    y[t]  <- yield_fun(c[t], f_vec[t]) #compute yield using current soil carbon and fertilizer.
    pr[t] <- profit_fun(c[t], f_vec[t], a_vec[t]) #compute profit: revenue from maize + residue - fertilizer cost - fixed cost.
    c[t+1] <- c_next(c[t], f_vec[t], a_vec[t], A = A) #update soil carbon using the equation from the paper
  }
  list(c = c, y = y, pr = pr)
}

## ----------------------------
## 7) Time-varying open-loop objective with tail SS penalty
## Decision vector x = [f_1..f_T, a_1..a_T] (unconstrained)
## Map: f_t in [0, 0.24] Mg N/ha; a_t in (0,1) via logistic
##The optimizer chooses the full 50-year path of fertilizer and residue decisions 
##that maximizes the present value of profits, while penalties ensure that the simulated path converges to the steady-state reported in the paper.
## ----------------------------
map_controls <- function(x) {
  f_raw <- x[1:T_h]                   #x is a long vector of decisions.
  a_raw <- x[(T_h+1):(2*T_h)]
  f <- pmin(0.24, pmax(0, f_raw))        # 0..0.24 Mg N/ha (0..240 kg)
  a <- 1/(1 + exp(-a_raw))               # (0,1)
  list(f = f, a = a)
}

objective <- function(x, c0, A = Ain) {
  # weights 
  K_tail    <- 5      # last K years to nudge to steady state
  lambda_g  <- 25     # weight on SS balance residual
  lambda_ca <- 1.0    # mild pull toward (f_ss, a_ss) in tail
  tau_term  <- 1.0    # terminal penalty multiplier
  
  ctrl <- map_controls(x)
  sim  <- simulate_path(ctrl$f, ctrl$a, c0, A)   #Simulates soil carbon, yields, and profits for all years
  pv   <- sum(disc * sim$pr)   #Computes present value of profits
  
  # terminal anchor
  term <- terminal_penalty(sim$c[T_h + 1], c_ss, delta, tau = tau_term)
  
  # tail residuals for steady-state balance: g_t = D*c_t - A*(aFk y)^B  (want ~0)
  idx  <- (T_h - K_tail + 1):T_h
  g    <- Dmin * sim$c[idx] - Ain * (ctrl$a[idx] * F_C * k_rg * yield_fun(sim$c[idx], ctrl$f[idx]))^Bexp
  pen_g  <- - lambda_g * sum(g^2)
  
  # gentle pull of tail controls toward reported (f*, a*) (numerical stabilizer)
  pen_ca <- - lambda_ca * sum( (ctrl$f[idx] - f_ss)^2 + (ctrl$a[idx] - a_ss)^2 )
  
  # minimize negative of (PV + penalties)
  -(pv + term + pen_g + pen_ca)
}

## ----------------------------
# Run optimizer to find the best time path; then map controls, simulate, and compute PV
## 8) Solve baseline (time-varying)
## ----------------------------
set.seed(1)
c0 <- 19.12  # medium initial SOC
x0 <- c(rep(f_ss, T_h), qlogis(rep(a_ss, T_h)))  # start near steady state

opt <- optim(
  par     = x0,
  fn      = objective,
  c0      = c0,
  A       = Ain,
  method  = "L-BFGS-B",
  control = list(maxit = 1000)
)

ctrl <- map_controls(opt$par)
sim  <- simulate_path(ctrl$f, ctrl$a, c0, A = Ain)
PV   <- sum(disc * sim$pr) + terminal_penalty(sim$c[T_h+1], c_ss)

cat("\n=== Baseline (time-varying) ===\n")
cat("Convergence (0=OK):", opt$convergence, "\n")
cat("PV ($/ha):", round(PV, 1), "\n")
cat("Tail means (yrs 41–50): f =", round(mean(ctrl$f[(T_h-9):T_h]), 3),
    " a =", round(mean(ctrl$a[(T_h-9):T_h]), 3),
    " c_T =", round(sim$c[T_h+1], 2),
    " y_T =", round(tail(sim$y,1), 2), "\n")

## ----------------------------
# Simulate current-practice path and compute its PV for comparison (welfare gain)
## 9) Counterfactual: fixed “current practices”
## ----------------------------
f_cur <- rep(0.06, T_h)   # ~60 kg N/ha
a_cur <- rep(0.25, T_h)
sim_cur <- simulate_path(f_cur, a_cur, c0, A = Ain)
PV_cur  <- sum(disc * sim_cur$pr) + terminal_penalty(sim_cur$c[T_h+1], c_ss)

cat("\n=== Counterfactual (fixed current practices) ===\n")
cat("PV_current ($/ha):", round(PV_cur, 1), "\n")
cat("ΔPV (baseline - current) ($/ha):", round(PV - PV_cur, 1), "\n")

## ----------------------------
## 10) Quick table
# Build and print summary table comparing optimized tail means vs current practice
## ----------------------------
tbl <- rbind(
  c("Baseline (tail mean)", round(mean(ctrl$f[(T_h-9):T_h]),3),
    round(mean(ctrl$a[(T_h-9):T_h]),3),
    round(sim$c[T_h+1],2),
    round(tail(sim$y,1),2),
    round(PV,1)),
  c("Current (fixed)", round(mean(f_cur),3),
    round(mean(a_cur),3),
    round(sim_cur$c[T_h+1],2),
    round(tail(sim_cur$y,1),2),
    round(PV_cur,1))
)
colnames(tbl) <- c("Case","f (Mg/ha)","a (share)","c_T (Mg/ha)","y_T (Mg/ha)","PV ($/ha)")
print(tbl, quote = FALSE, row.names = FALSE)

## ----------------------------
## 11) Plots — robust device + PNGs
# Plot SOC, yield, fertilizer, and residue-share paths and save PNGs for figures
## ----------------------------

# Helper to open a big interactive device
open_big_device <- function(w=10, h=7) {
  while (dev.cur() > 1) dev.off()
  if (.Platform$OS.type == "windows") {
    windows(width = w, height = h)
  } else if (Sys.info()[["sysname"]] == "Darwin") {
    quartz(width = w, height = h)
  } else {
    x11(width = w, height = h)
  }
}

# Interactive window
open_big_device(10, 7)
par(mfrow = c(2,2), mar = c(4,4,2,1), mgp = c(2.2,0.7,0))

plot(0:T_h, sim$c, type="l", lwd=2, xlab="Year", ylab="SOC c_t (Mg/ha)", main="Soil Carbon (baseline)")
abline(h = c_ss, lty = 2)

plot(1:T_h, sim$y, type="l", lwd=2, xlab="Year", ylab="Yield (Mg/ha)", main="Yield (baseline)")

plot(1:T_h, ctrl$f, type="l", lwd=2, xlab="Year", ylab="f_t (Mg N/ha)", main="N path (baseline)")
abline(h = f_ss, lty = 3)

plot(1:T_h, ctrl$a, type="l", lwd=2, xlab="Year", ylab="a_t (share)", main="Residue share (baseline)")
abline(h = a_ss, lty = 3)

# PNG fallback (works anywhere)
plot_png <- function(file, expr, width=1400, height=900, res=150) {
  png(file, width=width, height=height, res=res)
  par(mar = c(4,4,2,1), mgp = c(2.2,0.7,0))
  force(expr)
  dev.off()
}

plot_png("soil_carbon_path.png", {
  plot(0:T_h, sim$c, type="l", lwd=2, xlab="Year", ylab="SOC c_t (Mg/ha)", main="Soil Carbon Path")
  abline(h = c_ss, lty = 2)
})

plot_png("yield_path.png", {
  plot(1:T_h, sim$y, type="l", lwd=2, xlab="Year", ylab="Yield (Mg/ha)", main="Yield Path")
})

plot_png("n_path.png", {
  plot(1:T_h, ctrl$f, type="l", lwd=2, xlab="Year", ylab="f_t (Mg N/ha)", main="N path (baseline)")
  abline(h = f_ss, lty = 3)
})

plot_png("a_path.png", {
  plot(1:T_h, ctrl$a, type="l", lwd=2, xlab="Year", ylab="a_t (share)", main="Residue share (baseline)")
  abline(h = a_ss, lty = 3)
})

## ----------------------------
## 12) Steady-state diagnostic at T
# Check steady-state balance at T: inflow A*(aFk*y)^B vs outflow D*c (gap should be ~0)
## ----------------------------
c_last <- sim$c[T_h+1]; f_last <- ctrl$f[T_h]; a_last <- ctrl$a[T_h]
lhs <- Ain * (a_last * F_C * k_rg * yield_fun(c_last, f_last))^Bexp
rhs <- Dmin * c_last
cat("\n=== Steady-state diagnostic (year T) ===\n")
cat("A*(aFk*y)^B =", round(lhs,4), " ; D*c_T =", round(rhs,4),
    " ; gap =", round(lhs - rhs,4), "\n")
