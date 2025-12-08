rm(list = ls())

## ================================================================
## Carbon-Flow Subsidy Extension for Berazneva et al. (2019)
## ================================================================

## ----------------------------
## 0) Horizon
## ----------------------------
T_h   <- 100
delta <- 0.10
beta  <- 1/(1 + delta)
disc  <- beta^(0:(T_h-1))

## ----------------------------
## 1) Prices & Costs
## ----------------------------
p_maize <- 331
nP      <- 4434
qR      <- 58
m_fix   <- 375

## ----------------------------
## 2) Biophysical parameters 
## ----------------------------
k_rg <- 1.50
F_C  <- 0.43
Dmin <- 0.11
Bexp <- 0.52
Ain  <- 2.40     

## ----------------------------
## 3) Steady-state values from paper 
## ----------------------------
f_ss <- 0.13
a_ss <- 0.54
c_ss <- 25.63
y_ss_target <- 3.91

## ----------------------------
## 4) Quadratic yield function
## ----------------------------
b0   <- -1.461
bc   <-  0.113
bcc  <- -0.000413
bf   <-  27.04
bff  <- -41.30
bcf  <- -0.218

yield_raw <- function(c,f){
  b0 + bc*c + bcc*c^2 + bf*f + bff*f^2 + bcf*c*f
}

## Adjust b0 so yield matches paper steady state
b0 <- b0 + (y_ss_target - yield_raw(c_ss,f_ss))

yield_fun <- function(c,f){
  pmax(1e-9, b0 + bc*c + bcc*c^2 + bf*f + bff*f^2 + bcf*c*f)
}

## ----------------------------
## 5) SOC dynamic equation
## ----------------------------
c_next <- function(c,f,a){
  inflow <- (a * F_C * k_rg * yield_fun(c,f))^Bexp
  pmax(c - Dmin*c + Ain*inflow, 0)
}

## ----------------------------
## 6) Profit function WITH CARBON SUBSIDY
## ----------------------------
profit_fun <- function(c,f,a,c2,y,pC){
  ## CORRECT revenue (paper): p*y + q*R*(1-a)
  base <- p_maize*y +
    qR * k_rg * y * (1 - a) -
    nP*f - m_fix
  
  ## Carbon payment
  carbon <- pC * (c2 - c)
  
  base + carbon
}

## ----------------------------
## 7) Simulator
## ----------------------------
simulate_path <- function(f_vec,a_vec,c0,pC){
  c <- numeric(T_h+1); c[1] <- c0
  y <- numeric(T_h)
  pr <- numeric(T_h)
  
  for(t in 1:T_h){
    y[t] <- yield_fun(c[t],f_vec[t])
    c2   <- c_next(c[t],f_vec[t],a_vec[t])
    pr[t] <- profit_fun(c[t],f_vec[t],a_vec[t],c2,y[t],pC)
    c[t+1] <- c2
  }
  
  list(c=c,y=y,pr=pr)
}

sum_profits <- function(sim) sum(disc * sim$pr)

## ----------------------------
## 8) Control mapping
## ----------------------------
map_controls <- function(x){
  f_raw <- x[1:T_h]
  a_raw <- x[(T_h+1):(2*T_h)]
  f <- pmin(0.24,pmax(0,f_raw))
  a <- 1/(1+exp(-a_raw))
  list(f=f,a=a)
}

## ----------------------------
## 9) EXTENSION Objective
##    ONLY SOC Balance Penalty:
##    D*c_t - A*(a F k y)^B = 0   in last 10 years
## ----------------------------
lambda_s <- 20
K_tail   <- 10

objective_extension <- function(x,c0,pC){
  ctrl <- map_controls(x)
  sim  <- simulate_path(ctrl$f,ctrl$a,c0,pC)
  pv   <- sum_profits(sim)
  
  ## Stability penalty
  idx <- (T_h-K_tail+1):T_h
  g <- Dmin*sim$c[idx] -
    Ain*(ctrl$a[idx]*F_C*k_rg*yield_fun(sim$c[idx],ctrl$f[idx]))^Bexp
  
  pen_stab <- -lambda_s * sum(g^2)
  
  -(pv + pen_stab)
}

## ----------------------------
## 10) Initialization
## ----------------------------
set.seed(1)
c0 <- 19.12
x0 <- c(rep(f_ss,T_h), qlogis(rep(a_ss,T_h)))

## ----------------------------
## 11) Sweep over carbon prices
## ----------------------------
pC_grid <- c(0,50,100,150,200)

run_case <- function(pC_val){
  opt <- optim(
    par = x0,
    fn  = objective_extension,
    c0  = c0,
    pC  = pC_val,
    method = "L-BFGS-B",
    control = list(maxit = 2000)
  )
  
  ## warm start next run
  x0 <<- opt$par
  
  ctrl <- map_controls(opt$par)
  sim  <- simulate_path(ctrl$f,ctrl$a,c0,pC_val)
  
  data.frame(
    pC     = pC_val,
    f_tail = mean(ctrl$f[(T_h-9):T_h]),
    a_tail = mean(ctrl$a[(T_h-9):T_h]),
    c_T    = sim$c[T_h+1],
    y_T    = sim$y[T_h],
    PV     = sum_profits(sim)
  )
}

results <- do.call(rbind, lapply(pC_grid, run_case))
print(results)

## Add deltas vs pC=0 (baseline extension)
results$dc   <- results$c_T - results$c_T[results$pC==0]
results$dy   <- results$y_T - results$y_T[results$pC==0]
results$dPV  <- results$PV  - results$PV[results$pC==0]
results$CO2e <- results$dc * 3.67

cat("\n=== DELTA TABLE ===\n")
print(results)
