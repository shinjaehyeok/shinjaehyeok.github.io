# Exact binomial power calculation
n_to_power <- function(n, p0 = 0.5, p1 = 0.6, alpha = 0.05) {
  n <- ceiling(n)
  thres <- qbinom(alpha, n, p0, lower.tail = FALSE)
  pwr <- pbinom(thres, n, p1, lower.tail = FALSE)
  return(pwr)
}

power_to_n <- function(beta, p0 = 0.5, p1 = 0.6, alpha = 0.05) {
  f <- function(n) {
    beta - n_to_power(n, p0, p1, alpha)
  }
  root <- stats::uniroot(f, c(1, 1e+4), tol = 1e-6)
  return(ceiling(root$root))
}

# Binomial test 
run_binom_test <- function(n, p_true,
                           p0 = 0.5, p1 = 0.6,
                           alpha = 0.05) {
  thres <- qbinom(alpha, n, p0, lower.tail = FALSE)
  num_head <- rbinom(1, size = n, prob = p_true)
  return(num_head >= thres)
}

# SPRT
run_wald_sprt <- function(p_true, p0 = 0.5, p1 = 0.6,
                      alpha = 0.05, beta = 0.95, max_iter = 1e+3L) {
  a <- log(1 - beta)
  b <- log(1 / alpha)
  s <- 0
  for (n in 1:max_iter) {
    # Observe a new sample
    x <- rbinom(1, 1, p_true)
    # Update S_n
    s <- s + ifelse(x == 1, log(p1/p0), log((1-p1)/(1-p0)))
    # Make a decision
    if (s >= b) {
      decision <- "Reject the null"
      is_reject_null <- TRUE
      break 
    } else if (s <= a) {
      decision <- "Reject the alternative"
      is_reject_null <- FALSE
      break
    } 
    if (n == max_iter) {
      decision <- "Test reached the max interation"
      is_reject_null <- FALSE
    }
  } 
  return(list(
    stopped_n = n,
    decision = decision,
    is_reject_null = is_reject_null
  ))
}

# Binomial vs SPRT
alpha <- 0.05
p0 <- 0.5
p1 <- 0.6
beta <- 0.95
n <- power_to_n(beta, p0, p1, alpha) # n = 280
max_iter <- 1e+4L

# Under null p_true = 0.5
binom_null_err <- mean(replicate(max_iter, {
  run_binom_test(n, p_true = p0,
                 p0, p1,
                 alpha)
}))

sprt_null_out <- replicate(max_iter, {
  out <- run_wald_sprt(p_true = p0, p0, p1,
                alpha, beta, max_iter = 1e+3L)
  return(c(stopped_n = out$stopped_n, is_reject_null = out$is_reject_null))
})
sprt_null_err <- mean(sprt_null_out["is_reject_null", ])
sprt_null_sample_size <- mean(sprt_null_out["stopped_n", ])

# Under alternative p_true = 0.6
binom_power <- mean(replicate(max_iter, {
  run_binom_test(n, p_true = p1,
                 p0, p1,
                 alpha)
}))

sprt_power_out <- replicate(max_iter, {
  out <- run_wald_sprt(p_true = p1, p0, p1,
                       alpha, beta, max_iter = 1e+3L)
  return(c(stopped_n = out$stopped_n, is_reject_null = out$is_reject_null))
})
sprt_power_err <- mean(sprt_power_out["is_reject_null", ])
sprt_power_sample_size <- mean(sprt_power_out["stopped_n", ])

# Mis-specified 1. small p_true = 0.55
binom_power <- mean(replicate(max_iter, {
  run_binom_test(n, p_true = 0.55,
                 p0, p1,
                 alpha)
}))

sprt_power_out <- replicate(max_iter, {
  out <- run_wald_sprt(p_true = 0.55, p0, p1,
                       alpha, beta, max_iter = 1e+3L)
  return(c(stopped_n = out$stopped_n, is_reject_null = out$is_reject_null))
})
sprt_power_err <- mean(sprt_power_out["is_reject_null", ])
sprt_power_sample_size <- mean(sprt_power_out["stopped_n", ])

# Mis-specified 2. small p_true = 0.65
binom_power <- mean(replicate(max_iter, {
  run_binom_test(n, p_true = 0.65,
                 p0, p1,
                 alpha)
}))

sprt_power_out <- replicate(max_iter, {
  out <- run_wald_sprt(p_true = 0.65, p0, p1,
                       alpha, beta, max_iter = 1e+3L)
  return(c(stopped_n = out$stopped_n, is_reject_null = out$is_reject_null))
})
sprt_power_err <- mean(sprt_power_out["is_reject_null", ])
sprt_power_sample_size <- mean(sprt_power_out["stopped_n", ])



# If STCP is not installed then run below commands
# install.packages("devtools")
# devtools::install_github("shinjaehyeok/STCP")

library(STCP)


max_sample <- 2*n
min_sample <- n / 10
psi_fn_list <- generate_sub_B_fn(p = p0)


baseline_out <- compute_baseline_for_sample_size(alpha,
                                 n_upper = max_sample,
                                 n_lower = min_sample,
                                 psi_fn_list = psi_fn_list,
                                 v_min = 1,
                                 k_max = 200,
                                 tol = 1e-6)

# Compute optimal delta star
delta_star <- p1 - p0
delta_upper <- baseline_out$delta_upper
delta_lower <- baseline_out$delta_lower

v_min <- 1
k_max <- 1e+3

# Build CP detectors
# When delta_lower = delta_upper = delta_star
stcp_star <- build_stcp_exp(
  alpha,
  p0,
  delta_star,
  delta_star,
  is_test = TRUE,
  psi_fn_list,
  s_fn = function(x) {
    x - p0
  },
  v_fn = function(x) {
    1
  },
  v_min,
  k_max,
  tol = 1e-6
)

# When delta_lower < delta_star < delta_upper
stcp_mix <- build_stcp_exp(
  alpha,
  p0,
  delta_lower,
  delta_upper,
  is_test = TRUE,
  psi_fn_list
)


max_simul_num <- 1e+3L
reject_ind <- function(p){
  x_vec <- rbinom(max_sample, 1, p)
  run_star <- run_stcp(x_vec, stcp_mix)
  return(run_star$stopped_ind)
}

# type 1 error
simul_under_null <- sapply(rep(p0, max_simul_num), reject_ind)
type_1_err <- mean(ifelse(is.infinite(simul_under_null), 0, 1))

# power
simul_under_alter <- sapply(rep(p1, max_simul_num), reject_ind)
pwr_sprt <- mean(ifelse(is.infinite(simul_under_alter), 0, 1))
avg_sample_size <- mean(ifelse(is.infinite(simul_under_alter), max_sample, simul_under_alter))

# Fixed
thres <- qbinom(alpha, ceiling(avg_sample_size), p0, lower.tail = FALSE)
pwr <- pbinom(thres, ceiling(avg_sample_size), p1, lower.tail = FALSE)





run_mix <- run_stcp(x_vec, stcp_mix)

plot(run_star)
lines(1:max_sample, run_mix$log_mix_e_vec, col = 2)
print(run_star)
print(run_mix)