n <- 200
alpha <- 0.05
p0 <- 0.5
p1 <- 0.6

thres <- qbinom(alpha, n, p0, lower.tail = FALSE)
pwr <- pbinom(thres, n, p1, lower.tail = FALSE)



# If STCP is not installed then run below commands
# install.packages("devtools")
# devtools::install_github("shinjaehyeok/STCP")

library(STCP)

# Bernoulli case ----
# Null : Ber(0.5)
# Alternative: Ber(0.6)

max_sample <- 2*n

# Compute optimal delta star
delta_star <- p1 - p0
delta_upper <- 0.2
delta_lower <- 0.01

psi_fn_list <- generate_sub_B_fn(p = p0)
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


max_simul_num <- 1e+3L
reject_ind <- function(p){
  x_vec <- rbinom(max_sample, 1, p)
  run_star <- run_stcp(x_vec, stcp_star)
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