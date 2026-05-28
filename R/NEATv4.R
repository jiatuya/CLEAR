# ============================================================
# NEAT full Bayes with Rcpp
# Gaussian-mixture alternative
#
# Modification:
#   In p-value mode, obs = -log10(p), so support is x >= 0.
#   The Gaussian-mixture f1 is now truncated at 0 and renormalized:
#
#     f1_trunc(x) = f1_raw(x) / P_raw(X >= 0),  x >= 0
#
#   where
#
#     P_raw(X >= 0) = sum_k w_k * Phi(mu_k / sigma_k)
#
# Outputs traces for:
#   log posterior
#   covered gene proportion
#   p
#   mixture weights
#   mixture mus
#   mixture sigmas
#
# 04/30/26
# ============================================================

suppressPackageStartupMessages({
  library(Rcpp)
})

Rcpp::sourceCpp(code = '
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix cpp_component_logdens(NumericVector x,
                                    NumericVector w,
                                    NumericVector mu,
                                    NumericVector sigma2) {
  int n = x.size();
  int K = w.size();
  NumericMatrix out(n, K);

  for (int k = 0; k < K; k++) {
    double logw = std::log(std::max(w[k], 1e-300));
    double s2 = std::max(sigma2[k], 1e-12);
    double sd = std::sqrt(s2);
    double logconst = -0.5 * std::log(2.0 * M_PI) - std::log(sd);

    for (int i = 0; i < n; i++) {
      double z = (x[i] - mu[k]) / sd;
      out(i, k) = logw + logconst - 0.5 * z * z;
    }
  }

  return out;
}

// [[Rcpp::export]]
NumericVector cpp_rowLogSumExps(NumericMatrix mat) {
  int n = mat.nrow();
  int K = mat.ncol();
  NumericVector out(n);

  for (int i = 0; i < n; i++) {
    double m = mat(i, 0);

    for (int k = 1; k < K; k++) {
      if (mat(i, k) > m) {
        m = mat(i, k);
      }
    }

    double s = 0.0;

    for (int k = 0; k < K; k++) {
      s += std::exp(mat(i, k) - m);
    }

    out[i] = m + std::log(s);
  }

  return out;
}

// [[Rcpp::export]]
double cpp_delta_ll(IntegerVector old_states,
                    IntegerVector new_states,
                    NumericVector log_lr_obs) {
  int n = old_states.size();
  double out = 0.0;

  for (int i = 0; i < n; i++) {
    out += (new_states[i] - old_states[i]) * log_lr_obs[i];
  }

  return out;
}

// [[Rcpp::export]]
IntegerVector cpp_sample_allocations(NumericMatrix log_prob) {
  RNGScope scope;

  int n = log_prob.nrow();
  int K = log_prob.ncol();
  IntegerVector z(n);

  for (int i = 0; i < n; i++) {
    double m = log_prob(i, 0);

    for (int k = 1; k < K; k++) {
      if (log_prob(i, k) > m) {
        m = log_prob(i, k);
      }
    }

    NumericVector p(K);
    double s = 0.0;

    for (int k = 0; k < K; k++) {
      p[k] = std::exp(log_prob(i, k) - m);
      s += p[k];
    }

    for (int k = 0; k < K; k++) {
      p[k] /= s;
    }

    double u = R::runif(0.0, 1.0);
    double csum = 0.0;
    int draw = K - 1;

    for (int k = 0; k < K; k++) {
      csum += p[k];

      if (u <= csum) {
        draw = k;
        break;
      }
    }

    z[i] = draw + 1;
  }

  return z;
}
')


# ============================================================
# Main NEAT function
# ============================================================

NEAT <- function(
    genes,
    stats,
    GO,
    n_iterations = 1000000,
    burn_in = 500000,
    record_likelihoods = 1000,
    random_initial = FALSE,
    initial_set_states = NULL,
    data_type = c("p-values", "test_stats"),
    null = c("normal", "t"),
    t_df = Inf,
    setup = NULL,
    
    # MGSA-style p prior
    p_grid_size = 20,
    initial_p_idx = NULL,
    sample_p_every = 1,
    
    # Gaussian mixture prior/settings
    n_components = 3,
    dirichlet_alpha = 1,
    mu0 = NULL,
    kappa0 = 0.01,
    sigma2_a0 = 2,
    sigma2_b0 = NULL,
    mixture_update_every = 5,
    initialize_from_all_genes = TRUE,
    
    verbose = TRUE
) {
  data_type <- match.arg(data_type)
  null <- match.arg(null)
  
  if (is.null(setup)) {
    setup <- NEAT_build_fullBayes_setup(
      genes = genes,
      stats = stats,
      GO = GO,
      data_type = data_type,
      null = null,
      t_df = t_df
    )
  }
  
  go <- setup$go_index
  n_genes <- setup$n_genes
  m_sets  <- setup$m_sets
  
  if (m_sets < 2L) {
    stop("At least two gene sets are required to construct a valid p grid in (0, 1).")
  }
  
  # ------------------------------------------------------------
  # 1. Build MGSA-style discrete p grid
  # ------------------------------------------------------------
  p_grid <- NEAT_build_p_grid(
    m_sets = m_sets,
    p_grid_size = p_grid_size
  )
  
  n_p_grid <- length(p_grid)
  
  if (is.null(initial_p_idx)) {
    current_p_idx <- n_p_grid
  } else {
    current_p_idx <- as.integer(initial_p_idx)
    
    if (current_p_idx < 1L || current_p_idx > n_p_grid) {
      stop("initial_p_idx must be between 1 and length(p_grid).")
    }
  }
  
  current_p <- p_grid[current_p_idx]
  
  # Uniform prior over the discrete p grid.
  log_p_prior_weights <- rep(-log(n_p_grid), n_p_grid)
  
  # ------------------------------------------------------------
  # 2. Initialize gene-set states
  # ------------------------------------------------------------
  if (!is.null(initial_set_states)) {
    current_set_states <- as.integer(initial_set_states)
    
    if (length(current_set_states) != m_sets) {
      stop("initial_set_states must have length equal to number of gene sets.")
    }
    
    current_set_states[current_set_states != 0L] <- 1L
    
  } else if (isTRUE(random_initial)) {
    current_set_states <- sample(0:1, m_sets, replace = TRUE)
    
  } else {
    current_set_states <- rep(0L, m_sets)
  }
  
  # ------------------------------------------------------------
  # 3. Build gene coverage counts
  # ------------------------------------------------------------
  coverage_counts <- integer(n_genes)
  on_sets <- which(current_set_states == 1L)
  
  if (length(on_sets) > 0L) {
    for (j in on_sets) {
      idxs <- go[[j]]
      
      if (length(idxs) > 0L) {
        coverage_counts[idxs] <- coverage_counts[idxs] + 1L
      }
    }
  }
  
  current_gene_states <- as.integer(coverage_counts > 0L)
  current_n_active_genes <- sum(current_gene_states)
  current_n_on_sets <- sum(current_set_states)
  
  # ------------------------------------------------------------
  # 4. Initialize Gaussian mixture alternative f1
  # ------------------------------------------------------------
  theta <- NEAT_initialize_mixture(
    setup = setup,
    gene_states = current_gene_states,
    n_components = n_components,
    dirichlet_alpha = dirichlet_alpha,
    mu0 = mu0,
    kappa0 = kappa0,
    sigma2_a0 = sigma2_a0,
    sigma2_b0 = sigma2_b0,
    initialize_from_all_genes = initialize_from_all_genes
  )
  
  # ------------------------------------------------------------
  # 5. Helper functions
  # ------------------------------------------------------------
  
  compute_log_prior_T_from_n_on <- function(n_on, p) {
    n_on * log(p) + (m_sets - n_on) * log1p(-p)
  }
  
  sample_p_given_n_on <- function(n_on) {
    log_prob <- n_on * log(p_grid) +
      (m_sets - n_on) * log1p(-p_grid) +
      log_p_prior_weights
    
    idx <- sample_from_logprob(log_prob)
    
    list(
      idx = idx,
      p = p_grid[idx],
      log_prior_T = compute_log_prior_T_from_n_on(n_on, p_grid[idx])
    )
  }
  
  refresh_mixture_cache <- function(x, theta, log_f0_obs, truncate_f1_at_zero = FALSE) {
    log_comp_mat <- cpp_component_logdens(
      x = x,
      w = theta$w,
      mu = theta$mu,
      sigma2 = theta$sigma2
    )
    
    log_f1_obs_raw <- cpp_rowLogSumExps(log_comp_mat)
    
    # ----------------------------------------------------------
    # Important modification:
    #
    # In p-value mode:
    #   x = -log10(p), so x has support [0, Inf).
    #
    # The raw Gaussian mixture puts some density below 0.
    # We truncate that negative-support mass and renormalize:
    #
    #   f1_trunc(x) = f1_raw(x) / P_raw(X >= 0), x >= 0.
    #
    # For component k:
    #   P_k(X >= 0) = Phi(mu_k / sigma_k).
    #
    # For the mixture:
    #   P_mix(X >= 0) = sum_k w_k * Phi(mu_k / sigma_k).
    # ----------------------------------------------------------
    if (isTRUE(truncate_f1_at_zero)) {
      positive_mass <- NEAT_positive_mass_mixture(theta)
      
      if (!is.finite(positive_mass) || positive_mass <= 0) {
        stop("Invalid positive-domain mass for truncated f1.")
      }
      
      log_f1_obs <- log_f1_obs_raw - log(positive_mass)
      log_f1_obs[x < 0] <- -Inf
      
    } else {
      positive_mass <- NA_real_
      log_f1_obs <- log_f1_obs_raw
    }
    
    log_lr_obs <- log_f1_obs - log_f0_obs
    log_prior_theta <- NEAT_log_prior_theta(theta)
    
    list(
      log_comp_mat = log_comp_mat,
      log_f1_obs_raw = log_f1_obs_raw,
      log_f1_obs = log_f1_obs,
      log_lr_obs = log_lr_obs,
      log_prior_theta = log_prior_theta,
      f1_positive_mass = positive_mass
    )
  }
  
  # ------------------------------------------------------------
  # 6. Initial posterior components
  # ------------------------------------------------------------
  mix_cache <- refresh_mixture_cache(
    x = setup$obs,
    theta = theta,
    log_f0_obs = setup$log_f0_obs,
    truncate_f1_at_zero = setup$truncate_f1_at_zero
  )
  
  current_log_lr_obs <- mix_cache$log_lr_obs
  current_prior_theta <- mix_cache$log_prior_theta
  current_f1_positive_mass <- mix_cache$f1_positive_mass
  
  current_prior_T <- compute_log_prior_T_from_n_on(
    n_on = current_n_on_sets,
    p = current_p
  )
  
  current_ll <- sum(current_gene_states * current_log_lr_obs)
  current_pi_hat <- current_n_active_genes / n_genes
  
  current_lp <- current_prior_T + current_prior_theta + current_ll
  
  # ------------------------------------------------------------
  # 7. Allocate traces
  # ------------------------------------------------------------
  n_post <- max(0L, n_iterations - burn_in)
  n_records <- min(record_likelihoods, max(1L, n_post))
  
  likelihood_interval <- if (n_records > 0L) {
    floor(n_post / n_records)
  } else {
    0L
  }
  
  if (likelihood_interval < 1L && n_records > 0L) {
    likelihood_interval <- 1L
  }
  
  log_posterior_trace <- numeric(n_records)
  covered_gene_prop_trace <- rep(NA_real_, n_records)
  p_trace <- rep(NA_real_, n_records)
  
  weight_trace <- matrix(NA_real_, nrow = n_records, ncol = theta$K)
  mu_trace     <- matrix(NA_real_, nrow = n_records, ncol = theta$K)
  sigma_trace  <- matrix(NA_real_, nrow = n_records, ncol = theta$K)
  
  colnames(weight_trace) <- paste0("w", seq_len(theta$K))
  colnames(mu_trace)     <- paste0("mu", seq_len(theta$K))
  colnames(sigma_trace)  <- paste0("sigma", seq_len(theta$K))
  
  # Useful diagnostic for p-value mode.
  f1_positive_mass_trace <- rep(NA_real_, n_records)
  
  on_frequency <- numeric(m_sets)
  recording_index <- 1L
  
  flip01 <- integer(m_sets)
  flip10 <- integer(m_sets)
  
  # ============================================================
  # 8. MCMC
  # ============================================================
  
  for (iter in seq_len(n_iterations)) {
    if (isTRUE(verbose) && iter %% 5000L == 0L) {
      message(
        "iter = ", iter,
        " | active_sets = ", current_n_on_sets,
        " | covered_gene_prop = ", signif(current_pi_hat, 4),
        " | p = ", signif(current_p, 4)
      )
    }
    
    # ----------------------------------------------------------
    # 8.1 Propose flipping one gene set
    # ----------------------------------------------------------
    idx <- sample.int(m_sets, 1L)
    g_idx <- go[[idx]]
    
    proposed_n_on_sets <- current_n_on_sets +
      if (current_set_states[idx] == 0L) 1L else -1L
    
    proposed_prior_T <- compute_log_prior_T_from_n_on(
      n_on = proposed_n_on_sets,
      p = current_p
    )
    
    delta_prior_T <- proposed_prior_T - current_prior_T
    
    if (length(g_idx) > 0L) {
      proposed_coverage_local <- coverage_counts[g_idx] +
        if (current_set_states[idx] == 0L) 1L else -1L
      
      old_gene_states_local <- current_gene_states[g_idx]
      new_gene_states_local <- as.integer(proposed_coverage_local > 0L)
      
      delta_ll <- cpp_delta_ll(
        old_states = old_gene_states_local,
        new_states = new_gene_states_local,
        log_lr_obs = current_log_lr_obs[g_idx]
      )
      
      delta_lp <- delta_prior_T + delta_ll
      
      if (log(stats::runif(1)) < delta_lp) {
        old_val <- current_set_states[idx]
        
        current_set_states[idx] <- 1L - current_set_states[idx]
        coverage_counts[g_idx] <- proposed_coverage_local
        current_gene_states[g_idx] <- new_gene_states_local
        
        current_n_active_genes <- current_n_active_genes +
          sum(new_gene_states_local - old_gene_states_local)
        
        current_n_on_sets <- proposed_n_on_sets
        
        current_prior_T <- proposed_prior_T
        current_ll <- current_ll + delta_ll
        current_pi_hat <- current_n_active_genes / n_genes
        
        current_lp <- current_prior_T + current_prior_theta + current_ll
        
        if (iter > burn_in) {
          if (old_val == 0L) {
            flip01[idx] <- flip01[idx] + 1L
          } else {
            flip10[idx] <- flip10[idx] + 1L
          }
        }
      }
      
    } else {
      # Empty gene set: only prior changes.
      delta_lp <- delta_prior_T
      
      if (log(stats::runif(1)) < delta_lp) {
        old_val <- current_set_states[idx]
        
        current_set_states[idx] <- 1L - current_set_states[idx]
        current_n_on_sets <- proposed_n_on_sets
        
        current_prior_T <- proposed_prior_T
        current_lp <- current_prior_T + current_prior_theta + current_ll
        
        if (iter > burn_in) {
          if (old_val == 0L) {
            flip01[idx] <- flip01[idx] + 1L
          } else {
            flip10[idx] <- flip10[idx] + 1L
          }
        }
      }
    }
    
    # ----------------------------------------------------------
    # 8.2 Gibbs sample p from p | T
    # ----------------------------------------------------------
    if (sample_p_every >= 1L && iter %% sample_p_every == 0L) {
      p_draw <- sample_p_given_n_on(current_n_on_sets)
      
      current_p_idx <- p_draw$idx
      current_p <- p_draw$p
      current_prior_T <- p_draw$log_prior_T
      
      current_lp <- current_prior_T + current_prior_theta + current_ll
    }
    
    # ----------------------------------------------------------
    # 8.3 Update Gaussian mixture f1 given currently active genes
    # ----------------------------------------------------------
    if (mixture_update_every >= 1L && iter %% mixture_update_every == 0L) {
      theta <- NEAT_sample_mixture_given_active_rcpp(
        x = setup$obs,
        gene_states = current_gene_states,
        theta = theta
      )
      
      mix_cache <- refresh_mixture_cache(
        x = setup$obs,
        theta = theta,
        log_f0_obs = setup$log_f0_obs,
        truncate_f1_at_zero = setup$truncate_f1_at_zero
      )
      
      current_log_lr_obs <- mix_cache$log_lr_obs
      current_prior_theta <- mix_cache$log_prior_theta
      current_f1_positive_mass <- mix_cache$f1_positive_mass
      
      current_ll <- sum(current_gene_states * current_log_lr_obs)
      current_pi_hat <- current_n_active_genes / n_genes
      
      current_lp <- current_prior_T + current_prior_theta + current_ll
    }
    
    # ----------------------------------------------------------
    # 8.4 Record traces after burn-in
    # ----------------------------------------------------------
    if (iter > burn_in && n_records > 0L) {
      since_burn <- iter - burn_in
      
      if (since_burn %% likelihood_interval == 0L &&
          recording_index <= n_records) {
        
        log_posterior_trace[recording_index] <- current_lp
        covered_gene_prop_trace[recording_index] <- current_pi_hat
        p_trace[recording_index] <- current_p
        
        weight_trace[recording_index, ] <- theta$w
        mu_trace[recording_index, ] <- theta$mu
        sigma_trace[recording_index, ] <- sqrt(theta$sigma2)
        
        f1_positive_mass_trace[recording_index] <- current_f1_positive_mass
        
        recording_index <- recording_index + 1L
      }
      
      on_frequency <- on_frequency + current_set_states
    }
  }
  
  # ============================================================
  # 9. Summaries
  # ============================================================
  
  T_steps <- max(0L, n_iterations - burn_in)
  denom <- max(1L, T_steps)
  
  on_counts <- on_frequency
  on_frequency <- on_counts / denom
  names(on_frequency) <- names(go)
  
  N01 <- flip01
  N10 <- flip10
  N11 <- pmax(0L, as.integer(round(on_counts - N10)))
  
  off_counts <- T_steps - on_counts
  N00 <- pmax(0L, as.integer(round(off_counts - N01)))
  
  transition_counts <- cbind(
    `1->0` = N10,
    `1->1` = N11,
    `0->0` = N00,
    `0->1` = N01
  )
  
  rownames(transition_counts) <- names(go)
  
  transition_freq <- if (T_steps > 0L) {
    transition_counts / T_steps
  } else {
    transition_counts
  }
  
  # Trim unused trace rows if recording did not exactly fill n_records.
  n_filled <- recording_index - 1L
  
  if (n_filled < n_records) {
    keep_idx <- seq_len(max(0L, n_filled))
    
    log_posterior_trace <- log_posterior_trace[keep_idx]
    covered_gene_prop_trace <- covered_gene_prop_trace[keep_idx]
    p_trace <- p_trace[keep_idx]
    
    weight_trace <- weight_trace[keep_idx, , drop = FALSE]
    mu_trace <- mu_trace[keep_idx, , drop = FALSE]
    sigma_trace <- sigma_trace[keep_idx, , drop = FALSE]
    
    f1_positive_mass_trace <- f1_positive_mass_trace[keep_idx]
  }
  
  list(
    on_frequency = on_frequency,
    transition_counts = transition_counts,
    transition_freq = transition_freq,
    
    log_posterior_trace = log_posterior_trace,
    covered_gene_prop_trace = covered_gene_prop_trace,
    
    p_trace = p_trace,
    weight_trace = weight_trace,
    mu_trace = mu_trace,
    sigma_trace = sigma_trace,
    
    # Diagnostic:
    # In p-value mode, this is P_raw(X >= 0) for the current raw Gaussian mixture.
    # In test-stat mode, this is NA.
    f1_positive_mass_trace = f1_positive_mass_trace,
    
    final_p = current_p,
    final_p_idx = current_p_idx,
    p_grid = p_grid,
    
    final_theta = theta,
    final_set_states = current_set_states,
    final_gene_states = current_gene_states,
    
    data_type = setup$data_type,
    null = setup$null,
    t_df = setup$t_df,
    setup = setup
  )
}


# ============================================================
# Setup builder
# ============================================================

NEAT_build_fullBayes_setup <- function(
    genes,
    stats,
    GO,
    data_type = c("p-values", "test_stats"),
    null = c("normal", "t"),
    t_df = Inf
) {
  data_type <- match.arg(data_type)
  null <- match.arg(null)
  
  genes <- as.character(genes)
  stats <- as.numeric(stats)
  
  keep <- is.finite(stats) & !is.na(genes)
  genes <- genes[keep]
  stats <- stats[keep]
  
  if (data_type == "p-values") {
    p_values <- stats
    
    p_values[p_values < 0] <- 0
    p_values[p_values > 1] <- 1
    p_values[p_values == 0] <- .Machine$double.xmin
    
    test_stats <- NULL
    
    obs <- p_to_x10(p_values)
    
    # Null density for x = -log10(p), x >= 0.
    f0_obs <- f0_exp_on_x10(obs)
    
    obs_scale <- "x10"
    
    # Important:
    # f1 should also respect x >= 0 in p-value mode.
    truncate_f1_at_zero <- TRUE
    
  } else {
    test_stats <- stats
    p_values <- NULL
    
    obs <- test_stats
    
    f0_obs <- if (null == "normal" || is.infinite(t_df)) {
      stats::dnorm(obs)
    } else {
      stats::dt(obs, df = t_df)
    }
    
    obs_scale <- "test_stats"
    
    # Test statistics can be negative, so do not truncate f1.
    truncate_f1_at_zero <- FALSE
  }
  
  f0_obs <- pmax(f0_obs, 1e-300)
  log_f0_obs <- log(f0_obs)
  
  gene_index <- setNames(seq_along(genes), genes)
  
  go <- lapply(GO, function(x) {
    idx <- unname(gene_index[as.character(x)])
    unique(idx[!is.na(idx)])
  })
  
  names(go) <- names(GO)
  
  structure(
    list(
      genes = genes,
      GO = GO,
      go_index = go,
      n_genes = length(genes),
      m_sets = length(go),
      stats = stats,
      p_values = p_values,
      test_stats = test_stats,
      data_type = data_type,
      null = null,
      t_df = t_df,
      obs = obs,
      obs_scale = obs_scale,
      f0_obs = f0_obs,
      log_f0_obs = log_f0_obs,
      truncate_f1_at_zero = truncate_f1_at_zero
    ),
    class = "neat_fullBayes_setup"
  )
}


# ============================================================
# Mixture initialization
# ============================================================

NEAT_initialize_mixture <- function(
    setup,
    gene_states,
    n_components = 3,
    dirichlet_alpha = 1,
    mu0 = NULL,
    kappa0 = 0.01,
    sigma2_a0 = 2,
    sigma2_b0 = NULL,
    initialize_from_all_genes = TRUE
) {
  x_all <- setup$obs
  x_active <- x_all[gene_states == 1L]
  
  x_init <- if (length(x_active) >= n_components) {
    x_active
  } else if (isTRUE(initialize_from_all_genes)) {
    x_all
  } else {
    x_active
  }
  
  if (length(x_init) == 0L) {
    x_init <- x_all
  }
  
  if (is.null(mu0)) {
    mu0 <- mean(x_all)
  }
  
  if (is.null(sigma2_b0)) {
    vx <- stats::var(x_all)
    
    if (!is.finite(vx) || vx <= 0) {
      vx <- 1
    }
    
    sigma2_b0 <- vx
  }
  
  mu <- as.numeric(stats::quantile(
    x_init,
    probs = seq(0.2, 0.8, length.out = n_components),
    na.rm = TRUE,
    names = FALSE
  ))
  
  vx <- stats::var(x_init)
  
  if (!is.finite(vx) || vx <= 0) {
    vx <- 1
  }
  
  list(
    K = as.integer(n_components),
    w = rep(1 / n_components, n_components),
    mu = mu,
    sigma2 = rep(vx, n_components),
    alpha = rep(dirichlet_alpha, n_components),
    mu0 = mu0,
    kappa0 = kappa0,
    sigma2_a0 = sigma2_a0,
    sigma2_b0 = sigma2_b0
  )
}


# ============================================================
# Mixture prior density
# ============================================================

NEAT_log_prior_theta <- function(theta) {
  w <- pmax(theta$w, 1e-300)
  sigma2 <- pmax(theta$sigma2, 1e-12)
  
  logp_w <- lgamma(sum(theta$alpha)) -
    sum(lgamma(theta$alpha)) +
    sum((theta$alpha - 1) * log(w))
  
  logp_rest <- 0
  
  for (k in seq_len(theta$K)) {
    logp_sigma2 <- dinvgamma_log(
      sigma2[k],
      shape = theta$sigma2_a0,
      rate = theta$sigma2_b0
    )
    
    logp_mu <- stats::dnorm(
      theta$mu[k],
      mean = theta$mu0,
      sd = sqrt(sigma2[k] / theta$kappa0),
      log = TRUE
    )
    
    logp_rest <- logp_rest + logp_sigma2 + logp_mu
  }
  
  logp_w + logp_rest
}


# ============================================================
# Mixture update given active genes
# ============================================================

NEAT_sample_mixture_given_active_rcpp <- function(x, gene_states, theta) {
  x_active <- x[gene_states == 1L]
  K <- theta$K
  
  if (length(x_active) == 0L) {
    theta$w <- as.numeric(rdirichlet1(theta$alpha))
    
    for (k in seq_len(K)) {
      theta$sigma2[k] <- rinvgamma1(
        shape = theta$sigma2_a0,
        rate = theta$sigma2_b0
      )
      
      theta$mu[k] <- stats::rnorm(
        1L,
        mean = theta$mu0,
        sd = sqrt(theta$sigma2[k] / theta$kappa0)
      )
    }
    
    return(theta)
  }
  
  log_prob <- cpp_component_logdens(
    x = x_active,
    w = theta$w,
    mu = theta$mu,
    sigma2 = theta$sigma2
  )
  
  z <- cpp_sample_allocations(log_prob)
  nk <- tabulate(z, nbins = K)
  
  theta$w <- as.numeric(rdirichlet1(theta$alpha + nk))
  
  for (k in seq_len(K)) {
    xk <- x_active[z == k]
    nk_k <- length(xk)
    
    if (nk_k == 0L) {
      a_n <- theta$sigma2_a0
      b_n <- theta$sigma2_b0
      kappa_n <- theta$kappa0
      mu_n <- theta$mu0
      
    } else {
      xbar <- mean(xk)
      ss <- sum((xk - xbar)^2)
      
      kappa_n <- theta$kappa0 + nk_k
      
      mu_n <- (
        theta$kappa0 * theta$mu0 +
          nk_k * xbar
      ) / kappa_n
      
      a_n <- theta$sigma2_a0 + nk_k / 2
      
      b_n <- theta$sigma2_b0 +
        0.5 * ss +
        0.5 * (theta$kappa0 * nk_k / kappa_n) *
        (xbar - theta$mu0)^2
    }
    
    theta$sigma2[k] <- rinvgamma1(
      shape = a_n,
      rate = b_n
    )
    
    theta$mu[k] <- stats::rnorm(
      1L,
      mean = mu_n,
      sd = sqrt(theta$sigma2[k] / kappa_n)
    )
  }
  
  theta
}


# ============================================================
# Helper functions
# ============================================================

NEAT_build_p_grid <- function(m_sets, p_grid_size = 20) {
  if (m_sets < 2L) {
    stop("m_sets must be at least 2.")
  }
  
  k_max <- min(as.integer(p_grid_size), m_sets - 1L)
  
  if (k_max < 1L) {
    stop("No valid p values in (0, 1).")
  }
  
  seq_len(k_max) / m_sets
}


sample_from_logprob <- function(log_prob) {
  log_prob <- as.numeric(log_prob)
  
  if (all(!is.finite(log_prob))) {
    stop("All log probabilities are non-finite.")
  }
  
  m <- max(log_prob, na.rm = TRUE)
  prob <- exp(log_prob - m)
  prob[!is.finite(prob)] <- 0
  
  s <- sum(prob)
  
  if (!is.finite(s) || s <= 0) {
    stop("Invalid probability vector.")
  }
  
  prob <- prob / s
  
  sample.int(length(prob), size = 1L, prob = prob)
}


p_to_x10 <- function(p) {
  -log10(p)
}


f0_exp_on_x10 <- function(x) {
  lam <- log(10)
  lam * exp(-lam * pmax(x, 0))
}


# ------------------------------------------------------------
# New helper:
# positive-domain mass of the raw Gaussian mixture.
#
# For each component:
#   P_k(X >= 0) = 1 - Phi((0 - mu_k) / sigma_k)
#               = Phi(mu_k / sigma_k)
#
# For the mixture:
#   P_mix(X >= 0) = sum_k w_k * P_k(X >= 0)
# ------------------------------------------------------------
NEAT_positive_mass_mixture <- function(theta) {
  sigma <- sqrt(pmax(theta$sigma2, 1e-12))
  
  component_mass <- stats::pnorm(
    q = theta$mu / sigma,
    mean = 0,
    sd = 1
  )
  
  mass <- sum(theta$w * component_mass)
  
  pmin(pmax(mass, 1e-300), 1)
}


rinvgamma1 <- function(shape, rate) {
  1 / stats::rgamma(1L, shape = shape, rate = rate)
}


dinvgamma_log <- function(x, shape, rate) {
  if (x <= 0 || !is.finite(x)) {
    return(-Inf)
  }
  
  shape * log(rate) -
    lgamma(shape) -
    (shape + 1) * log(x) -
    rate / x
}


rdirichlet1 <- function(alpha) {
  y <- stats::rgamma(length(alpha), shape = alpha, rate = 1)
  y / sum(y)
}