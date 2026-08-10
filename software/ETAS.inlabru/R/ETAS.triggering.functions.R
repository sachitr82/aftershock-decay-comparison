################################################################################
## Functions for ETAS conditional intensity and integrated conditional intensity
################################################################################

#' Validate and standardise the temporal kernel name
#'
#' Checks that the requested temporal triggering kernel is one of the
#' kernels currently implemented in the package.
#'
#' @param kernel Character string specifying the temporal triggering kernel.
#'   One of `"ou"`, `"mse"`, or `"rate_state"`.
#'
#' @return A single character string containing the validated kernel name.

validate_temporal_kernel <- function(kernel) {
  valid.kernels <- c("ou", "mse", "rate_state")
  
  if (length(kernel) != 1 || !kernel %in% valid.kernels) {
    stop(
      "kernel must be one of: ",
      paste(valid.kernels, collapse = ", ")
    )
  }
  
  kernel
}

#' Return the parameter specification for an ETAS temporal kernel
#'
#' Constructs a table describing the physical ETAS parameters, corresponding
#' latent Gaussian component names, and link-function names required by a
#' selected temporal triggering kernel.
#'
#' The background rate \eqn{\mu}, productivity parameter \eqn{K}, and
#' magnitude productivity parameter \eqn{\alpha} are common to all kernels.
#' Kernel-specific parameters are:
#'    Omori--Utsu: \eqn{(c,p)};
#'    modified stretched exponential: \eqn{(d,\lambda,\gamma)};
#'    rate-state: \eqn{(B,t_a)}.
#' }
#'
#' @param kernel Character string specifying the temporal triggering kernel.
#'   One of `"ou"`, `"mse"`, or `"rate_state"`.
#'
#' @return A data frame with columns:
#'    physical: Physical ETAS parameter name;
#'    latent: Corresponding latent Gaussian component name used by inlabru;
#'    link: Name of the transformation from the latent Gaussian scale
#'          to the physical parameter scale. 
#'          
#' @details 
#'  `link` and `physical` match with the exception of the inherited 
#'  package conventions of physical parameter `c` being associated with
#'  link-function `c_`. The source does not document why this naming convention 
#'  was chosen but it likely avoids ambiguity with R's base `c()` function.
#' }

etas_parameter_spec <- function(kernel = "ou") {
  kernel <- validate_temporal_kernel(kernel)
  
  switch(
    kernel,
    
    ou = data.frame(
      physical = c("mu", "K", "alpha", "c", "p"),
      latent   = c("th.mu", "th.K", "th.alpha", "th.c", "th.p"),
      link     = c("mu", "K", "alpha", "c_", "p"),
      stringsAsFactors = FALSE
    ),
    
    mse = data.frame(
      physical = c("mu", "K", "alpha", "d", "lambda", "gamma"),
      latent   = c(
        "th.mu", "th.K", "th.alpha",
        "th.d", "th.lambda", "th.gamma"
      ),
      link     = c("mu", "K", "alpha", "d", "lambda", "gamma"),
      stringsAsFactors = FALSE
    ),
    
    rate_state = data.frame(
      physical = c("mu", "K", "alpha", "B", "ta"),
      latent   = c("th.mu", "th.K", "th.alpha", "th.B", "th.ta"),
      link     = c("mu", "K", "alpha", "B", "ta"),
      stringsAsFactors = FALSE
    )
  )
}

#' Validate temporal-kernel parameters
#'
#' Checks that the kernel-specific parameters required by the selected
#' temporal triggering function are present and lie within their admissible
#' parameter space.
#'
#' @param theta Named list containing the physical temporal-kernel parameters.
#' @param kernel Character string specifying the temporal triggering kernel.
#'   One of `"ou"`, `"mse"`, or `"rate_state"`.
#'
#' @return Invisibly returns \code{TRUE} when all required parameters are
#'   valid. Otherwise provides an informative error.
#'
#' @details
#' The admissible parameter spaces are
#'    Omori--Utsu: \eqn{c> 0, p > 1};
#'    modified stretched exponential: \eqn{d > 0, \lambda > 0, 0 < \gamma < 1};
#'    rate-state: \eqn{0 < B < 1, t_a > 0}.

validate_temporal_kernel_theta <- function(theta, kernel = "ou") {
  
  if (kernel == "ou") {
    if (is.null(theta$c) || is.null(theta$p)) {
      stop("OU kernel requires theta$c and theta$p.")
    }
    
    if (theta$c <= 0 || theta$p <= 1) {
      stop("OU kernel requires c > 0 and p > 1.")
    }
  }
  
  if (kernel == "mse") {
    if (is.null(theta$d) ||
        is.null(theta$lambda) ||
        is.null(theta$gamma)) {
      stop(
        "MSE kernel requires theta$d, theta$lambda and theta$gamma."
      )
    }
    
    if (theta$d <= 0 ||
        theta$lambda <= 0 ||
        theta$gamma <= 0 ||
        theta$gamma >= 1) {
      stop(
        "MSE kernel requires d > 0, lambda > 0 and 0 < gamma < 1."
      )
    }
  }
  
  if (kernel == "rate_state") {
    if (is.null(theta$B) || is.null(theta$ta)) {
      stop("Rate-state kernel requires theta$B and theta$ta.")
    }
    
    if (theta$B <= 0 || theta$B >= 1 || theta$ta <= 0) {
      stop("Rate-state kernel requires 0 < B < 1 and ta > 0.")
    }
  }
  
  invisible(TRUE)
}

#' Evaluate an ETAS temporal triggering kernel
#'
#' Evaluates one of the temporal triggering kernels implemented for the
#' temporal ETAS model.
#'
#' @param dt Numeric vector of elapsed times since a potential parent event.
#'   Values must be expressed in the same time units as the kernel parameters.
#' @param theta Named list containing the physical kernel parameters.
#' @param kernel Character string specifying the temporal triggering kernel.
#'   One of `"ou"`, `"mse"`, or `"rate_state"`.
#'
#' @return A numeric vector of temporal triggering-function values with the
#'   same length as \code{dt}. Values for negative elapsed times are zero.
#'
#' @details
#' All kernels use the parameterisation \eqn{g_k(0)=1}.
#'
#' For the Omori--Utsu kernel,
#' \deqn{g(t) = (1+t/c)^{-p}}, 
#' with \eqn{c>0} and \eqn{p>1}.
#'
#' For the modified stretched-exponential kernel,
#' \deqn{g(t) = (1+t/d)^{\gamma-1}
#'              \exp\left\{ -\lambda\left[(d+t)^\gamma-d^\gamma\right]\right\}},
#' with \eqn{d>0}, \eqn{\lambda>0}, and \eqn{0<\gamma<1}.
#'
#' For the rate-state kernel,
#' \deqn{g(t) =\frac{(1-B)\exp(-t/t_a)}{1-B\exp(-t/t_a)}},
#' with \eqn{0<B<1} and \eqn{t_a>0}.
#'
#' This function evaluates  only the temporal kernel \eqn{g_k(t)}:
#' the ETAS productivity term \eqn{K\exp\{\alpha(m-M_0)\}} is not included here. 
#' 
#' For negative or missing time lags, the kernel evaluation is 0.
#' @export
#' @examples
#' temporal_kernel(
#'   dt = c(0, 1, 2),
#'   theta = list(c = 0.1, p = 1.2),
#'   kernel = "ou"
#' )
#'
#' temporal_kernel(
#'   dt = c(0, 1, 2),
#'   theta = list(d = 0.1, lambda = 0.5, gamma = 0.5),
#'   kernel = "mse"
#' )

temporal_kernel <- function(dt, theta, kernel = "ou") {
  kernel <- validate_temporal_kernel(kernel)
  validate_temporal_kernel_theta(theta, kernel)
  
  dt <- as.numeric(dt)
  out <- numeric(length(dt))
  
  use <- !is.na(dt) & dt >= 0
  
  if (!any(use)) {
    return(out)
  }
  
  x <- dt[use]
  
  out[use] <- switch(
    kernel,
  
    ou = {(1 + x / theta$c)^(-theta$p)},
    
    mse = {
      log.scale <- (theta$gamma - 1) * log1p(x / theta$d)
      power_diff <-
        theta$d^theta$gamma *
        expm1(theta$gamma * log1p(x / theta$d))
      exp(log.scale - theta$lambda * power_diff)
    },
    
    rate_state = {
      log_g <- log1p(-theta$B) - 
        x / theta$ta -
        log1p(-theta$B * exp(-x / theta$ta))
      exp(log_g)
    }
  )
  
  out
}

#' Integrate an ETAS temporal triggering kernel over an interval
#'
#' Evaluates the analytical interval mass
#' \eqn{G_k(a,b)=\int_a^b g_k(t)\,dt} for one of the implemented temporal
#' triggering kernels.
#'
#' @param a Numeric vector of lower integration limits.
#' @param b Numeric vector of upper integration limits.
#' @param theta Named list containing the physical kernel parameters.
#' @param kernel Character string specifying the temporal triggering kernel.
#'   One of `"ou"`, `"mse"`, or `"rate_state"`.
#'
#' @return A numeric vector containing the kernel mass over each requested
#'   interval.
#'
#' @details
#' Integration intervals are interpreted relative to the parent-event time.
#' Triggering does not occur before the parent event, so negative lower or
#' upper limits are truncated at zero.
#'
#' For the Omori--Utsu kernel,
#'  \eqn{G(a,b) = \frac{c}{p-1} \left[(1+a/c)^{1-p} - (1+b/c)^{1-p} \right]}.
#'
#' For the modified stretched-exponential kernel,
#' \deqn{G(a,b) =\frac{d^{1-\gamma}}{\lambda\gamma}
#'                \left[\exp\{-\lambda[(d+a)^\gamma-d^\gamma]\}-
#'                  \exp\{-\lambda[(d+b)^\gamma-d^\gamma]\}\right]}.
#'
#' For the rate-state kernel, 
#' \deqn{G(a,b) = \frac{(1-B)t_a}{B}
#'                \left[\log\{1-Be^{-b/t_a}\} - \log\{1-Be^{-a/t_a}\}\right]}.
#'
#' These analytical interval masses are used in the existing
#' \code{ETAS.inlabru} temporal-binning strategy. The binning therefore does
#' not perform numerical quadrature when these closed-form expressions are
#' used.
#'
#' @examples
#' temporal_kernel_integral(
#'   a = 0,
#'   b = 1,
#'   theta = list(c = 0.1, p = 1.2),
#'   kernel = "ou"
#' )
#'
#' @export
temporal_kernel_integral <- function(a, b, theta, kernel = "ou") {
  
  kernel <- validate_temporal_kernel(kernel)
  validate_temporal_kernel_theta(theta, kernel)
  
  a <- as.numeric(a)
  b <- as.numeric(b)
  
  if (length(a) != length(b)) {
    stop("a and b must have the same length.")
  }
  
  a <- pmax(a, 0)
  b <- pmax(b, 0)
  
  out <- numeric(length(a))
  
  use <- !is.na(a) & !is.na(b) & b > a
  
  if (!any(use)) {
    return(out)
  }
  
  aa <- a[use]
  bb <- b[use]
  
  out[use] <- switch(
    kernel,
    
    ou = {
      lower <- (1 + aa / theta$c)^(1 - theta$p)
      upper <- (1 + bb / theta$c)^(1 - theta$p)
      
      theta$c / (theta$p - 1) * (lower - upper)
    },
    
    mse = {
      A <- theta$d^(1 - theta$gamma) / (theta$lambda * theta$gamma)
      
      mse_survival <- function(x) {
        power_diff <-
          theta$d^theta$gamma *
          expm1(theta$gamma * log1p(x / theta$d))
        
        exp(-theta$lambda * power_diff)
      }
      
      A * (mse_survival(aa) - mse_survival(bb))
    },
    
    rate_state = {
      exp_a <- exp(-aa / theta$ta)
      exp_b <- exp(-bb / theta$ta)
      
      log_ratio <- log1p(theta$B * (exp_a - exp_b) /(1 - theta$B * exp_a))
      
      (1 - theta$B) * theta$ta / theta$B *
        log_ratio
    }
  )
  
  out
}

#' Calculate the total mass of an ETAS temporal triggering kernel
#'
#' Calculates \eqn{G_k(0,\infty)=\int_0^\infty g_k(t)\,dt} 
#' for a selected temporal triggering kernel.
#'
#' @param theta Named list containing the physical kernel parameters.
#' @param kernel Character string specifying the temporal triggering kernel.
#'   One of `"ou"`, `"mse"`, or `"rate_state"`.
#'
#' @return A positive numeric scalar giving the total temporal kernel mass.
#'
#' @details
#' Under the \eqn{g_k(0)=1} parameterisation, \eqn{K} is the common initial
#' triggering amplitude. The expected number of triggered events of an
#' event of magnitude \eqn{m} is \eqn{K\exp\{\alpha(m-M_0)\}G_k(0,\infty)}.
#'
#' NB: In contrast to standard ETAS implementations, equal values of \eqn{K} 
#' do not imply equal total expected offspring across different temporal-kernel
#' families.
#'
#' @examples
#' temporal_kernel_total_mass(
#'   theta = list(c = 0.1, p = 1.2),
#'   kernel = "ou"
#' )
#'
#' @export
temporal_kernel_total_mass <- function(theta, kernel = "ou") {
  temporal_kernel_integral(a = 0, b = Inf, theta = theta, kernel = kernel)
}

#' Invert the cumulative mass of a temporal triggering kernel
#'
#' Returns the elapsed time corresponding to a specified cumulative kernel
#' mass. This function is used to generate daughter-event times by inverse
#' transform sampling.
#'
#' @param omega Numeric vector of cumulative kernel masses measured from time
#'   zero.
#' @param theta Named list containing the physical kernel parameters.
#' @param kernel Character string specifying the temporal triggering kernel.
#'   One of `"ou"`, `"mse"`, or `"rate_state"`.
#'
#' @return A numeric vector of elapsed times whose cumulative kernel masses
#'   equal the supplied values of \code{omega}.
#'
#' @details
#' Each element of \code{omega} must lie between zero and
#' \eqn{G_k(0,\infty)}. The function implements analytical inverse cumulative
#' kernel masses for all supported kernels.

temporal_kernel_cumint_inverse <- function(omega, theta, kernel = "ou") {
  kernel <- validate_temporal_kernel(kernel)
  validate_temporal_kernel_theta(theta, kernel)
  
  omega <- as.numeric(omega)
  
  total <- temporal_kernel_total_mass(theta, kernel)
  
  if (any(omega < 0) || any(omega > total)) {
    stop("omega must lie between 0 and the total kernel mass.")
  }
  
  switch(
    kernel,
    
    ou = {
      z <- 1 - omega * (theta$p - 1) / theta$c
      
      theta$c * (z^(-1 / (theta$p - 1)) -1)
    },
    
    mse = {
      A <- theta$d^(1 - theta$gamma) / (theta$lambda * theta$gamma)
      
      log_survival <- log1p(-omega / A)
      
      (theta$d^theta$gamma - log_survival / theta$lambda )^(1 / theta$gamma) -
        theta$d
    },
    
    rate_state = {
      C <- (1 - theta$B) * theta$ta / theta$B
      
      log_u <- log1p(-theta$B) + omega / C
      
      ratio <- -expm1(log_u) / theta$B
      
      -theta$ta * log(ratio)
    }
  )
}


#' Evaluate the triggering contribution from previous ETAS events
#'
#' Calculates the triggering contribution of one or more previous events to
#' the conditional intensity at a specified time.
#'
#' @param theta Named list containing physical ETAS parameters. For backwards
#'   compatibility with standard \code{ETAS.inlabru}, a numeric vector of the 
#'   form \code{c(mu, K, alpha, c, p)} is also accepted when \code{kernel = "ou"}.
#' @param t Numeric evaluation time.
#' @param th Numeric vector of previous event times.
#' @param mh Numeric vector of corresponding previous event magnitudes.
#' @param M0 Numeric magnitude threshold.
#' @param kernel Character string specifying the temporal triggering kernel.
#'   One of `"ou"`, `"mse"`, or `"rate_state"`.
#'
#' @return Numeric vector containing the triggering contribution of each
#'   previous event.
#'
#' @details
#' For parent event \eqn{i}, the contribution at time \eqn{t} is
#' \eqn{ K\exp\{\alpha(m_i-M_0)\}g_k(t-t_i)} for \eqn{t>t_i}, and zero otherwise.
#' Function retained from \code{ETAS.inlabru} but extended for other kernels
gt <- function(theta, t, th, mh, M0, kernel = "ou") {
  kernel <- validate_temporal_kernel(kernel)
  
  # Preserve the original numeric OU parameter interface.
  if (!is.list(theta)) {
  
    if (kernel != "ou") {
      stop("Numeric theta is only supported when kernel = 'ou'.")
    }
    
    if (length(theta) != 5) {
      stop("Numeric OU theta must have the form c(mu, K, alpha, c, p).")
    }
    
    theta <- list(
      mu = theta[1],
      K = theta[2],
      alpha = theta[3],
      c = theta[4],
      p = theta[5]
    )
  }
  
  if (length(th) != length(mh)) {
    stop("th and mh must have the same length.")
  }
  
  
  if (length(t) > 1 && length(th) > 1 && length(t) != length(th)) {
    stop(
      "If t and th are vectors, they must have the same length.")
  }

  t_diff <- t - th
  
  temporal_part <- temporal_kernel(t_diff, theta, kernel)
  
  output <- theta$K * exp(theta$alpha * (mh - M0)) * temporal_part

  output[!is.na(t_diff) & t_diff <= 0] <- 0
  
  output
}

#' Evaluate the temporal ETAS conditional intensity
#'
#' Calculates the conditional event intensity at a specified time given the
#' observed event history.
#'
#' @param theta Named list containing the physical ETAS parameters.
#' @param t Numeric time at which the conditional intensity is evaluated.
#' @param th Numeric vector of previous event times.
#' @param mh Numeric vector of corresponding previous event magnitudes.
#' @param M0 Numeric magnitude threshold.
#' @param kernel Character string specifying the temporal triggering kernel.
#'   One of `"ou"`, `"mse"`, or `"rate_state"`.
#'
#' @return A positive numeric scalar containing the ETAS conditional intensity.
#'
#' @details
#' The conditional intensity is
#' \deqn{\lambda(t\mid\mathcal H_t) =
#' \mu +\sum_{t_i<t} K\exp\{\alpha(m_i-M_0)\}g_k(t-t_i)}.
#'
cond_lambda <- function(theta, t, th, mh, M0) {
  if (is.null(th) || length(th) == 0 || all(th > t)) {
    return(theta$mu)
  }
  theta$mu + sum(gt(theta = theta, t = t, th = th, mh = mh, M0 = M0, 
                    kernel = kernel))
}


#' Calculate an integrated triggered ETAS contribution
#'
#' Calculates the logarithm of the expected number of directly triggered
#' events contributed by a parent event over a specified time interval.
#'
#' @param theta Named list containing the physical ETAS parameters.
#' @param th Numeric vector of parent-event times.
#' @param mh Numeric vector of corresponding parent-event magnitudes.
#' @param M0 Numeric magnitude threshold.
#' @param T1 Numeric lower bound of the integration interval.
#' @param T2 Numeric upper bound of the integration interval.
#' @param kernel Character string specifying the temporal triggering kernel.
#'   One of `"ou"`, `"mse"`, or `"rate_state"`.
#'
#' @return Numeric vector containing the log integrated triggered
#'   contributions.
#'
#' @details
#' For parent event \eqn{i}, the integrated contribution is
#' \eqn{\Lambda_i = K\exp\{\alpha(m_i-M_0)\} G_k(a_i,b_i)}, 
#' where \eqn{a_i=\max(T_1-t_i,0)} and \eqn{b_i=T_2-t_i}.
#'
#' The returned value is \eqn{\log \Lambda_i}.
log_Lambda_h <- function(theta, th, mh, M0, T1, T2, kernel = "ou") {
  kernel <- normalise_temporal_kernel(kernel)
  a <- pmax(T1 - th, 0)
  b <- T2 - th

  mass <- temporal_kernel_integral(a = a, b = b, theta = theta, kernel = kernel)
  
  log(theta$K) + theta$alpha * (mh - M0) + log(mass)
}
