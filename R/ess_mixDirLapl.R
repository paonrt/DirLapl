#' @title Elliptical Slice Sampling for Mixture of Dirichlet-Laplace Priors

#' @description Elliptical slice sampling for a generic likelihood under a
#'  mixture of Dirichlet-Laplace priors for the parameters. 

#' @usage ess_mixDirLapl = function(nsim, logL, start, burn = 0, thin = 1,
#'                               verbose = +Inf)

#' @details Generic likelihood, say
#' \deqn{
#'  y_i \overset{i.i.d.}{\sim} L(\cdot; \theta)
#' }
#' The prior on the vector parameter \eqn{\theta} is the Dirichlet-Laplace
#' prior:
#' \deqn{
#'  \theta_i \overset{ind}{\sim} N(0, \delta^2_i \psi_i) \, , \,
#'   \delta_i \overset{i.i.d.}{\sim} \mathrm{gamma}(\alpha, 1 / 2)
#'   \perp\!\!\!\perp \psi_i \overset{i.i.d.}{\sim} \exp(1 / 2)
#'  }
#' or, equivalently
#'  \deqn{
#'   \theta_i \overset{ind}{\sim} N(0, \tau^2 \phi^2_i \psi_i) \\
#'   \tau \sim \mathrm{gamma}(n \alpha, 1 / 2) \\
#'   \phi \sim \mathrm{Dirichlet}_n(\alpha, \alpha, \dots, \alpha) \\
#'   \psi_i \overset{i.i.d.}{\sim} \exp(1 / 2)
#'  }
#' and the prior of \eqn{\alpha} is uniform on the discrete support
#' \deqn{
#'  \frac{1}{p} + (\frac{1}{2} - \frac{1}{p}) \frac{i}{199} \, , \, i = 0, 1,
#'  \dots, 199
#' }
#' say, the support of \eqn{\alpha} is
#' \code{seq(from = 1 / p, to = 0.5, length.out = 200)}.
#' For the update of \eqn{p(\theta \vert y, \delta, \psi) \propto L(y; \theta)
#' N(\theta, 0, \delta^2 \psi)}, it is used the elliptical slice sampling
#' (Murray et al, 2010).
# See Onorati et al (2025) for the update of \eqn{delta, \psi}.
#' 
#' \code{start}: the starting point must be provided. It must be a numeric
#'  vector
#' 
#' \code{logL}: sum of the loglikelihood of the single observations.
#'
#' Only one value every \code{thin} values is kept in the chain, so the true
#'  number of complete scans will be \code{nsim * thin + burn}. By default
#'  \code{thin = 1}, that is no thinning.
#'
#' The current time and the current number of iteration are printed one every
#' \code{verbose} iterations. Furthermore:
#' \itemize{
#'  \item if \code{verbose == +-Inf} then there is no printing,
#'  \item if \code{verbose != +-Inf} then at least start and end of simulation
#'   are reported.
#' }

# @references Bhattacharya, A., Pati, D., Pillai, N. S., & Dunson, D. B. (2015).
# Dirichlet–Laplace priors for optimal shrinkage. Journal of the American
#  Statistical Association, 110(512), 1479-1490.
#
# Murray, I., Adams, R., & MacKay, D. (2010, March). Elliptical slice sampling.
# In Proceedings of the thirteenth international conference on artificial
# intelligence and statistics (pp. 541-548). JMLR Workshop and Conference
# Proceedings.
#
# Onorati, P., Canale, A., & Dunson, D.B. (2025). On the Posterior Computation
#  Under the Dirichlet-Laplace Prior. arXiv preprint arXiv:xxxx.xxxxx.

#' @param nsim The number of simulations.
#' @param logL The loglikelihood function.
#' @param start The starting point.
#' @param burn The number of draws to be discarded as burn-in.
#' @param thin The thinning parameter.
#' @param verbose The period for printing status of the chain.
#' @return Posterior sample of the parameters.
#' @import GIGrvg
#' @import statmod
#' @export
ess_mixDirLapl = function(
  nsim, logL, start = NULL, burn = 0, thin = 1, verbose = +Inf
) {

  ### print start time if required
  if (!is.infinite(verbose)) {
    print(paste(
      "ess_mixDirLapl: start time at ", Sys.time(), sep = ""
    ))
  }
  
  ### check start
  if (!is.numeric(start)) stop("start must be numeric")
  
  ### check logL
  if (length(logL(start)) != 1) stop("logL does not provide a single value")

  ### check verbose
  if (verbose == 0) verbose = nsim + 1

  ### check burn
  if (burn < 0) burn = 0




  ########## initialization

  ### number of parameters
  p = length(start)
  
  ### get support of alpha
  seq_alpha = seq(from = 1 / p, to = 0.5, length.out = 200)

  ### compute log-normalizing constant of full conditional of alpha
  log_const_alpha = - (0.5 * (1 + seq_alpha) * log(2) + lgamma(seq_alpha))

  ### get initial value
  theta = start
  
  ### compute loglikelihood at the initial value
  logL_curr = logL(theta)
  
  ### storing objects
  keep_theta = matrix(0, nrow = p, ncol = nsim)



  ########## draw the chain
  for (insim in (1 - burn):nsim) {

    ithin = 0



    ##### iterate up to thinning value
    while (ithin < thin) {
      
      ##########################################################################
      #                                                                        #
      #                      update alpha, delta, and psi                      #
      #                                                                        #
      ##########################################################################

      ### sample alpha
      log_wei_alpha = log_const_alpha + log_kernel_alpha(
        abs_theta = drop(abs(theta)),
        seq_alpha = seq_alpha,
        p = p
      )
      log_wei_alpha = log_wei_alpha - max(log_wei_alpha)
      alpha = sample(x = seq_alpha, size = 1, prob = exp(log_wei_alpha))

      ### sample delta
      delta = apply(
        matrix(2 * abs(theta)), 1, FUN = function(xappl) {
          tryCatch(
            expr = GIGrvg::rgig(
              n = 1, lambda = alpha - 1, chi = xappl, psi = 1
            ),
            error = function(err) {
              GIGrvg::rgig(
                n = 1, lambda = alpha - 1,
                chi = xappl + sqrt(.Machine$double.eps),
                psi = 1
              )
            }
          )
        }
      )

      ### sample psi
      psi = 1 / (statmod::rinvgauss(
        p, mean = delta / abs(theta), shape = 1
      )) #+ sqrt(.Machine$double.eps))
      
      

      ##########################################################################
      #                                                                        #
      #                              update theta                              #
      #                                                                        #
      ##########################################################################

      ### sample from the prior
      theta_prior = rnorm(n = p, mean = 0, sd = delta * sqrt(psi))
      
      ### loglikelihood threshold
      strike = logL_curr - rexp(n = 1, rate = 1)
      
      ### sample initial angle
      angle = runif(n = 1, min = 0, max = 2 * pi)
      
      ### get initial brackets
      angle_min = angle - 2 * pi
      angle_max = angle
      
      ### run up to success
      while (TRUE) {
        
        ### get proposal
        theta_prop = theta * cos(angle) + theta_prior * sin(angle)
        
        ### compute loglikelihood of the proposal
        logL_prop = logL(theta_prop)
        
        ### check success
        if (logL_prop > strike) {
          
          ### if success then update theta and current loglikelihood
          theta = theta_prop
          logL_curr = logL_prop
          
          ### stop while loop at the first success
          break
          
          ### unsuccess case
        } else {
          
          ### shrink brackets
          if (angle < 0) angle_min = angle else angle_max = angle
          
          ### get new angle
          angle = runif(1, min = angle_min, max = angle_max)
          
        }
        
      }



      ##### end of single complete scan
      if (insim > 0) ithin = ithin + 1 else ithin = thin

    }



    ### check if outside of burn-in
    if (insim > 0) {

      ### keep values
      keep_theta[, insim] = theta

      ### print status of the chain
      if (insim %% verbose == 0) {
        print(paste(
          "iteration ", insim, " of ", nsim, " completed at time ", Sys.time(),
          sep = ""
        ))
      }

    } else {

      ### print status of the chain during burn-in
      if ((insim + burn) %% verbose == 0) {
        print(paste(
          "iteration ", insim, " of ", nsim, " completed at time ", Sys.time(),
          sep = ""
        ))
      }

    }

  }




  ### print end time if required
  if (!is.infinite(verbose)) {
    print(paste(
      "ess_mixDirLapl: end time at ", Sys.time(), sep = ""
    ))
  }

  ### return results
  return(keep_theta)

}
