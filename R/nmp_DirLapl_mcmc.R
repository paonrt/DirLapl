#' @title Normal Means Problem with Dirichlet-Laplace Prior: Markov Chain Monte
#'  Carlo Algorithm

#' @description Gibbs sampler algorithm for the normal means problem with a
#'  Dirichlet-Laplace prior. 

#' @usage nmp_DirLapl_mcmc = function(nsim, y, alpha, burn = 0, thin = 1,
#'                                    verbose = +Inf)

#' @details The normal means problem is a statistical model in which a single
#' vector \eqn{y} is observed and is modeled as independent normal random
#' variables with known unit variance, say
#' \deqn{
#'  y \sim N_n (\theta, I_n)
#' }
#' where \eqn{I_n} is the \eqn{n}-dimensional identity matrix. The prior on the
#' mean vector \eqn{\theta} is the Dirichlet-Laplace prior:
#' \deqn{
#'  \theta_i \overset{ind}{\sim} N(0, \delta^2_i \psi_i) \, , \, \delta_i
#'   \overset{i.i.d.}{\sim} \mathrm{gamma}(\alpha, 1 / 2) \perp\!\!\!\perp
#'   \psi_i \overset{i.i.d.}{\sim} \exp(1 / 2)
#'  }
#' or, equivalently
#'  \deqn{
#'   \theta_i \overset{ind}{\sim} N(0, \tau^2 \phi^2_i \psi_i) \\
#'   \tau \sim \mathrm{gamma}(n \alpha, 1 / 2) \\
#'   \phi \sim \mathrm{Dirichlet}_n(\alpha, \alpha, \dots, \alpha) \\
#'   \psi_i \overset{i.i.d.}{\sim} \exp(1 / 2)
#'  }
#' For the Gibbs sampler, it is used the first representation.
# See Onorati et al (2025) for more details.
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

# @references
# Bhattacharya, A., Pati, D., Pillai, N. S., & Dunson, D. B. (2015).
#  Dirichlet–Laplace priors for optimal shrinkage. Journal of the American
#  Statistical Association, 110(512), 1479-1490.
#
# Onorati, P., Canale, A., & Dunson, D.B. (2025). On the Posterior
#  Computation Under the Dirichlet-Laplace Prior. arXiv preprint
#  arXiv:xxxx.xxxxx.

#' @param nsim The number of simulations.
#' @param y The response vector.
#' @param alpha The hyper-parameter controlling the amount of shrinkage.
#' @param burn The number of draws to be discarded as burn-in.
#' @param thin The thinning parameter.
#' @param verbose The period for printing status of the chain.
#' @return Posterior sample of the parameters.
#' @import GIGrvg
#' @import statmod
#' @export
nmp_DirLapl_mcmc = function(
  nsim, y, alpha, burn = 0, thin = 1, verbose = +Inf
) {

  ### print start time if required
  if (!is.infinite(verbose)) {
    print(paste(
      "nmp_DirLapl_mcmc: start time at ", Sys.time(), sep = ""
    ))
  }

  ### check verbose
  if (verbose == 0) verbose = nsim + 1

  ### check burn
  if (burn < 0) burn = 0




  ########## initialization

  ### sample size (and number of parameters)
  n = length(y)

  ### storing matrix
  keep_theta = matrix(0, nrow = n, ncol = nsim)

  ### initial value of theta (using JS estimator)
  theta = (1 - (n - 2) / sum(y^2)) * y
  theta[abs(theta) <= .Machine$double.eps] = .Machine$double.eps

  ### sample delta
  delta = apply(matrix(2 * abs(theta)), 1, FUN = function(xappl) {
    GIGrvg::rgig(n = 1, lambda = alpha - 1, chi = xappl, psi = 1)
  })

  ### sample psi
  psi = 1 / (statmod::rinvgauss(n, mean = delta / abs(theta), shape = 1))




  ########## draw the chain
  for (insim in (1 - burn):nsim) {

    ithin = 0

    ##### iterate up to thinning value
    while (ithin < thin) {

      ##########################################################################
      #                                                                        #
      #                              update theta                              #
      #                                                                        #
      ##########################################################################

      zetaSq = psi * delta^2 / (1 + psi * delta^2)
      mean_postTheta = zetaSq * y
      sd_postTheta = sqrt(zetaSq)
      tbs_rej = (.Machine$double.eps - abs(mean_postTheta)) / sd_postTheta < 1.3
      if (any(tbs_rej)) {
        theta[tbs_rej] = rng_rej(
          mean = mean_postTheta[tbs_rej], sd = sd_postTheta[tbs_rej]
        )
      }
      if (any(!tbs_rej)) {
        theta[!tbs_rej] = rng_trunc(
          mean = mean_postTheta[!tbs_rej], sd = sd_postTheta[!tbs_rej]
        )
      }


      ##########################################################################
      #                                                                        #
      #                              update delta                              #
      #                                                                        #
      ##########################################################################

      delta = apply(matrix(2 * abs(theta)), 1, FUN = function(xappl) {
        GIGrvg::rgig(n = 1, lambda = alpha - 1, chi = xappl, psi = 1)
      })



      ##########################################################################
      #                                                                        #
      #                               update psi                               #
      #                                                                        #
      ##########################################################################

      psi = 1 / (statmod::rinvgauss(n, mean = delta / abs(theta), shape = 1))



      ##### end of single complete scan
      if (insim > 0) ithin = ithin + 1 else ithin = thin

    }

    ### check if outside of burn-in
    if (insim > 0) {

      ### keep value
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
      "nmp_DirLapl_mcmc: end time at ", Sys.time(), sep = ""
    ))
  }

  ### return results
  return(keep_theta)

}
