#' @title Logit Model with Dirichlet-Laplace Prior: Markov Chain Monte Carlo
#'  Algorithm

#' @description Gibbs sampler algorithm for the logit model with no intercept.
#'  A Dirichlet-Laplace prior is used for the vector of the coefficients. 

#' @usage logit_DirLapl_mcmc = function(nsim, X, y, alpha, start = NULL,
#'                                      burn = 0, thin = 1, verbose = +Inf)

#' @details Logistic regression model, say
#' \deqn{
#'  y_i \overset{ind}{\sim} Bernoulli (X^\prime_i \beta)
#' }
#' where \eqn{X^\prime_i} is the \eqn{i}-th row of the design matrix \eqn{X}.
#' The prior on the coefficient vector \eqn{\beta} is the Dirichlet-Laplace
#' prior:
#' \deqn{
#'  \beta_i \overset{ind}{\sim} N(0, \delta^2_i \psi_i) \, , \,
#'   \delta_i \overset{i.i.d.}{\sim} \mathrm{gamma}(\alpha, 1 / 2)
#'   \perp\!\!\!\perp \psi_i \overset{i.i.d.}{\sim} \exp(1 / 2)
#'  }
#' or, equivalently
#'  \deqn{
#'   \beta_i \overset{ind}{\sim} N(0, \tau^2 \phi^2_i \psi_i) \\
#'   \tau \sim \mathrm{gamma}(n \alpha, 1 / 2) \\
#'   \phi \sim \mathrm{Dirichlet}_n(\alpha, \alpha, \dots, \alpha) \\
#'   \psi_i \overset{i.i.d.}{\sim} \exp(1 / 2)
#'  }
#' For the Gibbs sampler, it is used the first representation.
# See Onorati et al (2025).
#' 
#' \code{start}: the starting point is generated in the following way:
#' \itemize{
#'  \item if \code{start == NULL} then set the starting point of \code{beta}
#'        from elastic net.
#'  \item \code{else} set the starting point for \code{beta} to \code{start}.
#' }
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
# Onorati, P., Canale, A., & Dunson, D.B. (2025). On the Posterior Computation
#  Under the Dirichlet-Laplace Prior. arXiv preprint arXiv:xxxx.xxxxx.

#' @param nsim The number of simulations.
#' @param X The design matrix.
#' @param y The response vector.
#' @param alpha The hyper-parameter controlling the amount of shrinkage.
#' @param start The starting point.
#' @param burn The number of draws to be discarded as burn-in.
#' @param thin The thinning parameter.
#' @param verbose The period for printing status of the chain.
#' @return Posterior sample of the parameters.
#' @import glmnet
#' @import GIGrvg
#' @import statmod
#' @import BayesLogit
#' @export
logit_DirLapl_mcmc = function(
  nsim, X, y, alpha, start = NULL, burn = 0, thin = 1, verbose = +Inf
) {

  ### print start time if required
  if (!is.infinite(verbose)) {
    print(paste(
      "logit_DirLapl_mcmc: start time at ", Sys.time(), sep = ""
    ))
  }

  ### check verbose
  if (verbose == 0) verbose = nsim + 1

  ### check burn
  if (burn < 0) burn = 0




  ########## initialization

  ### sample size
  n = nrow(X)

  ### number of parameters
  p = ncol(X)
  
  ### compute t(X) %*% (y - 1/2)
  tXyOh = t(X) %*% (y - 0.5)

  ### storing objects
  keep_beta = matrix(0, nrow = p, ncol = nsim)



  ##### get initial values

  ### if initial values are not provided then use elastic net
  if (is.null(start)) {
    
    beta = as.numeric(coef(
      glmnet::cv.glmnet(
        x = X, y = factor(y), alpha = 0.25,
        intercept = FALSE, family = "binomial", 
      ),
      s = "lambda.min"
    ))[-1]

    ### otherwise get initial values
  } else {
    
    beta = start

  }




  ########## draw the chain
  for (insim in (1 - burn):nsim) {

    ithin = 0



    ##### iterate up to thinning value
    while (ithin < thin) {
      
      ##########################################################################
      #                                                                        #
      #                          update delta and psi                          #
      #                                                                        #
      ##########################################################################

      ### sample delta
      delta = apply(
        matrix(2 * abs(beta)), 1, FUN = function(xappl) {
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
        p, mean = delta / abs(beta), shape = 1
      )) #+ sqrt(.Machine$double.eps))
      
      

      ##########################################################################
      #                                                                        #
      #                         update beta and omega                          #
      #                                                                        #
      ##########################################################################

      ##### update omega
      omega = BayesLogit::rpg.devroye(n, 1, drop(X %*% beta))
      
      
      
      ##### update beta

      ### decomposition of variance-covariance matrix of beta given omega and y
      #    (V_omega)
      # upper Cholesky factor of the inverse
      decompV_omega = chol(
        t(X) %*% (X * omega) + diag(1 / (psi * delta^2), nrow = p)
      )
      # inverse of the upper Cholesky factor of the inverse
      decompV_omega = backsolve(
        r = decompV_omega,
        x = diag(1, nrow = p)
      )

      ### get V_omega
      V_omega = tcrossprod(decompV_omega)

      ### sample beta
      beta = V_omega %*% tXyOh + decompV_omega %*% rnorm(p)



      ##### end of single complete scan
      if (insim > 0) ithin = ithin + 1 else ithin = thin

    }



    ### check if outside of burn-in
    if (insim > 0) {

      ### keep values
      keep_beta[, insim] = beta

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
      "logit_DirLapl_mcmc: end time at ", Sys.time(), sep = ""
    ))
  }

  ### return results
  return(keep_beta)

}
