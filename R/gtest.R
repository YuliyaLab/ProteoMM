#' G Test for presence - absence analysis
#'
#' Log-likelihood test for independence & goodness of fit.
#' g.test() performs Williams' and Yates' correction;
#' Monte Carlo simulation of p-values, via gtestsim.c.
#' MC requires recompilation of R.
#' Written by Peter Hurd (V3.3 Pete Hurd Sept 29 2001, phurd AT ualberta.ca).
#' Yuliya Karpievitch added comments for ease of understanding and
#' incorporated into ProteoMM.
#' G & q calculation from Sokal & Rohlf (1995) Biometry 3rd ed.,
#' TOI Yates correction taken from Mike Camanns 2x2 G-test function,
#' GOF Yates correction as described in Zar (2000),
#' more stuff taken from ctest's chisq.test().
#'
#'
#'
#' @param x vector of boolean values corresponding to presence & absence
#'          eg: c(TRUE, TRUE, FALSE, FALSE) for present present absent absent
#'          values. Order of TRUE/FALSE does not matter, can be used
#'          interchangeably. Same length as parameter y
#'
#' @param y vector treatments (factor) corresponding to values in x,
#'            same length as x
#'            eg: as.factor(c('grp1;, 'grp1', 'grp2', 'grp2'))
#'
#' @param correct correction to apply, options: "yates", "williams", "none"
#'            default: "none"
#'            NOTE: in ProteoMM we only tested & used correction = "none"
#'
#' @param p default: rep(1/length(x), length(x)), used in Yates correction
#'           NOTE: in ProteoMM we only tested & used the default parameter value
#'
#' @return htest object the following variables
#' \describe{
#'   \item{statistic}{value of the G statistic produced by g test}
#'   \item{parameter}{degrees of freedom of the test}
#'   \item{p.value}{p-value}
#'   \item{method}{method used to produce statistic and p-value}
#'   \item{data.name}{data passed in to the function}
#'   \item{observed}{matrix of observed counts}
#'   \item{expected}{matrix of expected counts}
#' }
#' @examples
#' g.test(c(TRUE, TRUE, FALSE, FALSE),
#'        as.factor(c('grp1', 'grp1', 'grp2', 'grp2')))
#' @export
g.test = function(x, y = NULL, correct="none",
  p = rep(1/length(x), length(x) ) )
{
  DNAME <- deparse(substitute(x))
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.matrix(x)) {
    if (min(dim(x)) == 1)
      x <- as.vector(x)
  }

  if (!is.matrix(x) && !is.null(y)) {
    if (length(x) != length(y))
      stop("x and y must have the same length")
    DNAME <- paste(DNAME, "and", deparse(substitute(y)))
    OK <- stats::complete.cases(x, y)
    x <- as.factor(x[OK])
    y <- as.factor(y[OK])
    if ((nlevels(x) < 2) || (nlevels(y) < 2))
      stop("x and y must have at least 2 levels")
    x <- table(x, y)
  }
  if (any(x < 0) || any(is.na(x)))
    stop("all entries of x must be nonnegative and finite")
  if ((n <- sum(x)) == 0)
    stop("at least one entry of x must be positive")
  # If x is matrix, do test of independence
  if (is.matrix(x)) {
    # Test of Independence
    nrows<-nrow(x)
    ncols<-ncol(x)
    if (correct=="yates"){ # Do Yates' correction?
      if(dim(x)[1]!=2 || dim(x)[2]!=2) # check for 2x2 matrix
        stop("Yates' correction requires a 2 x 2 matrix")
      if((x[1,1]*x[2,2])-(x[1,2]*x[2,1]) > 0)
        {
          x[1,1] <- x[1,1] - 0.5
          x[2,2] <- x[2,2] - 0.5
          x[1,2] <- x[1,2] + 0.5
          x[2,1] <- x[2,1] + 0.5
        }
      else
        {
          x[1,1] <- x[1,1] + 0.5
          x[2,2] <- x[2,2] + 0.5
          x[1,2] <- x[1,2] - 0.5
          x[2,1] <- x[2,1] - 0.5
        }
    }
    sr <- apply(x,1,sum)
    sc <- apply(x,2,sum)
    E <- outer(sr,sc, "*")/n
    # we are not doing a monte-carlo, calculate G

    # no monte-carlo
    # calculate G
    g <- 0
    for (i in seq_len(nrows)) {
      for (j in seq_len(ncols)) {
        if (x[i,j] != 0) g <- g + x[i,j] * log(x[i,j]/E[i,j])
      }
    }
    q <- 1
    if (correct=="williams"){ # Do Williams' correction
      row.tot = 0 
      col.tot = 0
      # vectorizing the for-loops
      # for (i in seq_len(nrows)) { row.tot <- row.tot + 1/(sum(x[i,])) }
      row.tot = sum( 1 / (rowSums(x) ) )
      # for (j in seq_len(ncols)) { col.tot <- col.tot + 1/(sum(x[,j])) }
      col.tot = sum( 1 / (colSums(x) ) )
      q <- 1+ ((n*row.tot-1)*(n*col.tot-1))/(6*n*(ncols-1)*(nrows-1))
    }
    STATISTIC <- G <- 2 * g / q
    PARAMETER <- (nrow(x)-1)*(ncol(x)-1)
    PVAL <- 1-stats::pchisq(STATISTIC,df=PARAMETER)
    if(correct=="none")
      METHOD =
      "Log likelihood ratio/G test of independence without correction"
    if(correct=="williams")
      METHOD =
      "Log likelihood ratio/G test of independence with Williams' correction"
    if(correct=="yates")
      METHOD =
      "Log likelihood ratio/G test of independence with Yates' correction"
  } else {
    # x is not a matrix, so we do Goodness of Fit
    METHOD = "Log likelihood ratio/G goodness of fit test"
    if (length(x) == 1)
      stop("x must at least have 2 elements")
    if (length(x) != length(p))
      stop("x and p must have the same number of elements")
    E <- n * p

    if (correct=="yates"){ # Do Yates' correction
      if(length(x)!=2)
        stop("Yates' correction requires 2 data values")
      if ( (x[1]-E[1]) > 0.25) {
        x[1] <- x[1]-0.5
        x[2] <- x[2]+0.5
      }
      else if ( (E[1]-x[1]) > 0.25){
        x[1] <- x[1]+0.5
        x[2] <- x[2]-0.5
      }
    }
    names(E) <- names(x)
    g = 0
    # for (i in seq_len(length(x)) ) {
    #  if (x[i] != 0) g <- g + x[i] * log(x[i]/E[i])
    # } # vectorized code for the above loop!
    ppos = x != 0
    g = sum(x[ppos] * log(x[ppos]/E[ppos]))
    
    q <- 1
    if (correct=="williams"){ # Do Williams' correction
      q <- 1+(length(x)+1)/(6*n)
    }
    STATISTIC <- G <- 2*g/q
    PARAMETER <- length(x) - 1
    PVAL <- stats::pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)
  }
  names(STATISTIC) = "Log likelihood ratio statistic (G)"
  names(PARAMETER) = "X-squared df"
  names(PVAL) = "p.value"
  structure(list(statistic=STATISTIC,parameter=PARAMETER,p.value=PVAL,
            method=METHOD,data.name=DNAME, observed=x, expected=E),
            class="htest")
}
