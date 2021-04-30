#' Liptak weighted Z-test
#'
#' \code{liptak.ztest} returns the p-value associated with Liptak's weighted
#' z-test for partially matched pairs
#' 
#' Liptak's weighted Z-test computes the p-values of a paired sample t-test on
#' the n1 paired entries and of a two-sample t-test on the n2+n3 unpaired
#' entries (see vignette for further details on the paired vs. unpaired
#' distinction). The two p-values are then weighted and combined as detailed in
#' [Kuan & Huang, 2013].
#' 
#' If proper sample size conditions are not met, then \code{liptak.ztest} may
#' exit or perform a paired or unpaired two-sample t.test, depending on the
#' nature of the sample size issue.
#' 
#' If the variance of input data is close to zero, \code{liptak.ztest} will
#' return an error message.
#'
#' @param x a non-empty numeric vector of data values
#' @param y a non-empty numeric vector of data values
#' @param alternative specification of the alternative hypothesis.
#' Takes values: \code{two.sided}, \code{greater}, or \code{less}.
#'
#' @return p-value associated with the hypothesis test
#'
#' @examples
#' In the following, the true means are not equal:
#' 
#' x = rnorm(400, 0, 1)
#' x[sample(1:400, size=75, replace=FALSE)] = NA
#' y = rnorm(400, 0.4, 3)
#' y[sample(1:400, size=75, replace=FALSE)] = NA
#' liptak.ztest(x, y, alternative = 'two.sided')
#' 
#' @references
#' Kuan, Pei Fen, and Bo Huang. "A simple and robust method for partially
#' matched samples using the p‐values pooling approach." Statistics in 
#' medicine 32.19 (2013): 3247-3259.
#'
#' @export
liptak.ztest = function(x, y, alternative=c('two.sided', 'greater', 'less')) {
     # check whether length(x)==length(y)
     if (length(x)!=length(y)) {
          if (sum(!is.na(x))<3 | sum(!is.na(y))<3) {
               stop('Sample sizes are too small and length of x ',
                    'should equal length of y.')
          } else {
               warning('Length of x should equal length of y. ',
                       'Two sample t-test performed.')
               return (t.test(x[!is.na(x)], y[!is.na(y)])$p.value)
          }
     }
     pair.inds = !is.na(x) & !is.na(y)
     only.x = !is.na(x) & is.na(y)
     only.y = !is.na(y) & is.na(x)
     pair.x = x[pair.inds]
     pair.y = y[pair.inds]
     # test whether appropriate sample size conditons are met
     n1 = sum(pair.inds)
     n2 = sum(only.x)
     n3 = sum(only.y)
     if (n1<3 & n2+n3<5) {
          stop('Sample sizes are too small or too much missing data.')
     } else if (n1>=3 & n2+n3<5) {
          warning('Not enough missing data for Liptak z-test. ',
                  'Matched pairs t-test performed.')
          return (t.test(pair.x, pair.y,
                         alternative = alternative, paired = TRUE)$p.value)
     } else if (n1<3 & n2+n3>=5) {
          warning('Not enough matched pairs for Liptak z-test. ',
                  'Two sample t-test performed.')
          return (t.test(x[only.x], y[only.y],
                         alternative = alternative)$p.value)
     }
     # else, n1>=3 and n2+n3>=5 is met, Liptak's z-test is performed.
     # check whether variance of data is approx. zero
     if ((sd(x[only.x]) < 10 * .Machine$double.eps * abs(mean(x[only.x]))  &
          sd(y[only.y]) < 10 * .Machine$double.eps * abs(mean(y[only.y]))) |
          sd(pair.x-pair.y) < 10 *.Machine$double.eps *
                              max(abs(mean(pair.x)), abs(mean(pair.y))))   {
          stop('Variance of data is too close to zero.')
     }
     alternative = match.arg(alternative)
     if (alternative == 'greater') {
          p1 = t.test(pair.x, pair.y,
                      paired = TRUE, alternative = 'greater')$p.value
          p2 = t.test(x[only.x], y[only.y], alternative = 'greater',
                      var.equal = FALSE)$p.value
     } else {
          p1 = t.test(pair.x, pair.y,
                      paired = TRUE, alternative = 'less')$p.value
          p2 = t.test(x[only.x], y[only.y], alternative = 'less',
                      var.equal = FALSE)$p.value
     }
     w1 = sqrt(2*n1)
     w2 = sqrt(n2+n3)
     Z1 = qnorm(1-p1)
     Z2 = qnorm(1-p2)
     p.comb = 1 - pnorm( (w1*Z1+w2*Z2)/sqrt(w1^2+w2^2) )
     if (alternative == 'two.sided') {
          if (p.comb<0.5) {p.comb=2*p.comb} else {p.comb=2*(1-p.comb)}
     }
     return (p.comb)
}
