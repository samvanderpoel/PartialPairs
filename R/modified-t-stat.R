#' Kim et al.'s modified t-statistic
#'
#' \code{modified.tstat} uses Kim et al.'s modified t-statistic to obtain a
#' p-value for a partially matched pairs test.
#' 
#' Kim et al.’s modified t-statistic follows an approximately standard Gaussian
#' distribution under the null hypothesis. Mathematical details are provided in
#' [Kuan & Huang, 2013].
#' 
#' If proper sample size conditions are not met, then \code{modified.tstat} may
#' exit or perform a paired or unpaired two-sample t.test, depending on the
#' nature of the sample size issue.
#' 
#' If the variance of input data is close to zero, \code{modified.tstat} will
#' return an error message.
#'
#' @param x a non-empty numeric vector containing some NA values
#' @param y a non-empty numeric vector containing some NA values
#' @param alternative specification of the alternative hypothesis.
#' Takes values: \code{"two.sided"}, \code{"greater"}, or \code{"less"}.
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
#' modified.tstat(x, y, alternative = 'two.sided')
#' 
#' @references
#' Kuan, Pei Fen, and Bo Huang. "A simple and robust method for partially
#' matched samples using the p‐values pooling approach." Statistics in 
#' medicine 32.19 (2013): 3247-3259.
#'
#' @export
modified.tstat = function(x, y,
                           alternative = c('two.sided', 'greater', 'less')) {
     # check whether length(x)==length(y)
     if (length(x)!=length(y)) {
          if (sum(!is.na(x))<3 | sum(!is.na(y))<3) {
               stop('Sample sizes are too small and length of x ',
                    'should equal length of y.')
          } else {
               warning('Length of x should equal length of y. ',
                       'Two sample t-test attempted')
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
     if (n1<4 & n2+n3<5) {
             stop('Sample sizes are too small')
     } else if (n1>=4 & n2+n3<5) {
             warning('Not enough missing data for modified t-test. ',
                     'Matched pairs t-test attempted')
             return (t.test(pair.x, pair.y,
                            alternative = alternative, paired = TRUE)$p.value)
     } else if (n1<4 & n2+n3>=5) {
             warning('Not enough matched pairs for modified t-test. ',
                     'Two sample t-test attempted')
             return (t.test(x[only.x], y[only.y],
                            alternative = alternative)$p.value)
     }
     # else, n1>=4 and n2+n3>=5 is met, modified t-test is executed
     SD    = sd(pair.x-pair.y)
     t.bar = mean(x[only.x])
     n.bar = mean(y[only.y])
     ST = sd(x[only.x])
     SN = sd(y[only.y])
     # check whether variance of data is approx. zero
     if (ST < .Machine$double.eps * abs(t.bar) &
         SN < .Machine$double.eps * abs(n.bar) &
         SD < .Machine$double.eps *
              max(abs(mean(pair.x)), abs(mean(pair.y)))) {
          stop('Variance of data is too close to zero.')
     }
     nh = 2/(1/n2+1/n3)
     d.bar = mean(pair.x-pair.y)
     t3 = (n1*d.bar+nh*(t.bar-n.bar)) / sqrt(n1*SD^2 + nh^2*(ST^2/n2+SN^2/n3))
     alternative = match.arg(alternative)
     if (alternative == 'greater') {
          p.value = pnorm(t3, lower.tail = FALSE)
     } else if (alternative == 'less') {
          p.value = pnorm(t3, lower.tail = TRUE)
     } else if (alternative == 'two.sided') {
          p.value = 2*pnorm(abs(t3), lower.tail = FALSE)
     }
     return (p.value)
}
