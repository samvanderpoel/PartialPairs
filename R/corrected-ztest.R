#' Looney and Jones corrected Z-test
#' 
#' \code{corrected.ztest} uses the Looney and Jones corrected z-test
#' to obtain a p-value for a partially matched pairs test.
#' 
#' Looney and Jones’s corrected Z-test is ``corrected'' in the sense that it
#' adjusts the standard error of the difference of the two samples by accounting
#' for the correlation between the $n_1$ paired observations. Under the null
#' hypothesis, the resulting test statistic Z_corr has an asymptotic N(0,1)
#' distribution. 
#' 
#' If proper sample size conditions are not met, then \code{corrected.ztest} may
#' exit or perform a paired or unpaired two-sample t.test, depending on the
#' nature of the sample size issue.
#' 
#' If the variance of input data is close to zero, \code{corrected.ztest} will
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
#' corrected.ztest(x, y, alternative = 'two.sided')
#' 
#' @references
#' Kuan, Pei Fen, and Bo Huang. "A simple and robust method for partially
#' matched samples using the p‐values pooling approach." Statistics in 
#' medicine 32.19 (2013): 3247-3259.
#'
#' @export
corrected.ztest = function(x, y,
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
     # if n1>=4 and n2+n3>=5, modified t-test is executed
     T.bar = mean(x[!is.na(x)])
     N.bar = mean(y[!is.na(y)])
     ST = sd(x[!is.na(x)])
     SN = sd(y[!is.na(y)])
     # check whether variance of data is approx. zero
     if (ST < .Machine$double.eps * abs(T.bar) &
         SN < .Machine$double.eps * abs(N.bar)){
          stop('Variance of data is almost zero')
     }
     STN1 = cov(pair.x, pair.y)
     se = sqrt( ST^2/(n1+n2) + SN^2/(n1+n3) - 2*n1*STN1/((n1+n2)*(n1+n3)) )
     z.corr = (T.bar - N.bar) / se
     if (all(alternative == 'greater')) {
          p.value = pnorm(z.corr, lower.tail = FALSE)
     } else if (all(alternative == 'less')) {
          p.value = pnorm(z.corr, lower.tail = TRUE)
     } else if (all(alternative == 'two.sided')) {
          p.value = 2*pnorm(abs(z.corr), lower.tail = FALSE)
     }
     return (p.value)
}
