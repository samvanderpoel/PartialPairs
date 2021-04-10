#' Lin and Stivers’s MLE-based test under heteroscedasticity
#' 
#' \code{lin.mle.test} uses the Lin and Stivers’s MLE-based test 
#' under heteroscedasticity to obtain a p-value for a partially matched
#' pairs test.
#' 
#' These are the details
#'
#' @param x a non-empty numeric vector of data values
#' @param y a non-empty numeric vector of data values
#' @param alternative specification of the alternative hypothesis.
#' Takes values: "two.sided", "greater", or "less".
#'
#' @return p-value corresponding with the hypothesis test
#'
#' @examples
#' This is an example.
#' 
#' @references
#' Kuan, Pei Fen, and Bo Huang. "A simple and robust method for partially
#' matched samples using the p‐values pooling approach." Statistics in 
#' medicine 32.19 (2013): 3247-3259.
#'
#' @export
lin.mle.test = function(x, y,
                        alternative = c('two.sided', 'greater', 'less')) {
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
     if (n1<4 & n2+n3<5) {
             stop('Sample sizes are too small')
     } else if (n1>=4 & n2+n3<5) {
             warning('Not enough missing data for modified t-test. ',
                     'Matched pairs t-test executed.')
             return (t.test(pair.x, pair.y,
                            alternative = alternative, paired = TRUE)$p.value)
     } else if (n1<4 & n2+n3>=5) {
             warning('Not enough matched pairs for modified t-test. ',
                     'Two sample t-test executed.')
             return (t.test(x[only.x], y[only.y],
                            alternative = alternative)$p.value)
     }
     # else, n1>=4 and n2+n3>=5 is met, Lin MLE test is executed
     T1.bar = mean(pair.x)
     N1.bar = mean(pair.y)
     ST1 = sd(pair.x)
     SN1 = sd(pair.y)
     # check whether variance of data is approx. zero
     if (ST1 < 10*.Machine$double.eps * abs(T1.bar) |
         SN1 < 10*.Machine$double.eps * abs(N1.bar)){
          stop('Variance of data is too close to zero')
     }
     T.bar = mean(x[!is.na(x)])
     N.bar = mean(y[!is.na(y)])
     STN1 = cov(pair.x, pair.y)
     r = STN1 / (ST1 * SN1)
     f = n1*(n1+n3+n2*STN1/ST1^2) / ((n1+n2)*(n1+n3)-n2*n3*r^2)
     g = n1*(n1+n2+n3*STN1/SN1^2) / ((n1+n2)*(n1+n3)-n2*n3*r^2)
     V1 = ((f^2/n1 + (1-f)^2/n2)*(n1-1)*ST1^2 + 
           (g^2/n1 + (1-g)^2/n3)*(n1-1)*SN1^2 -
            2*f*g*STN1*(n1-1)/n1) / (n1-1)
     Z.ls = (f*(T1.bar-T.bar) - g*(N1.bar-N.bar) + T.bar - N.bar) / sqrt(V1)
     if (all(alternative == 'greater')) {
          p.value = pt(Z.ls, n1, lower.tail = FALSE)
     } else if (all(alternative == 'less')) {
          p.value = pt(Z.ls, n1, lower.tail = TRUE)
     } else if (all(alternative == 'two.sided')) {
          p.value = 2*pt(abs(Z.ls), n1, lower.tail = FALSE)
     }
     return (p.value)
}   
