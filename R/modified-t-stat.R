#' Kim et al.'s modified t-statistic
#'
#' Uses Kim et al.'s modified t-statistic to obtain a p-value for
#' partially matched pairs.
#'
#' @param x a non-empty numeric vector of data values
#' @param y a non-empty numeric vector of data values
#' @param alternative specification of the alternative hypothesis.
#' Takes values: "two.sided", "greater", or "less".
#'
#' @return p-value corresponding with the hypothesis test
#'
#' @examples
#' 
#'
#' @export
modified.t.stat = function(x, y,
                           alternative = c('two-sided', 'greater', 'less')) {
     full.sample.inds = !is.na(x) & !is.na(y)
     only.x = !is.na(x) & is.na(y)
     only.y = !is.na(y) & is.na(x)
     n1 = sum(full.sample.inds)
     n2 = sum(only.x)
     n3 = sum(only.y)
     if (n1<4 & n2+n3<5) {
             stop('Sample sizes are too small')
     } else if (n1>=4 & n2+n3<5) {
             warning('Not enough missing data for modified t-test. ',
                     'Matched pairs t-test executed.')
             return (t.test(x[full.sample.inds], y[full.sample.inds],
                            alternative = alternative, paired = TRUE)$p.value)
     } else if (n1<4 & n2+n3>=5) {
             warning('Not enough matched pairs for modified t-test. ',
                     'Two sample t-test executed.')
             return (t.test(x[only.x], y[only.y],
                            alternative = alternative)$p.value)
     }
     # if n1>=4 and n2+n3>=5, modified t-test is executed
     nh = 2/(1/n1+1/n2)
     d.bar = mean(x[full.sample.inds]-y[full.sample.inds])
     SD    = sd(x[full.sample.inds]-y[full.sample.inds])
     t.bar = mean(x[only.x]); ST = sd(x[only.x])
     n.bar = mean(y[only.y]); SN = sd(y[only.y])
     t3 = (n1*d.bar+nh*(t.bar-n.bar)) / sqrt(n1*SD^2 + nh^2*(ST^2/n2+SN^2/n3))
     if (alternative == 'greater') {
          p.value = pnorm(t3, lower.tail = FALSE)
     } else if (alternative == 'less') {
          p.value = pnorm(t3, lower.tail = TRUE)
     } else if (alternative == 'two-sided') {
          p.value = 2*pnorm(abs(t3), lower.tail = FALSE)
     }
     return (p.value)
}
