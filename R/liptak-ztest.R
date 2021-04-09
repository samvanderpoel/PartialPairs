# Liptak weighted Z-test

liptak.ztest = function(x, y, alternative=c('two-sided', 'greater', 'less')) {
     if (length(x)!=length(y) & (sum(!is.na(x))<3 | sum(!is.na(y))<3)) {
          stop('Sample sizes are too small, and length of x should ',
               'equal length of y.')
     } else if (length(x)!=length(y) & sum(!is.na(x))>=3 & sum(!is.na(y))>=3) {    
          warning('Length of x should equal length of y. ',
                  'Two sample t-test performed.')
          return (t.test(x[!is.na(x)], y[!is.na(y)])$p.value)
     }
     
     full.sample.inds = !is.na(x) & !is.na(y); n1 = sum(full.sample.inds)
     only.x = !is.na(x) & is.na(y);            n2 = sum(only.x)
     only.y = !is.na(y) & is.na(x);            n3 = sum(only.y)
     if (n1<3 & n2+n3<5) {
          stop('Sample sizes are too small')
     } else if (n1>=3 & n2+n3<5) {
          warning('Not enough missing data for Liptak z-test. ',
                  'Matched pairs t-test performed.')
          return (t.test(x[full.sample.inds], y[full.sample.inds],
                         alternative = alternative, paired = TRUE)$p.value)
     } else if (n1<3 & n2+n3>=5) {
          warning('Not enough matched pairs for Liptak z-test. ',
                  'Two sample t-test performed.')
          return (t.test(x[only.x], y[only.y],
                         alternative = alternative)$p.value)
     }
     # if n1>=3 and n2+n3>=5, Liptak's z-test is performed.
     
     len = length(strsplit(alternative, split='')[[1]])
     alt.str.comp = strsplit('alternative', split='')[[1]][1:len]
     if (all(alternative == alt.str.comp)) {
          p1 = t.test(x[full.sample.inds], y[full.sample.inds],
                      paired = TRUE, alternative = 'greater')$p.value
          p2 = t.test(x[only.x], y[only.y], alternative = 'greater',
                      var.equal = FALSE)$p.value
     } else {
          p1 = t.test(x[full.sample.inds], y[full.sample.inds],
                      paired = TRUE, alternative = 'less')$p.value
          p2 = t.test(x[only.x], y[only.y], alternative = 'less',
                      var.equal = FALSE)$p.value
     }
     w1 = sqrt(2*n1); w2 = sqrt(n2+n3)
     Z1 = qnorm(1-p1); Z2 = qnorm(1-p2)
     p.comb = 1 - pnorm( (w1*Z1+w2*Z2)/sqrt(w1^2+w2^2) )
     if (all(alternative == 'two-sided')) {
          if (p.comb<0.5) {p.comb=2*p.comb} else {p.comb=2*(1-p.comb)}
     }
     return (p.comb)
}
