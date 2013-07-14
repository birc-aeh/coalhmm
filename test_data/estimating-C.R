
posteriors <- cbind(read.table('marginal-posterior-0.txt',header=TRUE),
                    read.table('marginal-posterior-1.txt',header=TRUE)[,2],
                    read.table('marginal-posterior-2.txt',header=TRUE)[,2],
                    read.table('marginal-posterior-3.txt',header=TRUE)[,2],
                    read.table('marginal-posterior-4.txt',header=TRUE)[,2],
                    read.table('marginal-posterior-5.txt',header=TRUE)[,2],
                    read.table('marginal-posterior-6.txt',header=TRUE)[,2],
                    read.table('marginal-posterior-7.txt',header=TRUE)[,2],
                    read.table('marginal-posterior-8.txt',header=TRUE)[,2],
                    read.table('marginal-posterior-9.txt',header=TRUE)[,2])
names(posteriors) <- c('break.times',paste('posterior.',1:10,sep=''))

break.times <- posteriors$break.times
posteriors <- posteriors[,2:11]

nstates <- length(break.times)
delta.ts <- break.times[2:nstates] - break.times[1:(nstates-1)]

estimate.Cs <- function(post) {
  prob.left <- rev(cumsum(rev(post)))
  -log(1 - post[1:(nstates-1)] / prob.left[1:(nstates-1)]) / delta.ts
}

total.post <- apply(posteriors,1,mean)
total.Cs <- estimate.Cs(total.post)

bootstraps <- replicate(n=10, {
  boot.posteriors <- posteriors[,sample(10,replace=TRUE)]
  estimate.Cs(apply(boot.posteriors,1,mean))
})

plot(range(break.times), range(total.Cs), type='n', ylim=c(400,600),
     main='Estimation of coalescence rates over time',
     xlab='Time (mutations per basepair)', ylab='Coalescence rate (C)')
abline(v=break.times, lty='dashed',col='gray')
for (i in 1:ncol(bootstraps)) {
  segments(x0 = break.times[1:(nstates-1)], x1 = break.times[2:nstates], 
           y0 = as.vector(bootstraps[,i]), col='gray')  
}
segments(x0 = break.times[1:(nstates-1)], x1 = break.times[2:nstates], 
         y0 = total.Cs, col='blue', lwd=2)
abline(h=500, col='red')

plot(range(c(break.times,1.5*break.times[nstates])), range(total.post), type='n', ylim=c(0.05,0.15),
     main='Posterior probabilities',
     xlab='Time (mutations per basepair)', ylab='Posterior')
abline(v=break.times, lty='dashed',col='gray')
segments(x0 = break.times[1:(nstates)], x1 = c(break.times[2:nstates],1.5*break.times[nstates]), 
         y0 = total.post, col='blue', lwd=2)
abline(h=0.1, col='red')

