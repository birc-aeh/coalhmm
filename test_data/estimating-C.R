
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

posteriors <- cbind(read.table('marginal-posterior-20states-0.txt',header=TRUE),
                    read.table('marginal-posterior-20states-1.txt',header=TRUE)[,2],
                    read.table('marginal-posterior-20states-2.txt',header=TRUE)[,2],
                    read.table('marginal-posterior-20states-3.txt',header=TRUE)[,2],
                    read.table('marginal-posterior-20states-4.txt',header=TRUE)[,2],
                    read.table('marginal-posterior-20states-5.txt',header=TRUE)[,2],
                    read.table('marginal-posterior-20states-6.txt',header=TRUE)[,2],
                    read.table('marginal-posterior-20states-7.txt',header=TRUE)[,2],
                    read.table('marginal-posterior-20states-8.txt',header=TRUE)[,2],
                    read.table('marginal-posterior-20states-9.txt',header=TRUE)[,2])
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

opar <- par(mfrow=c(2,1))
par(mar=c(1, 4, 2, 2) + 0.1)
plot(range(break.times), range(total.Cs), type='n', ylim=c(400,600),
     #main='Estimation for individual chunks',
     xlab='Time (mutations per basepair)', ylab='Coalescence rate (C)', xaxt='n')
abline(v=break.times, lty='dashed',col='gray')
for (i in 1:ncol(posteriors)) {
  segments(x0 = break.times[1:(nstates-1)], x1 = break.times[2:nstates], 
           y0 = estimate.Cs(posteriors[,i]), col='darkgreen')  
}
abline(h=500, col='red')

par(mar=c(2, 4, 1, 2) + 0.1)
plot(range(break.times), range(total.Cs), type='n', ylim=c(400,600),
     #main='Estimation for combined chunks',
     xlab='Time (mutations per basepair)', ylab='Coalescence rate (C)')
abline(v=break.times, lty='dashed',col='gray')
for (i in 1:ncol(bootstraps)) {
  segments(x0 = break.times[1:(nstates-1)], x1 = break.times[2:nstates], 
           y0 = bootstraps[,i], col='gray')  
}
segments(x0 = break.times[1:(nstates-1)], x1 = break.times[2:nstates], 
         y0 = total.Cs, col='blue', lwd=2)
abline(h=500, col='red')
par(opar)

opar <- par(mfrow=c(2,1))
par(mar=c(1, 4, 2, 2) + 0.1)
plot(range(c(break.times,1.5*break.times[nstates])), range(posteriors), type='n',
     #     main='Posterior probabilities',
     xaxt='n',
     xlab='Time (mutations per basepair)', ylab='Posterior')
abline(v=break.times, lty='dashed',col='gray')
for (i in 1:ncol(posteriors)) {
  segments(x0 = break.times[1:(nstates)], x1 = c(break.times[2:nstates],1.5*break.times[nstates]), 
           y0 = posteriors[,i], col='darkgreen')  
}

par(mar=c(2, 4, 1, 2) + 0.1)
plot(range(c(break.times,1.5*break.times[nstates])), range(total.post), type='n',
#     main='Posterior probabilities',
     xlab='Time (mutations per basepair)', ylab='Posterior')
abline(v=break.times, lty='dashed',col='gray')
segments(x0 = break.times[1:(nstates)], x1 = c(break.times[2:nstates],1.5*break.times[nstates]), 
         y0 = total.post, col='blue', lwd=2)
par(opar)
