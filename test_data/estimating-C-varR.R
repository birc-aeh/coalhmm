
analyse.R <- function(R) {
  posteriors <- cbind(read.table(paste('marginal-posterior-0-R-',R,'.txt',sep=''),header=TRUE),
                      read.table(paste('marginal-posterior-1-R-',R,'.txt',sep=''),header=TRUE)[,2],
                      read.table(paste('marginal-posterior-2-R-',R,'.txt',sep=''),header=TRUE)[,2],
                      read.table(paste('marginal-posterior-3-R-',R,'.txt',sep=''),header=TRUE)[,2],
                      read.table(paste('marginal-posterior-4-R-',R,'.txt',sep=''),header=TRUE)[,2],
                      read.table(paste('marginal-posterior-5-R-',R,'.txt',sep=''),header=TRUE)[,2],
                      read.table(paste('marginal-posterior-6-R-',R,'.txt',sep=''),header=TRUE)[,2],
                      read.table(paste('marginal-posterior-7-R-',R,'.txt',sep=''),header=TRUE)[,2],
                      read.table(paste('marginal-posterior-8-R-',R,'.txt',sep=''),header=TRUE)[,2],
                      read.table(paste('marginal-posterior-9-R-',R,'.txt',sep=''),header=TRUE)[,2])
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
  
  list(break.times = break.times, posteriors = total.post, C.estimates = total.Cs)
}

R.0.2 <- analyse.R("0.2")
R.0.3 <- analyse.R("0.3")
R.0.4 <- analyse.R("0.4")
nstates = length(R.0.2$break.times)
break.times = R.0.2$break.times

opar <- par(mfrow=c(2,1))

plot(range(c(break.times,1.5*break.times[nstates-1])),
     range(R.0.2$posteriors), type='n', ylim=c(0.05,0.15),
     xlab='Time (mutations per basepair)', ylab='Posterior')
abline(v=break.times, lty='dashed',col='gray')
segments(x0 = break.times[1:nstates], 
         x1 = c(break.times[2:nstates],1.5*break.times[nstates-1]), 
         y0 = R.0.2$posteriors, col='black')  
segments(x0 = break.times[1:nstates], 
         x1 = c(break.times[2:nstates],1.5*break.times[nstates-1]), 
         y0 = R.0.3$posteriors, col='blue')
segments(x0 = break.times[1:nstates], 
         x1 = c(break.times[2:nstates],1.5*break.times[nstates-1]), 
         y0 = R.0.4$posteriors, col='darkgreen')
abline(h=0.1, col='red')

legend('topright',
       legend=c('R = 0.2','R = 0.3','R = 0.4'),
       fill=c('black','blue','darkgreen'))

plot(range(break.times), range(R.0.2$C.estimates), type='n', ylim=c(400,600),
     xlab='Time (mutations per basepair)', ylab='Coalescence rate (C)')
abline(v=break.times, lty='dashed',col='gray')
segments(x0 = break.times[1:(nstates-1)], x1 = break.times[2:nstates], 
         y0 = R.0.2$C.estimates, col='black')  
segments(x0 = break.times[1:(nstates-1)], x1 = break.times[2:nstates], 
         y0 = R.0.3$C.estimates, col='blue')  
segments(x0 = break.times[1:(nstates-1)], x1 = break.times[2:nstates], 
         y0 = R.0.4$C.estimates, col='darkgreen')  
abline(h=500, col='red')

par(opar)