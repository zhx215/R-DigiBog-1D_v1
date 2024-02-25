## moving average ##
x_scale <- c(0,t_extent)
window <- 50
yr <- seq(1+window/2,t_extent-window/2,1)
rm <- layer_remaining_mass[1:t_extent]
mb <- (prod_output-decay_output)*10000 # g m-2 yr-1, mass balance
mar <- layer_mass[1:t_extent]*10000 # g m-2 yr-1, mass accumulation rate
wtd_window <- rm_window <- mb_window <- mar_window <- vector()
for (i in 1:length(yr)){
  wtd_window[i] <- mean(wtd_output[(yr[i]-window/2):(yr[i]+window/2)])
  rm_window[i] <- mean(layer_remaining_mass[(yr[i]-window/2):(yr[i]+window/2)])
  mb_window[i] <- mean(mb[(yr[i]-window/2):(yr[i]+window/2)])
  mar_window[i] <- mean(mar[(yr[i]-window/2):(yr[i]+window/2)])
}

## plotting ##
par(mfrow=c(4,1),mar=c(4.5,5,1,2))
plot(NA,NA,xlim=x_scale,ylim=c(60,0),xlab="time yr",ylab="WTD cm")
lines(seq(1,t_extent,1),wtd_output,col="gray")
lines(yr,wtd_window,lwd=2)
plot(NA,NA,xlim=x_scale,ylim=c(1,0),xlab="time yr",ylab="remained mass")
lines(seq(1,t_extent,1),rm,col="gray")
lines(yr,rm_window,lwd=2)
plot(NA,NA,xlim=x_scale,ylim=c(200,0),xlab="time yr",ylab="mass balance g m-2 yr-1")
lines(seq(1,t_extent,1),mb,col="gray")
lines(yr,mb_window,lwd=2)
plot(NA,NA,xlim=x_scale,ylim=c(300,0),xlab="time yr",ylab="mass accmulation rate g m-2 yr-1")
lines(seq(1,t_extent,1),mar,col="gray")
lines(yr,mar_window,lwd=2)