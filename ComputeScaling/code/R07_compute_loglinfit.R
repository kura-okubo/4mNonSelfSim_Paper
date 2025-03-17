rm(list=ls())

setwd('/Users/kokubo/Dropbox/NIED_RESEARCH/4mBIAX_submission/4mNonSelfSim_Paper/ComputeScaling')
# imported file name
finame <- "./data/07_loglinearfit/logfitdata_fb03-087_G3_wlv_0.30_denoisemethod_detrend.csv"
mydata <- read.csv(finame)
mydata$X <- as.numeric(mydata$X)
mydata$Y <- as.numeric(mydata$Y)

quartz()

library(lmodel2)
fit_ma = lmodel2(formula=Y~X, data=mydata, range.y= "interval", range.x = "interval", nperm = 99)

lc <- c("black" ,"red", "blue", "green")

plot(fit_ma, "OLS", confidence=TRUE, xlab="log M0", ylab="log Tw", main = "", centr=FALSE, col=lc[1])
lines(fit_ma, "MA",  confidence=TRUE, centr=FALSE, col=c(lc[2], "pink"))
lines(fit_ma, "SMA", centr=FALSE, col=lc[3])
lines(fit_ma, "RMA", centr=FALSE, col=lc[4])

legend("topleft", c("OLS", "MA", "SMA", "RMA"), col=lc, lty=1)

# save the figure
quartz.save("./figure/07_loglinearfit/R_fit.png", dpi=150)

# dump the fitting result
sink("./data/07_loglinearfit/lmodel2_out_regression.txt")  
print(fit_ma$regression)
sink()

sink("./data/07_loglinearfit/lmodel2_out_confidence.txt")  
print(fit_ma$confidence)
sink()


