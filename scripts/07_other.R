#################################
#Load data and libraries
#################################
source("scripts/04_combine-MIC-flow.R")

plot(all10.ag.0$data, fitFlow.sd$t10.G1.1)
cor.test(all10.ag.0$data, fitFlow.sd$t10.G1.1)

plot(all10.ag.0$sd, fitFlow.sd$t10.G1.1, xlab="Variation in evolved YPD (OD)", ylab = "Variation in G1 mean")
abline(lm(fitFlow.sd$t10.G1.1~all10.ag.0$sd))
cor.test(all10.ag.0$sd, fitFlow.sd$t10.G1.1)
