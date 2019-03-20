#################################
#Load libraries and functions
#################################
library(Hmisc)
library(lme4)
library(car)
library(agricolae)
library(lmerTest)
library(beanplot)
library(RColorBrewer)
library(nlme)
library(MASS)
library(dplyr)
options(stringsAsFactors = FALSE, max.print=5000)

cv <-function(x) sd(x,na.rm=TRUE)/mean(x,na.rm=TRUE)
se <- function(x, na.rm=TRUE) sqrt(var(x, na.rm=TRUE)/(length(x) - 1))
source("/Users/acgerstein/Documents/Postdoc/Research/0CommonFiles/Rscript/mvee.R")

#################################
#Load data
#################################
source("/Users/acgerstein/Documents/Postdoc/Research/Clinical_isolates/environments/mutAccum/allMA/Rcode/2017fall/180115setup-lowFLCEvol.R")

#################################
#Fitness in FLC1
#################################
#Figure 1
pdf("/Users/acgerstein/Documents/Postdoc/Papers/MutAccum/Figures/Figure1-OD1.pdf", width=7, height=5.5)
par(mfrow=c(2, 1),mar=c(1,1 , 1, 1), oma=c(3, 3, 1, 1))
plot(all10.1$place, all10.1$data, col=all10.1$col, ylim=c(0, 2), xaxt="n", yaxt="n", ann=F, cex=0.7)
points(1:20, all0.ag.1$data[place], pch="-", cex=2.25, col=grey(0.3))
points(1:20, all10.ag.1$data[place], pch="-", cex=2.25, col=colours)
axis(1, 1:20, labels=FALSE)
axis(2, las=2)
mtext("OD at 24h", side=3, adj=0.01)

plot(all10.1.72$place, all10.1.72$data, col=all10.1.72$col, ylim=c(0, 2), xaxt="n", yaxt="n", ann=F, cex=0.7)
points(1:20, all0.ag.1.72$data[place], pch="-", cex=2.25, col=grey(0.4))
points(1:20, all10.ag.1.72$data[place], pch="-", cex=2.25, col=colours)
axis(1, seq(1, 20, 2), paste0("A", all0.ag.1$line[place][seq(1, 20, 2)]), cex.axis=0.8)
axis(1, seq(2, 20, 2), paste0("A", all0.ag.1$line[place][seq(2, 20, 2)]), cex.axis=0.8)
axis(2, las=2)
mtext("strain", side=1, line=2)
mtext("OD in FLC1", side=2, line=1.5, outer=TRUE)
mtext("OD at 72h", side=3, adj=0.01)
dev.off()
system("open /Users/acgerstein/Documents/Postdoc/Papers/MutAccum/Figures/Figure1-OD1.pdf")

#Figure 2
pdf("/Users/acgerstein/Documents/Postdoc/Papers/MutAccum/Figures/Figure2-OD1cor.pdf", width=3, height=5.5)
par(mfrow=c(2, 1),mar=c(1,1 , 1, 1), oma=c(3, 3, 1, 1))
plot(fitFlow.ag$all0.ag.1, fitFlow.ag$change24ave, xaxt="n", yaxt="n", pch=19, col=fitFlow.ag$col, ann=F, xlim=c(0, 2))
axis(1, labels=FALSE)
axis(2, las=2)
mtext("Assessed at 24h", side=3, adj=0.01)
abline(lm(fitFlow.ag$change24ave~fitFlow.ag$all0.ag.1))
text(-0.1, -0.18, "cor = -0.58", pos=4)

plot(fitFlow.ag$all0.ag.1.72, fitFlow.ag$all10.ag.1.72-fitFlow.ag$all0.ag.1.72, xaxt="n", yaxt="n", pch=19, col=fitFlow.ag$col, ann=F, xlim=c(0, 2))
axis(1)
axis(2, las=2)
mtext("Assessed at 72h", side=3, adj=0.01)
abline(lm(fitFlow.ag$change72ave~fitFlow.ag$all0.ag.1.72))
mtext("initial fitness", side=1, outer=FALSE, line=2.5)
txt <- expression(paste("change in fitness in 1" ,mu, "g/mL fluconazole"))
mtext(txt, side=2, outer=TRUE, line=2)
text(-0.1, 0.02, "cor = -0.94", pos=4)
dev.off()
system("open /Users/acgerstein/Documents/Postdoc/Papers/MutAccum/Figures/Figure2-OD1cor.pdf")

#################################
#Fitness in YPD
#################################

#dashes
pdf("/Users/acgerstein/Documents/Postdoc/Papers/MutAccum/Figures/Figure3-YPD.pdf", width=7, height=5.5)
par(mfrow=c(2, 1), mar=c(1, 1, 1, 1), oma=c(3, 3, 1, 1))
plot(all10.0$place, all10.0$data, col=all10.0$col, ylim=c(0, 2), xaxt="n", yaxt="n", ann=F, cex=0.9)
points(1:20, all0.ag.0$data[place], col=grey(0.3), pch="-", cex=2.25)
points(1:20, all10.ag.0$data[place], col=colours, pch="-", cex=2.25)
axis(1, 1:20, labels=FALSE)
axis(2, las=2)
mtext("Growth ability at 24 h", side=3, line=2, outer=FALSE)

plot(all10.0.72$place, all10.0.72$data, ylim=c(0, 2), xaxt="n", yaxt="n", col=all10.0.72$col, ann=F, cex=0.9)
points(1:20, all10.ag.0.72$data[place], col=colours, pch="-", cex=2.25)
points(1:20, all0.ag.0.72$data[place], col=grey(0.3), pch="-", cex=2.25)
axis(1, seq(1, 20, 2), paste0("A", all0.ag.1$line[place][seq(1, 20, 2)]), cex.axis=0.8)
axis(1, seq(2, 20, 2), paste0("A", all0.ag.1$line[place][seq(2, 20, 2)]), cex.axis=0.8)
axis(2, las=2)
mtext("strain", side=1, line=2, outer=TRUE)
mtext("Growth ability in drug-free environment ", side=2, line=1.5, outer=TRUE)
mtext("Growth ability at 72 h", side=3, adj=0.01)
dev.off()
system("open /Users/acgerstein/Documents/Postdoc/Papers/MutAccum/Figures/Figure3-YPD.pdf")

pdf("/Users/acgerstein/Documents/Postdoc/Papers/MutAccum/Figures/FigureS1-YPD-FL1-cor.pdf", width=5.5, height=5.5)
plot(subset(all10.ag, enviro==0)$data, subset(all10.ag, enviro==1)$data, yaxt="n", xlim=c(0.4, 1.7), ylim=c(0.4, 1.7), xlab="Growth ability in drug-free environment", ylab="Growth ability in low drug environment", cex=1.2, pch=19, cex.lab=1.2)
axis(2, las=2)
abline(lm( subset(all10.ag, enviro==1)$data~ subset(all10.ag, enviro==0)$data))
dev.off()

#significant correlation among replicates though in multiple lines at 24h for YPD
YPDplace <- c(2, 3, 4, 6, 7, 8, 9, 13, 14, 16, 17, 20, 10, 11, 15, 19, 1, 5, 12, 18)
pdf("/Users/acgerstein/Documents/Postdoc/Papers/MutAccum/Figures/FigureS2-YPD-FLC1-24h-cor-sigOrder.pdf", width=7, height=7)
par(mfrow=c(5, 4), mar=c(0.5, 0.5, 0.5, 0.5), oma=c(4, 4, 1, 1))
j <- 0
for(i in YPDplace){
  j <- j+1
  plot(subset(all10.0, line==i)$data, subset(all10.1, line==i)$data, xlim=c(0, 2), ylim=c(0, 2), xaxt="n", yaxt="n", col="black")
  if(j < 17){
  test <- cor.test(subset(all10.0, line==i)$data, subset(all10.1, line==i)$data, method="spearman")
  if(test$p.value < 0.05) abline(lm(subset(all10.1, line==i)$data~subset(all10.0, line==i)$data), col="red", lwd=1.5)
  else abline(lm(subset(all10.1, line==i)$data~subset(all10.0, line==i)$data), lty=2, col="red")
}
  text(0.15, 1.8, paste0("A", i), adj=0.01, font=2, cex=1.2)
  if(j %%4==1) axis(2, las=2, at=c(0, 0.5, 1,1.5, 2), labels=c(0, "0.5", 1, "1.5", 2))
  else axis(2, labels=FALSE)
  if(j > 16) axis(1, at=c(0, 0.5, 1,1.5, 2), labels=c(0, "0.5", 1, "1.5", 2))
  else axis(1, labels=FALSE)
}
mtext("Evolved growth ability in drug-free enviornment", side=1, outer=TRUE, line=2)
mtext("Evolved growth ability in low drug environment", side=2, outer=TRUE, line=2)
dev.off()

####################
#resistance
####################

pdf("/Users/acgerstein/Documents/Postdoc/Papers/MutAccum/Figures/Figure3-BMD_resist.pdf", width=6.5, height=4)
plot(1:20, log(fitFlow.plot.ag$MIC24)[place], ylim=c(log(0.125), log(256)), yaxt="n", xaxt="n", xlab="strain", ylab="", pch="-", cex=2.25, col=grey(0.3))
points(jitter(fitFlow.plot$place), log(fitFlow.plot$MIC24.10), pch=21, col=fitFlow.plot$col)
axis(2, at=c(log(0.25),log(1), log(4), log(16), log(64), log(256)), labels=c("<1", "1", "4", "16", "64", ">128"), las=2)
axis(1, 1:20, labels=FALSE, cex.axis=0.5)
abline(h=log(8), lty=2)
mtext(expression(MIC[50]), side=2, line=2)
axis(1, 1:20, labels=FALSE)
text(1:20, -3, paste0("A", fitFlow.plot.ag$line[place]), srt=-45, adj=0.1, xpd=NA)
dev.off()
system("open /Users/acgerstein/Documents/Postdoc/Papers/MutAccum/Figures/Figure3-BMD_resist.pdf")

pdf("/Users/acgerstein/Documents/Postdoc/Papers/MutAccum/Figures/FigureS1-tolerance.pdf", width=6, height=5.5)
par(mfrow=c(3, 1), mar=c(1, 1, 1, 1), oma=c(3, 4, 1, 1), mgp=c(1.5, 0.75, 0))
plot(1:20, fitFlow.ag$SMG24.2[place], ylim=c(0, 1), yaxt="n", xaxt="n", xlab="strain", ylab="", pch="-", cex=2.25, col=grey(0.3))
points(jitter(fitFlow$place), fitFlow$SMG24.10.2, pch=21, col=fitFlow$col, cex=1.2)
points(1:20, fitFlow.ag$SMG24.up[place], pch="-", cex=2.25, col=grey(0.5))
points(1:20, fitFlow.ag$SMG24.down[place], pch="-", cex=2.25, col=grey(0.5))
points(1:20, fitFlow.ag$SMG24.2[place], pch="-", cex=2.25, col=grey(0.3))
axis(2, las=2)
axis(1, 1:20, labels=FALSE, cex.axis=0.5)
mtext("a 24 h" , side=3, adj=0.01, font=2, cex=1)

plot(1:20, fitFlow.ag$SMG48.2[place], ylim=c(0, 1), yaxt="n", xaxt="n", xlab="strain", ylab="", pch="-", cex=2.25, col=grey(0.3))
points(jitter(fitFlow$place), fitFlow$SMG48.10.2, pch=21, col=fitFlow$col, cex=1.2)
points(1:20, fitFlow.ag$SMG48.up[place], pch="-", cex=2.25, col=grey(0.5))
points(1:20, fitFlow.ag$SMG48.down[place], pch="-", cex=2.25, col=grey(0.5))
points(1:20, fitFlow.ag$SMG48.2[place], pch="-", cex=2.25, col=grey(0.3))
axis(2, las=2)
axis(1, 1:20, labels=FALSE, cex.axis=0.5)
mtext("b 48 h" , side=3, adj=0.01, font=2, cex=1)

plot(1:20, fitFlow.ag$SMG72.2[place], ylim=c(0, 1), yaxt="n", xaxt="n", xlab="strain", ylab="", pch="-", cex=2.25, col=grey(0.3))
points(jitter(fitFlow$place), fitFlow$SMG72.10.2, pch=21, col=fitFlow.plot$col, cex=1.2)
points(1:20, fitFlow.ag$SMG72.up[place], pch="-", cex=2.25, col=grey(0.5))
points(1:20, fitFlow.ag$SMG72.down[place], pch="-", cex=2.25, col=grey(0.5))
points(1:20, fitFlow.ag$SMG72.2[place], pch="-", cex=2.25, col=grey(0.3))
axis(2, las=2)
axis(1, 1:20, labels=FALSE, cex.axis=0.5)
mtext("c 72 h" , side=3, adj=0.01, font=2, cex=1)
mtext("evolved tolerance" , side=2, outer=TRUE, line=2)
text(1:20, -0.15, paste0("A", fitFlow.plot.ag$line[place]), srt=-45, adj=0.1, xpd=NA)
dev.off()
system("open /Users/acgerstein/Documents/Postdoc/Papers/MutAccum/Figures/FigureS1-tolerance.pdf")

SMG24up <- subset(fitFlow, SMG24.10.2 > rep(SMG24.up, each=12)) #72 (13 strains)
SMG24down <- subset(fitFlow, SMG24.10.2 < rep(SMG24.down, each=12)) #79 (17 strains)

table(SMG72down$strain)
# 1  2 3 4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
# 1  8  0 8  1  4  1  5  4  6  1  1  6  9  8  9  7  6  2  4

SMG.mat.down <- data.frame(line = 1:20, fake = rev(c(11, 4, 12, 4, 11, 8, 11, 7, 8, 6, 11, 11, 6, 3, 4, 3, 5, 6, 10, 8)), increase = rev(c(1,  8, 0, 8, 1, 4, 1, 5, 4, 6, 1, 1, 6, 9, 8, 9, 7, 6, 2, 4)))




table(SMG72up$strain)
# 1 2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
# 3 0 12  3  1  2  4  2  6  2 11  8  3 0 2  0 2  0 7  5

SMG.mat.up <- data.frame(line = 1:20,  increase = rev(c( 4, 0, 12,  3,  1,  2,  4,  2,  6,  2, 11,  9,  3, 0, 2,  0, 2,  0, 7,  5)), fake = rev(c(8, 12, 0, 9, 11, 10, 8, 10, 6, 10, 1, 3, 9, 12, 10, 12, 10, 12, 5, 7)))

pdf("/Users/acgerstein/Documents/Postdoc/Papers/MutAccum/Figures/Figure4-SMG72_upDown.pdf", width=7.5, height=5.5)
par(mfrow=c(1, 2), mar=c(1, 1, 1, 2), oma=c(3, 3, 1, 1))
mp <- barplot(as.matrix(t(data.frame(rev(SMG.mat.down[,2][place]), rev(SMG.mat.down[,3][place])))), horiz=T, beside=FALSE, col=c("white", "darkred"), xlab="", yaxt="n", xaxt="n")
axis(1, at=0:12, 12:0, cex.axis=0.705, pos=0)
mtext("decreased tolerance", side=3, cex=0.9, line=-0.5)
axis(4, at=mp,  las=2, adj=0.5, labels=FALSE)
text(13.35,  mp, paste0("A", rev(place)), adj=0.5, xpd=NA)

mp <- barplot(as.matrix(t(data.frame(rev(SMG.mat.up[,2][place]), rev(SMG.mat.up[,3][place])))), horiz=T, beside=FALSE, col=c("navyblue", "white"), xlab="", yaxt="n", xaxt="n")
axis(1, at=0:12, 0:12, cex.axis=0.705, pos=0)
mtext("increased tolerance", side=3, cex=0.9, line=-0.5)
mtext("Number of replicates", side=1, cex=0.9, line=0.5, outer=TRUE)
axis(2, at=mp, labels=FALSE, pos=0)
dev.off()
system("open /Users/acgerstein/Documents/Postdoc/Papers/MutAccum/Figures/Figure4-SMG72_upDown.pdf")

################
#Ploidy
################
pdf("/Users/acgerstein/Documents/Postdoc/Papers/MutAccum/Figures/Figure6-genomeSize.pdf", width=3, height=6)
par(mfrow=c(2, 1), mar=c(1, 1, 1, 1), oma=c(3, 4, 1, 1))
plot(c(rep(1, 20), rep(2, 20)), c(flow.ag$t0.G1.mu, flow.ag$t10.G1.1), xlim=c(0.75, 2.25), xaxt="n", yaxt="n", ylim=c(140, 260))
axis(2, las=2)
for(i in 1:20){
  if (i %in% c(4, 9, 12))  lines(c(1, 2), c(flow.ag$t0.G1.mu[i], flow.ag$t10.G1.1[i]), lty=2)
  else lines(c(1, 2), c(flow.ag$t0.G1.mu[i], flow.ag$t10.G1.1[i]))
}
axis(1, c(1, 2), labels=FALSE)
mtext("median genome size \n (FITC intensity)", side=2, line=2.75)
mtext("a" , side=3, adj=0.01, font=2, cex=1)
mtext(" Genome size", side=3, adj=0.1)


plot(c(rep(1, 20), rep(2, 20)), c(flow.ag.cv$t0.G1.mu, flow.ag.cv$t10.G1.1), xlim=c(0.75, 2.25), xaxt="n", yaxt="n", ylim=c(0, 0.4))
axis(2, las=2)
for(i in 1:20){
  if (i %in% c(1, 5, 12, 18))  lines(c(1, 2), c(flow.ag.cv$t0.G1.mu[i], flow.ag.cv$t10.G1.1[i]), lty=2)
  else lines(c(1, 2), c(flow.ag.cv$t0.G1.mu[i], flow.ag.cv$t10.G1.1[i]))
}
axis(1, c(1, 2), labels=FALSE)
axis(1, c(1, 2), labels=c("ancestral", "evolved"))
mtext("b" , side=3, adj=0.01, font=2, cex=1)
mtext("Coefficient of variation \n (FITC intensity)", side=2, line=2.75)
mtext(" Genome size variation", side=3, adj=0.2)
dev.off()

par(mfrow=c(2, 1), mar=c(1, 1, 1, 1), oma=c(3, 4, 1, 1))
plot(fitFlow.ag$delta1.72, fitFlow.ag$t10.G1.1, xlab="", ylab="", yaxt="n", xaxt="n", pch=19, col=fitFlow.ag$col, cex=1.25)
axis(2, las=2)
axis(1, labels=FALSE)
mtext("evolved genome size", side=2, line=3)
abline(lm(fitFlow.ag$t10.G1.1~fitFlow.ag$delta1.72), lty=2)

plot(fitFlow.ag$delta1.72, fitFlow.ag$cv.t10.G1.1, xlab="", ylab="", yaxt="n", xaxt="n", pch=19, col=fitFlow.ag$col, cex=1.25)
axis(2, las=2)
mtext("evolved genome size variation", side=2, line=3)
abline(lm(fitFlow.ag$cv.t10.G1.1~fitFlow.ag$delta1.72), lty=2)
axis(1)
mtext("change in fitness", side=1, line=2)



#Evolved ploidy * evolved SMG72
pdf("/Users/acgerstein/Documents/Postdoc/Papers/MutAccum/Figures/FigureS3SMG72-evolPloidy.pdf", width=7, height=5)
k<-0
par(mar=c(1, 1, 1, 1), mfrow=c(4, 4), mar=c(1, 1, 1, 1), oma=c(4, 4, 1, 1))
for(i in c(place[1:16])){
  k <- k+1
plot(subset(fitFlow, line==i)$SMG72.10.2, subset(fitFlow, line==i)$t10.G1.1, ylim=c(100, 500), xlim=c(0, 2), xaxt="n", yaxt="n", pch=subset(fitFlow, line==i)$pch)
t <- cor.test(subset(fitFlow, line==i)$SMG72.10.2, subset(fitFlow, line==i)$t10.G1.1, method="spearman")
if(t$p.value < 0.05){
   abline(lm(subset(fitFlow, line==i)$t10.G1.1~subset(fitFlow, line==i)$SMG72.10.2))
   print(t)
 }
if(k > 12) axis(1)
else axis(1, labels=FALSE)
if(k %% 4 == 1) axis(2, las=2)
else axis(2, labels=FALSE)
text(2, 450, paste0("A", i), pos=2, font=2, cex=1.1)
}
mtext("Evolved residual growth (SMG72)", side=1, outer=TRUE, line=2)
mtext("Evolved genome size (FITC intensity)", side=2, outer=TRUE, line=2)
dev.off()

#######HERE#######

for(i in c(2:4, 6:11, 13:17, 19, 20)) {
  print(i)
print(cor.test(subset(fitFlow, line==i)$data, subset(fitFlow, line==i)$flow.t10.G1.1, method="spearman"))
}


pdf("/Users/acgerstein/Documents/Postdoc/Papers/MutAccum/Figures/FigureS3-evolFit-evolPloidy-72h-SMG72.pdf", width=7, height=5)
k<-0
par(mar=c(1, 1, 1, 1), mfrow=c(4, 4), mar=c(1, 1, 1, 1), oma=c(4, 4, 1, 1))
for(i in c(place[1:16])){
  k <- k+1
  plot(subset(fitFlow, line==i)$all10.ag.1.72, subset(fitFlow, line==i)$t10.G1.1, ylim=c(100, 500), xlim=c(0, 2), xaxt="n", yaxt="n", pch=subset(fitFlow, line==i)$pch, bg=ifelse(subset(fitFlow, line ==i)$change72.SMG >0, "black", "white"))
  t <- cor.test(subset(fitFlow, line==i)$all10.ag.1.72, subset(fitFlow, line==i)$t10.G1.1, method="spearman")
  if(t$p.value < 0.05) abline(lm(subset(fitFlow, line==i)$t10.G1.1~subset(fitFlow, line==i)$all10.ag.1.72))
  if(k > 12) axis(1)
  else axis(1, labels=FALSE)
  if(k %% 4 == 1) axis(2, las=2)
  else axis(2, labels=FALSE)
  text(0, 450, paste0("A", i), pos=4, font=2, cex=1.1)
  }
mtext("Evolved fitness in low drug (72 h OD)", side=1, outer=TRUE, line=2)
mtext("Evolved genome size (FITC intensity)", side=2, outer=TRUE, line=2)
dev.off()

for(i in c(2:4, 6:11, 13:17, 19, 20)) {
  print(i)
print(cor.test(subset(fitFlow, line==i)$data.1, subset(fitFlow, line==i)$flow.t10.G1.1, method="spearman"))
}


pdf("/Users/acgerstein/Documents/Postdoc/Papers/MutAccum/Figures/FigureS4-SMG24-evolPloidy.pdf", width=7, height=5)
k<-0
par(mar=c(1, 1, 1, 1), mfrow=c(4, 4), mar=c(1, 1, 1, 1), oma=c(4, 4, 1, 1))
for(i in c(place[1:16])){
  k <- k+1
  plot(subset(fitFlow.variable, line==i)$SMG24.10, subset(fitFlow.variable, line==i)$t10.G1.1, ylim=c(100, 500), xlim=c(0, 1), xaxt="n", yaxt="n")
  t <- cor.test(subset(fitFlow.variable, line==i)$SMG24.10, subset(fitFlow.variable, line==i)$t10.G1.1, method="spearman")
  if(t$p.value < 0.05) abline(lm(subset(fitFlow.variable, line==i)$t10.G1.1~subset(fitFlow.variable, line==i)$SMG24.10))
  if(k > 12) axis(1)
  else axis(1, labels=FALSE)
  if(k %% 4 == 1) axis(2, las=2)
  else axis(2, labels=FALSE)
  text(0, 450, paste0("A", i), pos=4, font=2, cex=1.1)
  }
mtext("Evolved SMG24", side=1, outer=TRUE, line=2)
mtext("Evolved genome size (FITC intensity)", side=2, outer=TRUE, line=2)
dev.off()

for(i in c(2:4, 6:11, 13:17, 19, 20)) {
  print(i)
print(cor.test(subset(fitFlow.variable, line==i)$SMG24.10, subset(fitFlow.variable, line==i)$t10.G1.1, method="spearman"))
}



pdf("/Users/acgerstein/Documents/Postdoc/Papers/MutAccum/Figures/FigureS4-SMG48-evolPloidy.pdf", width=7, height=5)
k<-0
par(mar=c(1, 1, 1, 1), mfrow=c(4, 4), mar=c(1, 1, 1, 1), oma=c(4, 4, 1, 1))
for(i in c(place[1:16])){
  k <- k+1
  plot(subset(fitFlow.variable, line==i)$SMG48.10, subset(fitFlow.variable, line==i)$t10.G1.1, ylim=c(100, 500), xlim=c(0, 1), xaxt="n", yaxt="n")
  t <- cor.test(subset(fitFlow.variable, line==i)$SMG48.10, subset(fitFlow.variable, line==i)$t10.G1.1, method="spearman")
  if(t$p.value < 0.05) abline(lm(subset(fitFlow.variable, line==i)$t10.G1.1~subset(fitFlow.variable, line==i)$SMG48.10))
  if(k > 12) axis(1)
  else axis(1, labels=FALSE)
  if(k %% 4 == 1) axis(2, las=2)
  else axis(2, labels=FALSE)
  text(0, 450, paste0("A", i), pos=4, font=2, cex=1.1)
  }
mtext("Evolved SMG24", side=1, outer=TRUE, line=2)
mtext("Evolved genome size (FITC intensity)", side=2, outer=TRUE, line=2)
dev.off()

for(i in c(2:4, 6:11, 13:17, 19, 20)) {
  print(i)
print(cor.test(subset(fitFlow.variable, line==i)$SMG48.10, subset(fitFlow.variable, line==i )$t10.G1.1, method="spearman"))
}


pdf("/Users/acgerstein/Documents/Postdoc/Papers/MutAccum/Figures/FigureS4-SMG72-evolPloidy.pdf", width=7, height=5)
k<-0
par(mar=c(1, 1, 1, 1), mfrow=c(4, 4), mar=c(1, 1, 1, 1), oma=c(4, 4, 1, 1))
for(i in c(2:4, 6:11, 13:17, 19, 20)){
  k <- k+1
  plot(subset(fitFlow.variable, line==i)$SMG72.10, subset(fitFlow.variable, line==i)$t10.G1.1, ylim=c(100, 500), xlim=c(0, 1), xaxt="n", yaxt="n")
  if(i !=3){
  t <- cor.test(subset(fitFlow.variable, line==i)$SMG72.10, subset(fitFlow.variable, line==i)$t10.G1.1, method="spearman")
  if(t$p.value < 0.05) abline(lm(subset(fitFlow.variable, line==i)$t10.G1.1~subset(fitFlow.variable, line==i)$SMG72.10))
}
  if(k > 12) axis(1)
  else axis(1, labels=FALSE)
  if(k %% 4 == 1) axis(2, las=2)
  else axis(2, labels=FALSE)
  text(0, 450, paste0("A", i), pos=4, font=2, cex=1.1)
  }
mtext("Evolved FoG80", side=1, outer=TRUE, line=2)
mtext("Evolved genome size (FITC intensity)", side=2, outer=TRUE, line=2)
dev.off()

for(i in c(2:4, 6:11, 13:17, 19, 20)) {
  print(i)
print(cor.test(subset(fitFlow.variable, line==i)$SMG72.10, subset(fitFlow.variable, line==i)$t10.G1.1, method="spearman"))
}
