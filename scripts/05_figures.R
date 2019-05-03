#################################
#Load data and libraries
#################################
source("scripts/04_combine-MIC-flow.R")

#################################
#Fitness in FLC1
#################################
#Figure 1
pdf("manuscript/figures/Figure1-OD1.pdf", width=7, height=5.5)
par(mfrow=c(2, 1),mar=c(1,1 , 1, 1), oma=c(3, 3, 1, 1))
beeswarm(all10.1~place72, data = all10.1, corral = "wrap", col=coloursVa, ylim=c(0, 2), xaxt="n", yaxt="n", ann=F, cex=0.7, pch=19)
points(1:20, all0.ag.1$data[place72], pch="-", cex=2.25, col=grey(0.3))
points(1:20, all10.ag.1$data[place72], pch="-", cex=2.25, col=coloursV)
axis(1, 1:20, labels=FALSE)
axis(2, las=2)
mtext("Fitness at 24h", side=3, adj=0.01)

beeswarm(all10.1.72~place72, data = all10.1.72, corral = "wrap", col=coloursVa, ylim=c(0, 2), xaxt="n", yaxt="n", ann=F, cex=0.7, pch=19)
points(1:20, all0.ag.1.72$data[place72], pch="-", cex=2.25, col=grey(0.4))
points(1:20, all10.ag.1.72$data[place72], pch="-", cex=2.25, col=coloursV)
axis(1, seq(1, 20, 2), paste0("A", all0.ag.1$line[place72][seq(1, 20, 2)]), cex.axis=0.8)
axis(1, seq(2, 20, 2), paste0("A", all0.ag.1$line[place72][seq(2, 20, 2)]), cex.axis=0.8)
axis(2, las=2)
mtext("strain", side=1, line=2)
mtext(expression("Growth in YPD + 1 " ~ mu~ "M FLC (optical density)"), side=2, line=1.5, outer=TRUE)
mtext("Fitness at 72h", side=3, adj=0.01)

xl <- 1
yb <- 1
xr <- 1.5
yt <- 2
b <- coloursV[seq(1, length(coloursV), 2)]

par(mar=c(2.5,28,5,1))
par(new=TRUE)
plot(NA,type="n",ann=FALSE,xlim=c(1,2),ylim=c(1,2),xaxt="n",yaxt="n",bty="n")
rect(
  xl,
  head(seq(yb,yt,(yt-yb)/20),-1),
  xr,
  tail(seq(yb,yt,(yt-yb)/20),-1),
  col=coloursV, border=NA
)

mtext(expression("low anc.\n fitness"),side=1,cex=0.6, adj=-1, line=0.25)
mtext(expression("high anc.\n fitness"),side=3,cex=0.6, adj=-10, line=-0.25)

dev.off()
system("open manuscript/figures/Figure1-OD1.pdf")

#Figure 2
fitFlow.ag$change24ave <- fitFlow.ag$all10.1 - fitFlow.ag$all0.1
fitFlow.ag$change72ave <- fitFlow.ag$all10.1.72 - fitFlow.ag$all0.1.72

pdf("manuscript/figures/Figure2-OD1cor-mean_SD.pdf", width=7, height=6)
par(mfrow=c(2, 2),mar=c(1,1 , 1, 1), oma=c(3, 4, 1, 1))
plot(fitFlow.ag$all0.1, fitFlow.ag$change24ave, xaxt="n", yaxt="n", pch=19, col=fitFlow.ag$col72, ann=F, xlim=c(0, 2), ylim=c(-0.2, 1), cex=1.5)
axis(1, labels=FALSE)
axis(2, las=2)
mtext("Assessed at 24h", side=3, adj=0.01)
abline(lm(fitFlow.ag$change24ave~fitFlow.ag$all0.1))
text(-0.1, -0.18, "cor = -0.58", pos=4, cex=1.25)
txt <- expression(paste(Delta," mean fitness"))
mtext(txt, side=2, outer=FALSE, line=2.5)

par(mar=c(1,1 , 1, 1))
plot(fitFlow.ag$all0.1.72, fitFlow.ag$all10.1.72-fitFlow.ag$all0.1.72, xaxt="n", yaxt="n", pch=19, col=fitFlow.ag$col, ann=F, xlim=c(0, 2), ylim=c(-0.2, 1), cex=1.5)
axis(1, labels=FALSE)
axis(2, labels=FALSE)
mtext("Assessed at 72h", side=3, adj=0.01)
abline(lm(fitFlow.ag$change72ave~fitFlow.ag$all0.1.72))
text(-0.1, -0.18, "cor = -0.94", pos=4, cex=1.25)

plot(fitFlow.ag$all0.1, fitFlow.sd$all10.1-fitFlow.sd$all0.1, xaxt="n", yaxt="n", pch=19, col=fitFlow.ag$col, ann=F, xlim=c(0, 2), ylim=c(-0.1, 0.45), cex=1.5)
axis(1, labels=FALSE)
axis(2, las=2)
abline(lm(c(fitFlow.sd$all10.1-fitFlow.sd$all0.1)~fitFlow.ag$all0.1))
axis(2, las=2)
axis(1)
txt <- expression(paste(Delta," replicate variation"))
mtext(txt, side=2, outer=FALSE, line=2.5)
text(-0.1, -0.08, "cor = -0.77", pos=4, cex=1.25)

plot(fitFlow.ag$all0.1.72, fitFlow.sd$all10.1.72-fitFlow.sd$all0.1.72, xaxt="n", yaxt="n", pch=19, col=fitFlow.ag$col, ann=F, xlim=c(0, 2), ylim=c(-0.1, 0.45), cex=1.5)
axis(1)
axis(2, labels=FALSE)
abline(lm(c(fitFlow.sd$all10.1-fitFlow.sd$all0.1)~fitFlow.ag$all0.1))
text(-0.1, -0.08, "cor = -0.47", pos=4, cex=1.25)
mtext("Initial fitness", side=1, outer=TRUE, line=1.5)
dev.off()
system("open manuscript/figures/Figure2-OD1cor-mean_SD.pdf")


#################################
# Fitness in YPD
# Figure 3
#################################
#dashes
pdf("manuscript/figures/Figure3-YPD.pdf", width=7, height=5.5)
par(mfrow=c(2, 1), mar=c(1, 1, 1, 1), oma=c(3, 3, 1, 1))
beeswarm(all10.0~place72, data = all10.0, corral = "wrap", col=coloursVa, ylim=c(0, 2), xaxt="n", yaxt="n", ann=F, cex=0.7, pch=19)
points(1:20, all0.ag.0$data[place72], col=grey(0.3), pch="-", cex=2.25)
points(1:20, all10.ag.0$data[place72], col=coloursV, pch="-", cex=2.25)
axis(1, 1:20, labels=FALSE)
axis(2, las=2)
mtext("Growth ability at 24 h", side=3, line=2, outer=FALSE)

beeswarm(all10.0.72~place72, data = all10.0.72, corral = "wrap", col=coloursVa, ylim=c(0, 2), xaxt="n", yaxt="n", ann=F, cex=0.7, pch=19)
points(1:20, all10.ag.0.72$data[place72], col=coloursV, pch="-", cex=2.25)
points(1:20, all0.ag.0.72$data[place72], col=grey(0.3), pch="-", cex=2.25)
axis(1, seq(1, 20, 2), paste0("A", all0.ag.1$line[place72][seq(1, 20, 2)]), cex.axis=0.8)
axis(1, seq(2, 20, 2), paste0("A", all0.ag.1$line[place72][seq(2, 20, 2)]), cex.axis=0.8)
axis(2, las=2)
mtext("strain", side=1, line=2, outer=TRUE)
mtext("Growth ability in drug-free environment (optical density)", side=2, line=1.5, outer=TRUE)
mtext("Growth ability at 72 h", side=3, adj=0.01)

xl <- 1
yb <- 1
xr <- 1.5
yt <- 2
b <- coloursV[seq(1, length(coloursV), 2)]

par(mar=c(2.5,28,5,1))
par(new=TRUE)
plot(NA,type="n",ann=FALSE,xlim=c(1,2),ylim=c(1,2),xaxt="n",yaxt="n",bty="n")
rect(
  xl,
  head(seq(yb,yt,(yt-yb)/20),-1),
  xr,
  tail(seq(yb,yt,(yt-yb)/20),-1),
  col=coloursV, border=NA
)

mtext(expression("low anc.\n fitness"),side=1,cex=0.6, adj=-1, line=0.25)
mtext(expression("high anc.\n fitness"),side=3,cex=0.6, adj=-10, line=-0.25)
dev.off()
system("open manuscript/figures/Figure3-YPD.pdf")

pdf("manuscript/figures/FigureS1-YPD-FL1-cor.pdf", width=5.5, height=5.5)
plot(subset(all10.ag, enviro==0)$data, subset(all10.ag, enviro==1)$data, yaxt="n", xlim=c(0.4, 1.7), ylim=c(0.4, 1.7), xlab="Growth ability in drug-free environment", ylab="Growth ability in low drug environment", cex=1.2, pch=19, cex.lab=1.2, col = fitFlow.ag$col)
axis(2, las=2)
abline(lm( subset(all10.ag, enviro==1)$data~ subset(all10.ag, enviro==0)$data))
dev.off()

#significant correlation among replicates though in multiple lines at 24h for YPD
YPDplace <- c(3, 4, 6, 7, 8, 9, 13, 14, 16, 17, 2, 10, 11, 15, 19, 20, 1, 5, 12, 18)

pdf("manuscript/figures/FigureS2-YPD-FLC1-24h-cor-sigOrder.pdf", width=7, height=7)
par(mfrow=c(5, 4), mar=c(0.5, 0.5, 0.5, 0.5), oma=c(4, 4, 1, 1))
j <- 0
for(i in YPDplace){
  j <- j+1
  plot(subset(all10.0, line==i)$all10.0, subset(all10.1, line==i)$all10.1, xlim=c(0, 2), ylim=c(0, 2), xaxt="n", yaxt="n", col="black")
  if(j < 17){
  test <- cor.test(subset(all10.0, line==i)$all10.0, subset(all10.1, line==i)$all10.1, method="spearman")
  if(test$p.value < 0.05) abline(lm(subset(all10.1, line==i)$all10.1~subset(all10.0, line==i)$all10.0), col="red", lwd=1.5)
  else abline(lm(subset(all10.1, line==i)$all10.1~subset(all10.0, line==i)$all10.0), lty=2, col="red")
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
# Figure 4
####################

fitFlow.ag$MIC24[fitFlow.ag$MIC24==0.000125] <- 0.25
fitFlow.ag$MIC24.10[fitFlow.ag$MIC24.10==0.000125] <- 0.25
fitFlow$MIC24[fitFlow$MIC24==0.000125] <- 0.25
fitFlow$MIC24.10[fitFlow$MIC24.10==0.000125] <- 0.25

MIC24.ag.less1 <- subset(fitFlow.ag, MIC24 < 1)
MIC24.ag.less1 <- MIC24.ag.less1[order(MIC24.ag.less1$MIC24, MIC24.ag.less1$place72),]
MIC24.ag.1 <- subset(fitFlow.ag, MIC24 == 1)
MIC24.ag.1 <- MIC24.ag.1[order(MIC24.ag.1$place72),]
MIC24.ag.more1 <- subset(fitFlow.ag, MIC24 > 1)
MIC24.ag.more1 <- MIC24.ag.more1[order(MIC24.ag.more1$MIC24),]

subset(fitFlow, MIC24.10 > 64) #1 line in A5, 5 in A12
fitFlow$MIC24.10[fitFlow$MIC24.10 > 64] <- 256

MIC24.less1 <- subset(fitFlow, line %in% MIC24.ag.less1$line)
MIC24.less1 <- MIC24.less1[order(MIC24.less1$MIC24, MIC24.less1$place72),]
MIC24.1 <- subset(fitFlow, line %in% MIC24.ag.1$line)
MIC24.1 <- MIC24.1[order(MIC24.1$place72),]
MIC24.more1 <- subset(fitFlow, line %in% MIC24.ag.more1$line)
MIC24.more1 <- MIC24.more1[order(MIC24.more1$MIC24),]

order72_MIC <- c(MIC24.ag.less1$line, MIC24.ag.1$line, MIC24.ag.more1$line)

fitFlow$tabMIC10[fitFlow.plot$MIC24.10 < 1] <- "<1"
fitFlow$tabMIC10[fitFlow.plot$MIC24.10 == 1] <- "=1"
fitFlow$tabMIC10[fitFlow.plot$MIC24.10 > 1] <- ">1"

highMIC <- subset(fitFlow, MIC24.10 >=8)
table(highMIC$strain)

pdf("manuscript/figures/Figure4-BMD_resist.pdf", width=7.5, height=4.5)
plot(1:11, log2(MIC24.ag.less1$MIC24), ylim=c(log2(0.25), log2(256)), yaxt="n", xaxt="n", xlab="Strain", ylab="", pch="-", cex=2.25, col=grey(0.3), xlim=c(1, 22))
points(13:16, log2(MIC24.ag.1$MIC24), pch="-", cex=2.25, col=grey(0.3))
points(18:22, log2(MIC24.ag.more1$MIC24), pch="-", cex=2.25, col=grey(0.3))
points(jitter(rep(1:11, each=12)), log2(MIC24.less1$MIC24.10), pch=19, col=MIC24.less1$col72a)
abline(v=12, lty=2)
points(jitter(rep(13:16, each=12)), log2(MIC24.1$MIC24.10), pch=19, col=MIC24.1$col72a)
abline(v=17, lty=2)
points(jitter(rep(18:22, each=12)), log2(MIC24.more1$MIC24.10), pch=19, col=MIC24.more1$col72a)
axis(2, at=c(log2(0.25),log2(1), log2(4), log2(16), log2(64), log2(256)), labels=c("< 0.5", "1", "4", "16", "64", ">128"), las=2)
#abline(h=log(8), lty=2)
mtext(expression(MIC[50]), side=2, line=2)
axis(1, c(1:11, 13:16, 18:22), labels=c(paste0("A", c(MIC24.ag.less1$line, MIC24.ag.1$line, MIC24.ag.more1$line))), cex.axis=0.8)
axis(1, c(2, 4, 6, 8, 10, 14, 16, 19, 21) , labels=c(paste0("A", c(MIC24.ag.less1$line, MIC24.ag.1$line, MIC24.ag.more1$line)[c(2, 4, 6,8, 10, 13, 15, 17, 19)])),  cex.axis=0.8)
dev.off()

#################################
#Tolerance above 1ug
#################################


# SMG*.2 is above MIC
# SMG*.3 is above 1ug

fitFlow$place72_MIC <- c(rep(16, 12), rep(12, 12), rep(3, 12), rep(10, 12), rep(17, 12), rep(13, 12), rep(14, 12), rep(11, 12), rep(4, 12), rep(15, 12), rep(2, 12), rep(20, 12), rep(6, 12), rep(5, 12), rep(7, 12), rep(9, 12), rep(1, 12), rep(18, 12), rep(8, 12), rep(19, 12))

pdf("manuscript/figures/FigureS3-tolerance-aboveMIC.pdf", width=7, height=6.5)
par(mfrow=c(3, 1), mar=c(1, 1, 1, 1), oma=c(3, 4, 1, 1), mgp=c(1.5, 0.75, 0))
beeswarm(SMG24.10.2~place72_MIC, data = fitFlow, ylim=c(0, 1), yaxt="n", xaxt="n", xlab="strain", ylab="", pch=19, cex=1.2, col=coloursVa, corral="wrap", xlim=c(1, 20))
#points(1:20, fitFlow.ag$SMG24.2[place72], pch="-", cex=2.25, col=grey(0.2))
points(1:20, SMG24.up.2[order72_MIC], pch="-", cex=2.25, col=grey(0.2))
points(1:20, SMG24.down.2[order72_MIC], pch="-", cex=2.25, col=grey(0.2))
arrows(1:20, SMG24.up.2[order72_MIC], 1:20, SMG24.down.2[order72_MIC], col=grey(0.2), length=0)
axis(2, las=2)
axis(1, 1:20, labels=FALSE, cex.axis=0.5)
mtext("a 24 h" , side=3, adj=0.01, font=2, cex=1)

beeswarm(SMG48.10.2~place72_MIC, data = fitFlow, ylim=c(0, 1), yaxt="n", xaxt="n", xlab="strain", ylab="", pch=19, cex=1.2, col=coloursVa, corral="wrap")
#points(1:20, fitFlow.ag$SMG48.2[place72], pch="-", cex=2.25, col=grey(0.2))
points(1:20, SMG48.up.2[order72_MIC], pch="-", cex=2.25, col=grey(0.2))
points(1:20, SMG48.down.2[order72_MIC], pch="-", cex=2.25, col=grey(0.2))
arrows(1:20, SMG48.up.2[order72_MIC], 1:20, SMG48.down.2[order72_MIC], col=grey(0.2), length=0)
axis(2, las=2)
axis(1, 1:20, labels=FALSE, cex.axis=0.5)
mtext("b 48 h" , side=3, adj=0.01, font=2, cex=1)

beeswarm(SMG72.10.2~place72_MIC, data = fitFlow, ylim=c(0, 1), yaxt="n", xaxt="n", xlab="strain", ylab="", pch=19, cex=1.2, col=coloursVa, corral="wrap")
#points(1:20, fitFlow.ag$SMG72.2[place72], pch="-", cex=2.25, col=grey(0.2))
points(1:20, SMG72.up.2[order72_MIC], pch="-", cex=2.25, col=grey(0.2))
points(1:20, SMG72.down.2[order72_MIC], pch="-", cex=2.25, col=grey(0.2))
arrows(1:20, SMG72.up.2[order72_MIC], 1:20, SMG72.down.2[order72_MIC], col=grey(0.2), length=0)

axis(2, las=2)
axis(1, 1:20, labels=FALSE, cex.axis=0.5)
mtext("c 72 h" , side=3, adj=0.01, font=2, cex=1)
mtext("Evolved tolerance" , side=2, outer=TRUE, line=2)
axis(1, seq(1, 20, 2), paste0("A", all0.ag.1$line[as.numeric(order72_MIC)][seq(1, 20, 2)]), cex.axis=0.9)
axis(1, seq(2, 20, 2), paste0("A", all0.ag.1$line[as.numeric(order72_MIC)][seq(2, 20, 2)]), cex.axis=0.9)
#text(1:20, -0.15, paste0("A", fitFlow.plot.ag$line[place]), srt=-45, adj=0.1, xpd=NA)
dev.off()

order <- c(13, 17, 11, 3, 15, 9, 14, 4, 16, 19, 10, 8, 20, 7, 2, 6, 12, 5, 1, 18)

#All replicates, SMG measured as above 1
pdf("manuscript/figures/Figure5-SMG72-above1oug-flip.pdf", width=7.5, height=4.5)
plot(rep(1:11, each =12), MIC24.less1$SMG72.10.3, yaxt="n", xaxt="n", xlab="Strain", ylab="", cex=1.2, pch=19, xlim=c(1, 22), col =MIC24.less1$col72a, ylim=c(0,1.7))
points(rep(13:16, each=12), MIC24.1$SMG72.10.3, cex=1.2, col=MIC24.1$col72a, pch=19)
points(rep(18:22, each=12), MIC24.more1$SMG72.10.3, cex=1.2, col=MIC24.more1$col72a, pch=19)
abline(v=12, lty=2)
abline(v=17, lty=2)

axis(1, c(1:11, 13:16, 18:22), labels=c(paste0("A", c(MIC24.ag.less1$line, MIC24.ag.1$line, MIC24.ag.more1$line))), cex.axis=0.8)
axis(1, c(2, 4, 6, 8, 10, 14, 16, 19, 21) , labels=c(paste0("A", c(MIC24.ag.less1$line, MIC24.ag.1$line, MIC24.ag.more1$line)[c(2, 4, 6,8, 10, 13, 15, 17, 19)])),  cex.axis=0.8)

points(1:11, MIC24.ag.less1$SMG72.up.3, pch="-", col="black", cex=2)
points(1:11, MIC24.ag.less1$SMG72.down.3, pch="-", col="black", cex=2)
arrows(1:11, MIC24.ag.less1$SMG72.up.3, 1:11, MIC24.ag.less1$SMG72.down.3, length=0, col="black")
points(12:15, MIC24.ag.1$SMG72.up.3, pch="-", col="black", cex=2)
points(12:15, MIC24.ag.1$SMG72.down.3, pch="-", col="black", cex=2)
arrows(12:15, MIC24.ag.1$SMG72.up.3, 12:15, MIC24.ag.1$SMG72.down.3, length=0, col="black")
points(16:20, MIC24.ag.more1$SMG72.up.3, pch="-", col="black", cex=2)
points(16:20, MIC24.ag.more1$SMG72.down.3, pch="-", col="black", cex=2)
arrows(16:20, MIC24.ag.more1$SMG72.up.3, 16:20, MIC24.ag.more1$SMG72.down.3, length=0, col="black")
axis(2, las=2)
mtext("Tolerance above FLC1 (72h)", side=2, line=3)
dev.off()

MIC24.less1_1 <- subset(MIC24.less1, MIC24.10 == 1)
MIC24.1_1 <- subset(MIC24.1, MIC24.10 == 1)
MIC24.more1_1 <- subset(MIC24.more1, MIC24.10 == 1)

table(MIC24.less1_1$line)
table(MIC24.1_1$line)
table(MIC24.more1_1$line)

#JUST strains with MIC = 1
pdf("manuscript/figures/Figure5-SMG72-aboveMIC-flip_MIC24-1.pdf", width=7.5, height=4.5)
plot(rep(1:11, c(11, 3,5, 5, 10, 12, 10, 8, 12, 12, 9)), MIC24.less1_1$SMG72.10.2, yaxt="n", xaxt="n", xlab="Strain", ylab="", cex=1.2, pch=19, xlim=c(1, 22), col = MIC24.less1_1$col72a, ylim=c(0,1))
points(rep(13:16, c(11, 12, 11, 8)), MIC24.1_1$SMG72.10.2, cex=1.2, col=MIC24.1_1$col72a, pch=19)
points(rep(18:22, c(0, 2, 0, 12, 6)), MIC24.more1_1$SMG72.10.2, cex=1.2, col=MIC24.more1_1$col72a, pch=19)
abline(v=12, lty=2)
abline(v=17, lty=2)

axis(1, c(1:11, 13:16, 18:22), labels=c(paste0("A", c(MIC24.ag.less1$line, MIC24.ag.1$line, MIC24.ag.more1$line))), cex.axis=0.8)
axis(1, c(2, 4, 6, 8, 10, 14, 16, 19, 21) , labels=c(paste0("A", c(MIC24.ag.less1$line, MIC24.ag.1$line, MIC24.ag.more1$line)[c(2, 4, 6,8, 10, 13, 15, 17, 19)])),  cex.axis=0.8)

points(1:11, MIC24.ag.less1$SMG72.up, pch="-", col="black", cex=2)
points(1:11, MIC24.ag.less1$SMG72.down, pch="-", col="black", cex=2)
arrows(1:11, MIC24.ag.less1$SMG72.up, 1:11, MIC24.ag.less1$SMG72.down, length=0, col="black")
points(13:16, MIC24.ag.1$SMG72.up, pch="-", col="black", cex=2)
points(13:16, MIC24.ag.1$SMG72.down, pch="-", col="black", cex=2)
arrows(13:16, MIC24.ag.1$SMG72.up, 13:16, MIC24.ag.1$SMG72.down, length=0, col="black")
points(18:22, MIC24.ag.more1$SMG72.up, pch="-", col="black", cex=2)
points(18:22, MIC24.ag.more1$SMG72.down, pch="-", col="black", cex=2)
arrows(18:22, MIC24.ag.more1$SMG72.up, 18:22, MIC24.ag.more1$SMG72.down, length=0, col="black")
axis(2, las=2)
mtext("Tolerance above MIC (72h)", side=2, line=3)
dev.off()


#All replicates, SMG above MIC
pdf("manuscript/figures/Figure5-SMG72-aboveMIC-flip.pdf", width=7.5, height=5)
plot(rep(1:11, each =12), MIC24.less1$SMG72.10.2, yaxt="n", xaxt="n", xlab="Strain", ylab="", cex=1.2, pch=19, xlim=c(1, 22), col =MIC24.less1$col72a, ylim=c(0,1))
points(rep(13:16, each=12), MIC24.1$SMG72.10.2, cex=1.2, col=MIC24.1$col72a, pch=19)
points(rep(18:22, each=12), MIC24.more1$SMG72.10.2, cex=1.2, col=MIC24.more1$col72a, pch=19)
abline(v=12, lty=2)
abline(v=17, lty=2)

axis(1, c(1:11, 13:16, 18:22), labels=c(paste0("A", c(MIC24.ag.less1$line, MIC24.ag.1$line, MIC24.ag.more1$line))), cex.axis=0.8)
axis(1, c(2, 4, 6, 8, 10, 14, 16, 19, 21) , labels=c(paste0("A", c(MIC24.ag.less1$line, MIC24.ag.1$line, MIC24.ag.more1$line)[c(2, 4, 6,8, 10, 13, 15, 17, 19)])),  cex.axis=0.8)

points(1:11, MIC24.ag.less1$SMG72.up, pch="-", col="black", cex=2)
points(1:11, MIC24.ag.less1$SMG72.down, pch="-", col="black", cex=2)
arrows(1:11, MIC24.ag.less1$SMG72.up, 1:11, MIC24.ag.less1$SMG72.down, length=0, col="black")
points(12:15, MIC24.ag.1$SMG72.up, pch="-", col="black", cex=2)
points(12:15, MIC24.ag.1$SMG72.down, pch="-", col="black", cex=2)
arrows(12:15, MIC24.ag.1$SMG72.up, 12:15, MIC24.ag.1$SMG72.down, length=0, col="black")
points(16:20, MIC24.ag.more1$SMG72.up, pch="-", col="black", cex=2)
points(16:20, MIC24.ag.more1$SMG72.down, pch="-", col="black", cex=2)
arrows(16:20, MIC24.ag.more1$SMG72.up, 16:20, MIC24.ag.more1$SMG72.down, length=0, col="black")
axis(2, las=2)
mtext("Tolerance above MIC (72h)", side=2, line=3)
dev.off()


pdf("manuscript/figures/Figure5-SMG72-above1ug.pdf", width=6, height=6)
layout(matrix(1:2,nrow=1),widths=c(0.8,0.1))
plot(subset(fitFlow, strain==17)$SMG72.10.3, rep(1, 12), ylim=c(0.5, 20.5), xlim=c(0, 1.8), yaxt="n", ylab="", xlab="", col=coloursVa[1], cex=1.2, pch=19)
#points(median(subset(fitFlow, strain==2)$SMG72.10.2, na.rm=TRUE), 1, pch="l", col=as.character(subset(fitFlow, strain==2)$col), cex=1.4)
points(fitFlow.ag$SMG72.up.3[17], 1, pch="|", col="black")
points(fitFlow.ag$SMG72.down.3[17], 1,  pch="|",  col="black")
arrows(fitFlow.ag$SMG72.down.3[17], 1, fitFlow.ag$SMG72.up.3[17], 1, length=0, col="black")
#points(fitFlow.ag$SMG72.2[2], 1, pch="|", col="black", cex=1.1)
k <- 1
for(i in as.numeric(order72_MIC)[2:20]){
  k <- k+1
  points(subset(fitFlow, strain==i)$SMG72.10.3, rep(k, 12), col=coloursVa[k], cex=1.2, pch=19)
  points(fitFlow.ag$SMG72.up.3[i], k, pch="|", col="black")
  points(fitFlow.ag$SMG72.down.3[i], k,  pch="|",  col="black")
  #points(fitFlow.ag$SMG72.2[i], k, pch="|", col="black", cex=1.1)
  arrows(fitFlow.ag$SMG72.down.3[i], k, fitFlow.ag$SMG72.up.3[i], k, length=0, col="black")
}
axis(2, labels=paste0("A", order72_MIC), at=1:20, las=2)
mtext("Tolerance (72h)", side=1, line=2)

xl <- 1
yb <- 1
xr <- 1.5
yt <- 2

par(mar=c(10,0,15,0))
plot(NA,type="n",ann=FALSE,xlim=c(1,2),ylim=c(1,2),xaxt="n",yaxt="n",bty="n")
rect(
  xl,
  head(seq(yb,yt,(yt-yb)/20),-1),
  xr,
  tail(seq(yb,yt,(yt-yb)/20),-1),
  col=coloursV, border=NA
)

mtext(expression("low anc.\n fitness"),side=1,cex=0.6, adj=0, line=0.25)
mtext(expression("high anc.\n fitness"),side=3,cex=0.6, adj=0, line=-0.25)

dev.off()
system("open manuscript/figures/Figure5-SMG72-above1ug.pdf")

################
#Ploidy
################
pdf("manuscript/figures/Figure6-genomeSize.pdf", width=3.25, height=7)
par(mfrow=c(2, 1), mar=c(1, 1, 1, 1), oma=c(3, 4, 1, 1))
plot(c(rep(1, 20), rep(2, 20)), c(flow.ag$t0.G1.mu, flow.ag$t10.G1.1), xlim=c(0.75, 2.25), xaxt="n", yaxt="n", ylim=c(140, 260), col=fitFlow.ag$col72, pch=19)
axis(2, las=2)
for(i in 1:20){
  if (i %in% c(4, 9, 12))  lines(c(1, 2), c(flow.ag$t0.G1.mu[i], flow.ag$t10.G1.1[i]), lty=2)
  else lines(c(1, 2), c(flow.ag$t0.G1.mu[i], flow.ag$t10.G1.1[i]))
}
points(c(rep(1, 20), rep(2, 20)), c(flow.ag$t0.G1.mu, flow.ag$t10.G1.1), col=fitFlow.ag$col72, pch=19)

axis(1, c(1, 2), labels=FALSE)
mtext("median genome size \n (FITC intensity)", side=2, line=2.75)
mtext("a" , side=3, adj=0.01, font=2, cex=1)
mtext("   Genome size", side=3, adj=0.1)


plot(c(rep(1, 20), rep(2, 20)), c(flow.ag.cv$t0.G1.mu, flow.ag.cv$t10.G1.1), xlim=c(0.75, 2.25), xaxt="n", yaxt="n", ylim=c(0, 0.4), col=fitFlow.ag$col72, pch=19)
axis(2, las=2)
for(i in 1:20){
  if (i %in% c(1, 5, 12, 18))  lines(c(1, 2), c(flow.ag.cv$t0.G1.mu[i], flow.ag.cv$t10.G1.1[i]), lty=2)
  else lines(c(1, 2), c(flow.ag.cv$t0.G1.mu[i], flow.ag.cv$t10.G1.1[i]))
}
points(c(rep(1, 20), rep(2, 20)), c(flow.ag.cv$t0.G1.mu, flow.ag.cv$t10.G1.1), col=fitFlow.ag$col72, pch=19)
axis(1, c(1, 2), labels=FALSE)
axis(1, c(1, 2), labels=c("ancestral", "evolved"))
mtext("b" , side=3, adj=0.01, font=2, cex=1)
mtext("Coefficient of variation \n (FITC intensity)", side=2, line=2.75)
mtext("   Genome size variation", side=3, adj=0.2)
dev.off()

#################################
#(non)correlations
#################################

pdf("manuscript/figures/Figure6-evolFit_72h-evolPloidy.pdf", width=7, height=6)
k<-0
par(mar=c(1, 1, 1, 1), mfrow=c(5, 4), mar=c(1, 1, 1, 1), oma=c(4, 4, 1, 1))
for(i in c(order72_MIC[1:20])){
  k <- k+1
  plot(subset(fitFlow, line==i)$all10.1.72, subset(fitFlow, line==i)$t10.G1.1, ylim=c(100, 500), xlim=c(1, 2), xaxt="n", yaxt="n", pch=subset(fitFlow, line==i)$pch, bg=ifelse(subset(fitFlow, line ==i)$MIC24.10 >1, "black", "white"))
  if(k < 16){
  t <- cor.test(subset(fitFlow, line==i)$all10.1.72, subset(fitFlow, line==i)$t10.G1.1, method="spearman")
  if(t$p.value < 0.05) abline(lm(subset(fitFlow, line==i)$t10.G1.1~subset(fitFlow, line==i)$all10.1.72), col="red")
  if(t$p.value > 0.05) abline(lm(subset(fitFlow, line==i)$t10.G1.1~subset(fitFlow, line==i)$all10.1.72), col="red", lty=2)
  }
  if(k > 12) axis(1)
  else axis(1, labels=FALSE)
  if(k %% 4 == 1) axis(2, las=2)
  else axis(2, labels=FALSE)
  text(1, 450, paste0("A", i), pos=4, font=2, cex=1.1)
}
mtext("Evolved growth  in low drug (OD at 72 h)", side=1, outer=TRUE, line=2)
mtext("Evolved genome size (FITC intensity)", side=2, outer=TRUE, line=2)
dev.off()


pdf("manuscript/figures/FigureSx-evolFit_72h-evolMIC.pdf", width=7, height=6)
k<-0
par(mar=c(1, 1, 1, 1), mfrow=c(5, 4), mar=c(1, 1, 1, 1), oma=c(4, 4, 1, 1))
for(i in c(order72_MIC[1:20])){
  k <- k+1
  plot(subset(fitFlow, line==i)$all10.1.72, log2(subset(fitFlow, line==i)$MIC24.10), xlim=c(1, 2), ylim=c(log2(0.125), log2(256)), xaxt="n", yaxt="n", pch=21)
  t <- cor.test(subset(fitFlow, line==i)$all10.1.72, log2(subset(fitFlow, line==i)$MIC24.10), method="spearman")
  if(!is.na(t$p.value)){
    if(t$p.value < 0.05) abline(lm(log2(subset(fitFlow, line==i)$MIC24.10)~subset(fitFlow, line==i)$all10.1.72), col="red")
    if(t$p.value > 0.05) abline(lm(log2(subset(fitFlow, line==i)$MIC24.10)~subset(fitFlow, line==i)$all10.1.72), col="red", lty=2)
  }
  if(k > 16) axis(1)
  else axis(1, labels=FALSE)
  if(k %% 4 == 1) axis(2, at=c(log2(0.25),log2(1), log2(4), log2(16), log2(64), log2(256)), labels=c("< 0.5", "1", "4", "16", "64", ">128"), las=2)
  else axis(2, labels=FALSE)
  text(1, 7, paste0("A", i), pos=4, font=2, cex=1.1)
}
mtext("Evolved growth  in low drug (OD at 72 h)", side=1, outer=TRUE, line=2)
mtext(expression(MIC[50]), side=2, line=2, outer = TRUE)
dev.off()

pdf("manuscript/figures/FigureSx-evolFit_72h-evolSMG.pdf", width=7, height=5)
k<-0
par(mar=c(1, 1, 1, 1), mfrow=c(5, 4), mar=c(1, 1, 1, 1), oma=c(4, 4, 1, 1))
for(i in c(order72_MIC[1:20])){
  k <- k+1
  plot(subset(fitFlow, line==i)$all10.1.72, subset(fitFlow, line==i)$SMG72.10.2, xlim=c(1, 2), ylim=c(0, 1), xaxt="n", yaxt="n", pch=21, bg=ifelse(subset(fitFlow, line ==i)$MIC24.10 ==1 , "blue", "white"))
  t <- cor.test(subset(fitFlow, line==i)$all10.1.72, subset(fitFlow, line==i)$SMG72.10.2, method="spearman")
  if(!is.na(t$p.value)){
    if(t$p.value < 0.05) abline(lm(subset(fitFlow, line==i)$SMG72.10.2~subset(fitFlow, line==i)$all10.1.72), col="red")
    if(t$p.value > 0.05) abline(lm(subset(fitFlow, line==i)$SMG72.10.2~subset(fitFlow, line==i)$all10.1.72), col="red", lty=2)
  }
  if(length(subset(fitFlow, MIC24.10 == 1 & line==i)$all10.1.72) > 0){
    if(i == 1) print(length(subset(fitFlow, MIC24.10 == 1 & line==i)$all10.1.72) > 0)
    t2 <- cor.test(subset(fitFlow, MIC24.10 == 1 & line==i)$all10.1.72, subset(fitFlow, MIC24.10 == 1 & line==i)$SMG72.10.2, method="spearman")
    if(!is.na(t2$p.value)){
      if(t2$p.value < 0.05) abline(lm(subset(fitFlow, MIC24.10 == 1 & line==i)$SMG72.10.2~subset(fitFlow, MIC24.10 == 1 & line==i)$all10.1.72), col="blue")
      if(t2$p.value > 0.05) abline(lm(subset(fitFlow, MIC24.10 == 1 & line==i)$SMG72.10.2~subset(fitFlow, MIC24.10 == 1 & line==i)$all10.1.72), col="blue", lty=2)
    }
  }
  if(k > 16) axis(1)
  else axis(1, labels=FALSE)
  if(k %% 4 == 1) axis(2, las=2)
  else axis(2, labels=FALSE)
  text(1, 0.8, paste0("A", i), pos=4, font=2, cex=1.1)
}
mtext("Evolved growth  in low drug (OD at 72 h)", side=1, outer=TRUE, line=2)
mtext("Tolerance above MIC (72 h)", side=2, line=2, outer = TRUE)
dev.off()


#For a resistance barplot (not used)
#fitFlow.tab <- data.frame(table(as.numeric(fitFlow.plot$line), fitFlow.plot$tabMIC10))
#names(fitFlow.tab) <- c("line", "MIC24.10", "Freq")

# fitFlow.tab.wide <- as_tibble(fitFlow.tab) %>%
#   spread(MIC24.10, Freq)
# fitFlow.tab.wide.flip <- t(fitFlow.tab.wide)
# fitFlow.tab.wide.flip[2,] <- as.numeric(fitFlow.tab.wide.flip[2,])
#
# pdf("manuscript/figures/Figure4-BMD_resist-bar.pdf", width=8.5, height=4.5)
# par(xpd=NA, oma=c(1, 1, 1, 8))
# m <- barplot(as.matrix(fitFlow.tab.wide.flip[2:4,as.numeric(order72_MIC)]), beside=FALSE, xaxt="n", yaxt="n", col=c("darkgrey", "purple", "orange"),  border=NA, space=c(rep(0.2, 11), 1.5, rep(0.2, 3), 1.5, rep(0.2, 5)), xlim=c(0.7, 26.1))
# axis(1, m[1:20][seq(1, 20, 2)], labels=FALSE)
# text(m[1:20][seq(1, 20, 2)], -1, paste0("A", order72_MIC[seq(1, 20, 2)]), cex=0.6)
# axis(1, m[1:20][seq(2, 21, 2)], labels = FALSE)
# text(m[1:20][seq(2, 21, 2)], -1, paste0("A", order72_MIC[seq(2, 21, 2)]), cex=0.6)
# axis(2, las=2)
# box()
# legend(29, 12, legend=c(expression(MIC[50] ~ "< 1", MIC[50] ~ "= 1", MIC[50] ~ "> 1" )), inset=0.05, pch=22, col=c("darkgrey", "purple", "orange"), pt.bg =c("darkgrey", "purple", "orange"), cex=0.7)
# mtext("number of evolved replicates", side=2, line=2)
# mtext(expression("Ancestral \n" ~ MIC[50] ~ "< 1"), side= 1, adj = 0.2, line = 3, cex= 0.9)
# mtext(expression("Ancestral \n" ~ MIC[50] ~ "= 1"), side= 1, adj = 0.65, line = 3, cex= 0.9)
# mtext(expression("Ancestral \n" ~ MIC[50] ~ "\n> 1"), side= 1, adj = 1, line=3, cex= 0.9)
# dev.off()
#

#fitFlow.ag$change24ave <- fitFlow.ag$all10.1 - fitFlow.ag$all0.1
#fitFlow.ag$change72ave <- fitFlow.ag$all10.1.72 - fitFlow.ag$all0.1.72

pdf("manuscript/figures/FigureSx-ancOD1cor-Ploidy_mean_SD.pdf", width=7, height=6)
par(mfrow=c(2, 2),mar=c(1,1 , 1, 1), oma=c(3, 4, 1, 1))
plot(fitFlow.ag$all0.1, fitFlow.ag$t10.G1.1, xaxt="n", yaxt="n", pch=19, col=fitFlow.ag$col72, ann=F, xlim=c(0, 2), cex=1.5, ylim=c(140, 260))
axis(1, labels=FALSE)
axis(2, las=2)
mtext("Assessed at 24h", side=3, adj=0.01)
abline(lm(fitFlow.ag$t10.G1.1~fitFlow.ag$all0.1))
text(-0.1, 145, "cor = -0.45", pos=4, cex=1.25)
txt <- "Mean evolved genome size"
mtext(txt, side=2, outer=FALSE, line=2.5)
#cor.test(fitFlow.ag$t10.G1.1, fitFlow.ag$all0.1)
#t = -2.1448, df = 18, p-value = 0.04587

par(mar=c(1,1 , 1, 1))
plot(fitFlow.ag$all0.1.72, fitFlow.ag$t10.G1.1, xaxt="n", yaxt="n", pch=19, col=fitFlow.ag$col72, ann=F, xlim=c(0, 2), cex=1.5, ylim=c(140, 260))
axis(1, labels=FALSE)
axis(2, labels=FALSE)
mtext("Assessed at 72h", side=3, adj=0.01)
abline(lm(fitFlow.ag$t10.G1.1~fitFlow.ag$all0.1.72), lty = 2)
#cor.test(fitFlow.ag$t10.G1.1, fitFlow.ag$all0.1.72)
#t = -0.6966, df = 18, p-value = 0.495

par(mar=c(1,1 , 1, 1))
plot(fitFlow.ag$all0.1, fitFlow.sd$t10.G1.1, xaxt="n", yaxt="n", pch=19, col=fitFlow.ag$col, ann=F, xlim=c(0, 2), ylim=c(0, 90), cex=1.5)
axis(1, labels=FALSE)
axis(2, las=2)
abline(lm(fitFlow.sd$t10.G1.1~fitFlow.ag$all0.1))
#cor.test(fitFlow.sd$t10.G1.1, fitFlow.ag$all0.1)
#t = -3.2213, df = 18, p-value = 0.004735
axis(1)
txt <- "Evolved variation in genome size (sd)"
mtext(txt, side=2, outer=FALSE, line=2.5)
text(-0.1, 5, "cor = -0.60", pos=4, cex=1.25)

plot(fitFlow.ag$all0.1.72, fitFlow.sd$t10.G1.1, xaxt="n", yaxt="n", pch=19, col=fitFlow.ag$col, ann=F, xlim=c(0, 2), ylim=c(0, 90), cex=1.5)
axis(1, labels=TRUE)
axis(2, las=2, labels=FALSE)
abline(lm(fitFlow.sd$t10.G1.1~fitFlow.ag$all0.1.72), lty=2)
#cor.test(fitFlow.sd$t10.G1.1, fitFlow.ag$all0.1.72)
mtext("Initial growth ability", side=1, outer=TRUE, line=1.5)
dev.off()
system("open manuscript/figures/FigureSx-ancOD1cor-Ploidy_mean_SD.pdf")

pdf("manuscript/figures/FigureSx-SD-OD1_cor_SD-G1_SD-SMG72.pdf", width=7, height=6)
par(mfrow=c(2, 2),mar=c(1,1 , 1, 1), oma=c(4, 4, 1, 1))
par(mar=c(1,1 , 1, 1))
plot(fitFlow.sd$all10.1, fitFlow.sd$t10.G1.1, xaxt="n", yaxt="n", pch=19, col=fitFlow.ag$col, ann=F, xlim=c(0, 0.5), ylim=c(0, 90), cex=1.5)
axis(1, labels=FALSE)
axis(2, las=2)
abline(lm(fitFlow.sd$t10.G1.1~fitFlow.sd$all10.1))
#cor.test(fitFlow.sd$t10.G1.1, fitFlow.sd$all10.1)
#t = 5.6545, df = 18, p-value = 0.00002304
#txt <- "Evolved variation in genome size (sd)"
mtext(expression("Evolved variation\nin genome size (sd)"), side=2, outer=FALSE, line=2.5)
text(0.3, 5, "cor = 0.80", pos=4, cex=1.25)
mtext("Assessed at 24h", side=3, adj=0.01)

plot(fitFlow.sd$all10.1.72, fitFlow.sd$t10.G1.1, xaxt="n", yaxt="n", pch=19, col=fitFlow.ag$col, ann=F, xlim=c(0, 0.5), ylim=c(0, 90))
axis(1, labels=FALSE)
axis(2, las=2, labels=FALSE)
#cor.test(fitFlow.sd$t10.G1.1, fitFlow.sd$all10.1.72)
abline(lm(fitFlow.sd$t10.G1.1~fitFlow.sd$all0.1.72), lty=2)
mtext("Assessed at 72h", side=3, adj=0.01)

plot(fitFlow.sd$all10.1, fitFlow.sd$SMG72.10, xaxt="n", yaxt="n", pch=19, col=fitFlow.ag$col, ann=F, xlim=c(0, 0.5), ylim=c(0, 0.5), cex=1.5)
axis(1, labels=FALSE)
axis(2, las=2)
abline(lm(fitFlow.sd$SMG72.10~fitFlow.sd$all10.1))
#cor.test(fitFlow.sd$SMG72.10, fitFlow.sd$all10.1)
#t = 3.167, df = 18, p-value = 0.005336
axis(2, las=2)
axis(1)
#txt <- "Evolved variation\nin tolerance (SMG72)"
mtext(expression("Evolved variation\nin tolerance (SMG72)"), side=2, outer=FALSE, line=2.5)
text(0.3, 0.05, "cor = 0.60", pos=4, cex=1.25)

plot(fitFlow.sd$all10.1.72, fitFlow.sd$SMG72.10, xaxt="n", yaxt="n", pch=19, col=fitFlow.ag$col, ann=F, xlim=c(0, 0.5), ylim=c(0, 0.5), cex=1.5)
axis(1, labels=TRUE)
axis(2, las=2, labels=FALSE)
#cor.test(fitFlow.sd$SMG72.10, fitFlow.sd$all10.1.72)
#t = 1.6282, df = 18, p-value = 0.1209
#abline(lm(fitFlow.sd$SMG72.10~fitFlow.sd$all0.1.72), lty=2)
mtext("Evolved variance in growth (FLC1)", side=1, outer=TRUE, line=2.5)
dev.off()

pdf("manuscript/figures/FigureSx-ancOD1_cor_SMG-mean-SD.pdf", width=7, height=6)
par(mfrow=c(2, 2),mar=c(1,1 , 1, 1), oma=c(3, 4, 1, 1))
plot(fitFlow.ag$all0.1, fitFlow.ag$change72.SMG.2, xaxt="n", yaxt="n", pch=19, col=fitFlow.ag$col72, ann=F, xlim=c(0, 2), cex=1.5, ylim=c(-0.3, 0.6))
axis(1, labels=FALSE)
axis(2, las=2)
mtext("Assessed at 24h", side=3, adj=0.01)
abline(lm(fitFlow.ag$change72.SMG.2~fitFlow.ag$all0.1), lty=2)
#text(-0.1, 145, "cor = -0.45", pos=4, cex=1.25)
txt <- expression(paste(Delta," tolerance (SMG72)"))
mtext(txt, side=2, outer=FALSE, line=2.5)
#cor.test(fitFlow.ag$change72.SMG.2, fitFlow.ag$all0.1)
#t = -0.38601, df = 18, p-value = 0.704

par(mar=c(1,1 , 1, 1))
plot(fitFlow.ag$all0.1.72, fitFlow.ag$change72.SMG.2, xaxt="n", yaxt="n", pch=19, col=fitFlow.ag$col72, ann=F, xlim=c(0, 2), cex=1.5, ylim=c(-0.3, 0.6))
axis(1, labels=FALSE)
axis(2, labels=FALSE)
mtext("Assessed at 72h", side=3, adj=0.01)
abline(lm(fitFlow.ag$change72.SMG.2~fitFlow.ag$all0.1.72), lty = 2)
#cor.test(fitFlow.ag$t10.G1.1, fitFlow.ag$all0.1.72)
#t = -0.6966, df = 18, p-value = 0.495

par(mar=c(1,1 , 1, 1))
plot(fitFlow.ag$all0.1, fitFlow.sd$SMG72.10, xaxt="n", yaxt="n", pch=19, col=fitFlow.ag$col, ann=F, xlim=c(0, 2),  ylim=c(0, 0.6), cex=1.5)
axis(1, labels=FALSE)
axis(2, las=2)
abline(lm(fitFlow.sd$SMG72.10~fitFlow.ag$all0.1))
#cor.test(fitFlow.sd$SMG72.10, fitFlow.ag$all0.1)
#t = -2.9812, df = 18, p-value = 0.008006
axis(2, las=2)
axis(1)
txt <- "Evolved variation in tolerance (SMG72)"
mtext(txt, side=2, outer=FALSE, line=2.5)
text(-0.1, 0.05, "cor = -0.57", pos=4, cex=1.25)

par(mar=c(1,1 , 1, 1))
plot(fitFlow.ag$all0.1.72, fitFlow.sd$SMG72.10, xaxt="n", yaxt="n", pch=19, col=fitFlow.ag$col, ann=F, xlim=c(0, 2),  ylim=c(0, 0.6), cex=1.5)
axis(1, labels=TRUE)
axis(2, labels=FALSE)
abline(lm(fitFlow.sd$SMG72.10~fitFlow.ag$all0.1.72), lty=2)
#cor.test(fitFlow.sd$SMG72.10, fitFlow.ag$all0.1.72)
#t = -1.8053, df = 18, p-value = 0.08779
mtext("Initial growth ability", side=1, outer=TRUE, line=1.5)
dev.off()
system("open manuscript/figures/FigureSx-ancOD1_cor_SMG-mean-SD.pdf")
