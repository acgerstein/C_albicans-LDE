#################################
#Load libraries and functions
#################################
library(Hmisc)
library(lme4)
library(car)
library(agricolae)
library(lmerTest)
library(ggbeeswarm)
library(beeswarm)
library(RColorBrewer)
library(nlme)
library(MASS)
library(dplyr)
options(stringsAsFactors = FALSE, max.print=5000)

cv <-function(x) sd(x,na.rm=TRUE)/mean(x,na.rm=TRUE)
se <- function(x, na.rm=TRUE) sqrt(var(x, na.rm=TRUE)/(length(x) - 1))
source("data_in/general/mvee.R")

#################################
#Load data
#################################
source("scripts/04_combine-MIC-flow.R")

sub <- data.frame(place = as.factor(all10.1$place), data= all10.1$all10.1, col =all10.1$col)

summary(lmer(all10.1.72 ~ log2(MIC24.10) * SMG72.10.2 + (1|line), data = fitFlow))

fitFlow_1 <- subset(fitFlow, MIC24.10 == 1)

print(summary(lmer(all10.1.72 ~ SMG72.10.2 + (1|line), data = fitFlow_1)))

#################################
#Fitness in FLC1
#################################
#Figure 1
pdf("manuscript/figures/Figure1-OD1.pdf", width=7, height=5.5)
par(mfrow=c(2, 1),mar=c(1,1 , 1, 1), oma=c(3, 3, 1, 1))
#plot(all10.1$place, all10.1$data, col=all10.1$col, ylim=c(0, 2), xaxt="n", yaxt="n", ann=F, cex=0.7)
beeswarm(all10.1~place72, data = all10.1, corral = "wrap", col=coloursVa, ylim=c(0, 2), xaxt="n", yaxt="n", ann=F, cex=0.7, pch=19)
points(1:20, all0.ag.1$data[place72], pch="-", cex=2.25, col=grey(0.3))
points(1:20, all10.ag.1$data[place72], pch="-", cex=2.25, col=coloursV)
axis(1, 1:20, labels=FALSE)
axis(2, las=2)
mtext("OD at 24h", side=3, adj=0.01)

#plot(all10.1.72$place, all10.1.72$data, col=all10.1.72$col, ylim=c(0, 2), xaxt="n", yaxt="n", ann=F, cex=0.7)
beeswarm(all10.1.72~place72, data = all10.1.72, corral = "wrap", col=coloursVa, ylim=c(0, 2), xaxt="n", yaxt="n", ann=F, cex=0.7, pch=19)
points(1:20, all0.ag.1.72$data[place72], pch="-", cex=2.25, col=grey(0.4))
points(1:20, all10.ag.1.72$data[place72], pch="-", cex=2.25, col=coloursV)
axis(1, seq(1, 20, 2), paste0("A", all0.ag.1$line[place72][seq(1, 20, 2)]), cex.axis=0.8)
axis(1, seq(2, 20, 2), paste0("A", all0.ag.1$line[place72][seq(2, 20, 2)]), cex.axis=0.8)
axis(2, las=2)
mtext("strain", side=1, line=2)
mtext("Growth ability in FLC1 (optical density)", side=2, line=1.5, outer=TRUE)
mtext("OD at 72h", side=3, adj=0.01)

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
pdf("manuscript/figures/Figure2-OD1cor-mean_SD.pdf", width=7, height=6)
par(mfrow=c(2, 2),mar=c(1,1 , 1, 1), oma=c(3, 4, 1, 1))
plot(fitFlow.ag$all0.1, fitFlow.ag$change24ave, xaxt="n", yaxt="n", pch=19, col=fitFlow.ag$col72, ann=F, xlim=c(0, 2), ylim=c(-0.2, 1), cex=1.5)
axis(1, labels=FALSE)
axis(2, las=2)
mtext("Assessed at 24h", side=3, adj=0.01)
abline(lm(fitFlow.ag$change24ave~fitFlow.ag$all0.1))
text(-0.1, -0.18, "cor = -0.58", pos=4, cex=1.25)
txt <- expression(paste(Delta," mean growth improvement"))
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
mtext("Initial growth ability", side=1, outer=TRUE, line=1.5)


# xl <- 1
# yb <- 1
# xr <- 1.5
# yt <- 2
#
# par(mar=c(8,12.5,0.5,0.5))
# par(new=TRUE)
# plot(NA,type="n",ann=FALSE,xlim=c(1,2),ylim=c(1,2),xaxt="n",yaxt="n",bty="n")
# rect(xl, head(seq(yb,yt,(yt-yb)/20),-1), xr, tail(seq(yb,yt,(yt-yb)/20),-1), col=coloursV, border=NA)
#
# mtext(expression("low anc.\n fitness"),side=1,cex=0.6, adj=-1, line=0.25)
# mtext(expression("high anc.\n fitness"),side=3,cex=0.6, adj=-2, line=-0.25)

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
#plot(all10.0$place, all10.0$data, col=all10.0$col, ylim=c(0, 2), xaxt="n", yaxt="n", ann=F, cex=0.9)
points(1:20, all0.ag.0$data[place72], col=grey(0.3), pch="-", cex=2.25)
points(1:20, all10.ag.0$data[place72], col=coloursV, pch="-", cex=2.25)
axis(1, 1:20, labels=FALSE)
axis(2, las=2)
mtext("Growth ability at 24 h", side=3, line=2, outer=FALSE)

beeswarm(all10.0.72~place72, data = all10.0.72, corral = "wrap", col=coloursVa, ylim=c(0, 2), xaxt="n", yaxt="n", ann=F, cex=0.7, pch=19)
#plot(all10.0.72$place, all10.0.72$data, ylim=c(0, 2), xaxt="n", yaxt="n", col=all10.0.72$col, ann=F, cex=0.9)
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

fitFlow.tab <- data.frame(table(as.numeric(fitFlow.plot$line), fitFlow.plot$tabMIC10))
names(fitFlow.tab) <- c("line", "MIC24.10", "Freq")

fitFlow.tab.wide <- as_tibble(fitFlow.tab) %>%
  spread(MIC24.10, Freq)
fitFlow.tab.wide.flip <- t(fitFlow.tab.wide)
fitFlow.tab.wide.flip[2,] <- as.numeric(fitFlow.tab.wide.flip[2,])

pdf("manuscript/figures/Figure4-BMD_resist-bar.pdf", width=8.5, height=4.5)
par(xpd=NA, oma=c(1, 1, 1, 8))
m <- barplot(as.matrix(fitFlow.tab.wide.flip[2:4,as.numeric(order72_MIC)]), beside=FALSE, xaxt="n", yaxt="n", col=c("darkgrey", "purple", "orange"),  border=NA, space=c(rep(0.2, 11), 1.5, rep(0.2, 3), 1.5, rep(0.2, 5)), xlim=c(0.7, 26.1))
axis(1, m[1:20][seq(1, 20, 2)], labels=FALSE)
text(m[1:20][seq(1, 20, 2)], -1, paste0("A", order72_MIC[seq(1, 20, 2)]), cex=0.6)
axis(1, m[1:20][seq(2, 21, 2)], labels = FALSE)
text(m[1:20][seq(2, 21, 2)], -1, paste0("A", order72_MIC[seq(2, 21, 2)]), cex=0.6)
axis(2, las=2)
box()
legend(29, 12, legend=c(expression(MIC[50] ~ "< 1", MIC[50] ~ "= 1", MIC[50] ~ "> 1" )), inset=0.05, pch=22, col=c("darkgrey", "purple", "orange"), pt.bg =c("darkgrey", "purple", "orange"), cex=0.7)
mtext("number of evolved replicates", side=2, line=2)
mtext(expression("Ancestral \n" ~ MIC[50] ~ "< 1"), side= 1, adj = 0.2, line = 3, cex= 0.9)
mtext(expression("Ancestral \n" ~ MIC[50] ~ "= 1"), side= 1, adj = 0.65, line = 3, cex= 0.9)
mtext(expression("Ancestral \n" ~ MIC[50] ~ "\n> 1"), side= 1, adj = 1, line=3, cex= 0.9)
dev.off()



pdf("manuscript/figures/Figure4-BMD_resist.pdf", width=7.5, height=4.5)
plot(1:11, log2(MIC24.ag.less1$MIC24), ylim=c(log2(0.25), log2(256)), yaxt="n", xaxt="n", xlab="Strain", ylab="", pch="-", cex=2.25, col=grey(0.3), xlim=c(1, 20))
points(12:15, log2(MIC24.ag.1$MIC24), pch="-", cex=2.25, col=grey(0.3))
points(16:20, log2(MIC24.ag.more1$MIC24), pch="-", cex=2.25, col=grey(0.3))
points(jitter(rep(1:11, each=12)), log2(MIC24.less1$MIC24.10), pch=21, col=MIC24.less1$col72)
abline(v=11.5, lty=2)
points(jitter(rep(12:15, each=12)), log2(MIC24.1$MIC24.10), pch=21, col=MIC24.1$col72)
abline(v=15.5, lty=2)
points(jitter(rep(16:20, each=12)), log2(MIC24.more1$MIC24.10), pch=21, col=MIC24.more1$col72)
axis(2, at=c(log2(0.25),log2(1), log2(4), log2(16), log2(64), log2(256)), labels=c("0.25", "1", "4", "16", "64", ">128"), las=2)
#abline(h=log(8), lty=2)
mtext(expression(MIC[50]), side=2, line=2)
axis(1, 1:20, labels=c(paste0("A", c(MIC24.ag.less1$line, MIC24.ag.1$line, MIC24.ag.more1$line))), cex.axis=0.8)
axis(1, c(2, 6,8, 10, 12, 15, 18, 20) , labels=c(paste0("A", c(MIC24.ag.less1$line, MIC24.ag.1$line, MIC24.ag.more1$line)[c(2, 6,8, 10, 12, 15, 18, 20)])),  cex.axis=0.8)
#text(1:20, -3, paste0("A", fitFlow.plot.ag$line[place]), srt=-45, adj=0.1, xpd=NA)
dev.off()

# pdf("manuscript/figures/Figure4-BMD_resist.pdf", width=7, height=4.5)
# plot(1:20, log(fitFlow.plot.ag$MIC24)[place72], ylim=c(log(0.125), log(256)), yaxt="n", xaxt="n", xlab="Strain", ylab="", pch="-", cex=2.25, col=grey(0.3))
# points(jitter(as.numeric(fitFlow.plot$place)), log(fitFlow.plot$MIC24.10), pch=21, col=fitFlow.plot$col)
# axis(2, at=c(log(0.25),log(1), log(4), log(16), log(64), log(256)), labels=c("<1", "1", "4", "16", "64", ">128"), las=2)
# axis(1, 1:20, labels=FALSE, cex.axis=0.5)
# abline(h=log(8), lty=2)
# mtext(expression(MIC[50]), side=2, line=2)
# axis(1, 1:20, labels=FALSE)
# axis(1, seq(1, 20, 2), paste0("A", all0.ag.1$line[place72][seq(1, 20, 2)]), cex.axis=0.8)
# axis(1, seq(2, 20, 2), paste0("A", all0.ag.1$line[place72][seq(2, 20, 2)]), cex.axis=0.8)
# #text(1:20, -3, paste0("A", fitFlow.plot.ag$line[place]), srt=-45, adj=0.1, xpd=NA)
# dev.off()
#system("open /Users/acgerstein/Documents/Postdoc/Papers/MutAccum/Figures/Figure4-BMD_resist.pdf")

#################################
#Tolerance above 1ug
#################################

pdf("manuscript/figures/FigureS3-tolerance-above1ug.pdf", width=7, height=6.5)
fitFlow$place <- as.factor(fitFlow$place)
par(mfrow=c(3, 1), mar=c(1, 1, 1, 1), oma=c(3, 4, 1, 1), mgp=c(1.5, 0.75, 0))
beeswarm(SMG24.10.3~place72, data = fitFlow, ylim=c(0, 1.6), yaxt="n", xaxt="n", xlab="strain", ylab="", pch=19, cex=1.2, col=coloursVa, corral="wrap")
points(1:20, fitFlow.ag$SMG24.3[place72], pch="-", cex=2.25, col=grey(0.2))
points(1:20, fitFlow.ag$SMG24.up.3[place72], pch="-", cex=2.25, col=grey(0.2))
points(1:20, fitFlow.ag$SMG24.down.3[place72], pch="-", cex=2.25, col=grey(0.2))
axis(2, las=2)
axis(1, 1:20, labels=FALSE, cex.axis=0.5)
mtext("a 24 h" , side=3, adj=0.01, font=2, cex=1)

beeswarm(SMG48.10.3~place72, data = fitFlow, ylim=c(0, 1.6), yaxt="n", xaxt="n", xlab="strain", ylab="", pch=19, cex=1.2, col=coloursVa, corral="wrap")
points(1:20, fitFlow.ag$SMG48.3[place72], pch="-", cex=2.25, col=grey(0.2))
points(1:20, fitFlow.ag$SMG48.up.3[place72], pch="-", cex=2.25, col=grey(0.2))
points(1:20, fitFlow.ag$SMG48.down.3[place72], pch="-", cex=2.25, col=grey(0.2))
axis(2, las=2)
axis(1, 1:20, labels=FALSE, cex.axis=0.5)
mtext("b 48 h" , side=3, adj=0.01, font=2, cex=1)

beeswarm(SMG72.10.3~place72, data = fitFlow, ylim=c(0, 1.6), yaxt="n", xaxt="n", xlab="strain", ylab="", pch=19, cex=1.2, col=coloursVa, corral="wrap")
points(1:20, fitFlow.ag$SMG72.3[place72], pch="-", cex=2.25, col=grey(0.2))
points(1:20, fitFlow.ag$SMG72.up.3[place72], pch="-", cex=2.25, col=grey(0.2))
points(1:20, fitFlow.ag$SMG72.down.3[place72], pch="-", cex=2.25, col=grey(0.2))

axis(2, las=2)
axis(1, 1:20, labels=FALSE, cex.axis=0.5)
mtext("c 72 h" , side=3, adj=0.01, font=2, cex=1)
mtext("Evolved tolerance" , side=2, outer=TRUE, line=2)
axis(1, seq(1, 20, 2), paste0("A", all0.ag.1$line[place72][seq(1, 20, 2)]), cex.axis=0.9)
axis(1, seq(2, 20, 2), paste0("A", all0.ag.1$line[place72][seq(2, 20, 2)]), cex.axis=0.9)
#text(1:20, -0.15, paste0("A", fitFlow.plot.ag$line[place]), srt=-45, adj=0.1, xpd=NA)
dev.off()

order <- c(13, 17, 11, 3, 15, 9, 14, 4, 16, 19, 10, 8, 20, 7, 2, 6, 12, 5, 1, 18)

#ordered by ancestral SMG - Figure 5
#SMGorder <- order(fitFlow.ag$all0.1.72)
pdf("manuscript/figures/Figure5-SMG72-above1ug-flip.pdf", width=8, height=5)
plot(rep(1:11, each =12), MIC24.less1$SMG72.10.3, yaxt="n", xaxt="n", xlab="Strain", ylab="", cex=1.2, pch=19, xlim=c(1, 20), col =MIC24.less1$col72a, ylim=c(0,1.7))
points(rep(12:15, each=12), MIC24.1$SMG72.10.3, cex=1.2, col=MIC24.1$col72a, pch=19)
points(rep(16:20, each=12), MIC24.more1$SMG72.10.3, cex=1.2, col=MIC24.more1$col72a, pch=19)
abline(v=11.5, lty=2)
abline(v=15.5, lty=2)
axis(1, 1:20, labels=c(paste0("A", c(MIC24.ag.less1$line, MIC24.ag.1$line, MIC24.ag.more1$line))), cex.axis=0.8)
axis(1, c(2, 7, 10, 19) , labels=c(paste0("A", c(MIC24.ag.less1$line, MIC24.ag.1$line, MIC24.ag.more1$line)[c(2, 7, 10, 19)])),  cex.axis=0.8)

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

pdf("manuscript/figures/Figure5-SMG72-aboveMIC-flip_MIC24-1.pdf", width=8, height=5)
plot(rep(1:11, c(11, 3, 5, 5, 12, 10, 12, 9, 10, 8, 12)), MIC24.less1_1$SMG72.10.2, yaxt="n", xaxt="n", xlab="Strain", ylab="", cex=1.2, pch=19, xlim=c(1, 20), col = MIC24.less1_1$col72a, ylim=c(0,1))
points(rep(12:15, c(11, 12, 11, 8)), MIC24.1_1$SMG72.10.2, cex=1.2, col=MIC24.1_1$col72a, pch=19)
points(rep(16:20, c(0, 2, 0, 12, 6)), MIC24.more1_1$SMG72.10.2, cex=1.2, col=MIC24.more1_1$col72a, pch=19)
abline(v=11.5, lty=2)
abline(v=15.5, lty=2)
axis(1, 1:20, labels=c(paste0("A", c(unique(MIC24.less1_1$line), unique(MIC24.1_1$line), unique(MIC24.more1$line)))), cex.axis=0.8)
axis(1, c(2, 7, 10, 19) , labels=c(paste0("A", c(MIC24.ag.less1$line, MIC24.ag.1$line, MIC24.ag.more1$line)[c(2, 7, 10, 19)])),  cex.axis=0.8)

points(1:11, MIC24.ag.less1$SMG72.up, pch="-", col="black", cex=2)
points(1:11, MIC24.ag.less1$SMG72.down, pch="-", col="black", cex=2)
arrows(1:11, MIC24.ag.less1$SMG72.up, 1:11, MIC24.ag.less1$SMG72.down, length=0, col="black")
points(12:15, MIC24.ag.1$SMG72.up, pch="-", col="black", cex=2)
points(12:15, MIC24.ag.1$SMG72.down, pch="-", col="black", cex=2)
arrows(12:15, MIC24.ag.1$SMG72.up, 12:15, MIC24.ag.1$SMG72.down, length=0, col="black")
points(16:20, MIC24.ag.more1$SMG72.up, pch="-", col="black", cex=2)
points(16:20, MIC24.ag.more1$SMG72.down, pch="-", col="black", cex=2)
arrows(16:20, MIC24.ag.more1$SMG72.up, 16:20, MIC24.ag.more1$SMG72.down, length=0, col="black")


# points(16:18, MIC24.ag.more1$SMG72.up[c(2, 4, 5)], pch="-", col="black", cex=2)
# points(16:18, MIC24.ag.more1$SMG72.down[c(2, 4, 5)], pch="-", col="black", cex=2)
# arrows(16:18, MIC24.ag.more1$SMG72.up[c(2, 4, 5)], 16:18, MIC24.ag.more1$SMG72.down[c(2, 4, 5)], length=0, col="black")
axis(2, las=2)
mtext("Tolerance above MIC (72h)", side=2, line=3)
dev.off()



pdf("manuscript/figures/Figure5-SMG72-aboveMIC-flip.pdf", width=8, height=5)
plot(rep(1:11, each =12), MIC24.less1$SMG72.10.2, yaxt="n", xaxt="n", xlab="Strain", ylab="", cex=1.2, pch=19, xlim=c(1, 20), col =MIC24.less1$col72a, ylim=c(0,1))
points(rep(12:15, each=12), MIC24.1$SMG72.10.2, cex=1.2, col=MIC24.1$col72a, pch=19)
points(rep(16:20, each=12), MIC24.more1$SMG72.10.2, cex=1.2, col=MIC24.more1$col72a, pch=19)
abline(v=11.5, lty=2)
abline(v=15.5, lty=2)
axis(1, 1:20, labels=c(paste0("A", c(MIC24.ag.less1$line, MIC24.ag.1$line, MIC24.ag.more1$line))), cex.axis=0.8)
axis(1, c(2, 7, 10, 19) , labels=c(paste0("A", c(MIC24.ag.less1$line, MIC24.ag.1$line, MIC24.ag.more1$line)[c(2, 7, 10, 19)])),  cex.axis=0.8)
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
#Ploidy (non)correlations
#################################

pdf("manuscript/figures/Figure7-evolFit-evolPloidy-24h-SMG72.pdf", width=7, height=5)
k<-0
par(mar=c(1, 1, 1, 1), mfrow=c(4, 4), mar=c(1, 1, 1, 1), oma=c(4, 4, 1, 1))
for(i in c(place[1:16])){
  k <- k+1
  plot(subset(fitFlow, line==i)$all10.1, subset(fitFlow, line==i)$t10.G1.1, ylim=c(100, 500), xlim=c(0, 2), xaxt="n", yaxt="n", pch=subset(fitFlow, line==i)$pch, bg=ifelse(subset(fitFlow, line ==i)$change72.SMG.3 >0, "black", "white"))
  t <- cor.test(subset(fitFlow, line==i)$all10.1, subset(fitFlow, line==i)$t10.G1.1, method="spearman")
  if(t$p.value < 0.05) abline(lm(subset(fitFlow, line==i)$t10.G1.1~subset(fitFlow, line==i)$all10.1), col="red")
  if(t$p.value > 0.05) abline(lm(subset(fitFlow, line==i)$t10.G1.1~subset(fitFlow, line==i)$all10.1), col="red", lty=2)
  if(k > 12) axis(1)
  else axis(1, labels=FALSE)
  if(k %% 4 == 1) axis(2, las=2)
  else axis(2, labels=FALSE)
  text(0, 450, paste0("A", i), pos=4, font=2, cex=1.1)
}
mtext("Evolved growth  in low drug", side=1, outer=TRUE, line=2)
mtext("Evolved genome size (FITC intensity)", side=2, outer=TRUE, line=2)
dev.off()


pdf("manuscript/figures/FigureS4-evolFit-evolPloidy-72h-SMG72.pdf", width=7, height=5)
k<-0
par(mar=c(1, 1, 1, 1), mfrow=c(4, 4), mar=c(1, 1, 1, 1), oma=c(4, 4, 1, 1))
for(i in c(place[1:16])){
  k <- k+1
  plot(subset(fitFlow, line==i)$all10.1.72, subset(fitFlow, line==i)$t10.G1.1, ylim=c(100, 500), xlim=c(0, 2), xaxt="n", yaxt="n", pch=subset(fitFlow, line==i)$pch, bg=ifelse(subset(fitFlow, line ==i)$change72.SMG.3 >0, "black", "white"))
  t <- cor.test(subset(fitFlow, line==i)$all10.1.72, subset(fitFlow, line==i)$t10.G1.1, method="spearman")
  if(t$p.value < 0.05) abline(lm(subset(fitFlow, line==i)$t10.G1.1~subset(fitFlow, line==i)$all10.1.72), col="red")
  if(t$p.value > 0.05) abline(lm(subset(fitFlow, line==i)$t10.G1.1~subset(fitFlow, line==i)$all10.1.72), col="red", lty=2)
  if(k > 12) axis(1)
  else axis(1, labels=FALSE)
  if(k %% 4 == 1) axis(2, las=2)
  else axis(2, labels=FALSE)
  text(0, 450, paste0("A", i), pos=4, font=2, cex=1.1)
}
mtext("Evolved growth in drug (72 h OD)", side=1, outer=TRUE, line=2)
mtext("Evolved genome size (FITC intensity)", side=2, outer=TRUE, line=2)
dev.off()

#THESE SEEM UNNECESSARY
pdf("/Users/acgerstein/Documents/Postdoc/Papers/MutAccum/Figures/FigureS5-SMG24-evolPloidy.pdf", width=7, height=5)
k<-0
par(mar=c(1, 1, 1, 1), mfrow=c(4, 4), mar=c(1, 1, 1, 1), oma=c(4, 4, 1, 1))
for(i in c(place[1:16])){
  k <- k+1
  plot(subset(fitFlow.variable, line==i)$SMG24.10, subset(fitFlow.variable, line==i)$t10.G1.1, ylim=c(100, 500), xlim=c(0, 1), xaxt="n", yaxt="n")
  t <- cor.test(subset(fitFlow.variable, line==i)$SMG24.10, subset(fitFlow.variable, line==i)$t10.G1.1, method="spearman")
  if(t$p.value < 0.05) abline(lm(subset(fitFlow.variable, line==i)$t10.G1.1~subset(fitFlow.variable, line==i)$SMG24.10), col="red")
  if(t$p.value > 0.05) abline(lm(subset(fitFlow.variable, line==i)$t10.G1.1~subset(fitFlow.variable, line==i)$SMG24.10), col="red", lty=2)
  if(k > 12) axis(1)
  else axis(1, labels=FALSE)
  if(k %% 4 == 1) axis(2, las=2)
  else axis(2, labels=FALSE)
  text(0, 450, paste0("A", i), pos=4, font=2, cex=1.1)
}
mtext("Evolved tolerance (SMG24) ", side=1, outer=TRUE, line=2)
mtext("Evolved genome size (FITC intensity)", side=2, outer=TRUE, line=2)
dev.off()


#Evolved ploidy * evolved SMG72
pdf("/Users/acgerstein/Documents/Postdoc/Papers/MutAccum/Figures/FigureS6-SMG72-evolPloidy.pdf", width=7, height=5)
k<-0
par(mar=c(1, 1, 1, 1), mfrow=c(4, 4), mar=c(1, 1, 1, 1), oma=c(4, 4, 1, 1))
for(i in c(place[1:16])){
  k <- k+1
  plot(subset(fitFlow, line==i)$SMG72.10.2, subset(fitFlow, line==i)$t10.G1.1, ylim=c(100, 500), xlim=c(0, 2), xaxt="n", yaxt="n", pch=subset(fitFlow, line==i)$pch)
  t <- cor.test(subset(fitFlow, line==i)$SMG72.10.2, subset(fitFlow, line==i)$t10.G1.1, method="spearman")
  if(t$p.value < 0.05) abline(lm(subset(fitFlow, line==i)$t10.G1.1~subset(fitFlow, line==i)$SMG72.10.2), col="red")
  if(t$p.value > 0.05) abline(lm(subset(fitFlow, line==i)$t10.G1.1~subset(fitFlow, line==i)$SMG72.10.2), col="red", lty=2)
  if(k > 12) axis(1)
  else axis(1, labels=FALSE)
  if(k %% 4 == 1) axis(2, las=2)
  else axis(2, labels=FALSE)
  text(2, 450, paste0("A", i), pos=2, font=2, cex=1.1)
}
mtext("Evolved tolerance (SMG72)", side=1, outer=TRUE, line=2)
mtext("Evolved genome size (FITC intensity)", side=2, outer=TRUE, line=2)
dev.off()
