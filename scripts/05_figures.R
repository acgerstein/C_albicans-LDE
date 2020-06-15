#################################
#Load data and libraries
#################################
source("scripts/04_combine-MIC-flow.R")


#################################
#Fitness in FLC1
#################################
#Figure 1
pdf("manuscript/figures/Figure1-OD1.pdf", width=7, height=5.5, font = "Times")
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
mtext(expression("Growth in YPD + 1"~mu~"g/mL FLC (optical density)"), side=2, line=1.5, outer=TRUE)
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

mtext(expression("low parental\n fitness"),side=1,cex=0.6, adj=-1, line=0.25)
mtext(expression("high parental\n fitness"),side=3,cex=0.6, adj=-10, line=-0.25)

dev.off()
system("open manuscript/figures/Figure1-OD1.pdf")

#Figure 2
fitFlow.ag$change24ave <- fitFlow.ag$all10.1 - fitFlow.ag$all0.1
fitFlow.ag$change72ave <- fitFlow.ag$all10.1.72 - fitFlow.ag$all0.1.72

pdf("manuscript/figures/Figure2-OD1cor-mean_SD.pdf", width=7, height=6)
par(mfrow=c(2, 2),mar=c(1,1 , 1, 1), oma=c(3, 4, 1, 1))
plot(fitFlow.ag$all0.1, c(fitFlow.ag$all10.1-fitFlow.ag$all0.1), xaxt="n", yaxt="n", pch=19, col=fitFlow.ag$col72, ann=F, xlim=c(0, 2), ylim=c(-0.2, 1), cex=1.5)
axis(1, labels=FALSE)
axis(2, las=2)
mtext("Assessed at 24h", side=3, adj=0.01)
abline(lm(c(fitFlow.ag$all10.1-fitFlow.ag$all0.1)~fitFlow.ag$all0.1))
txt <- expression(paste(Delta," mean fitness"))
mtext(txt, side=2, outer=FALSE, line=2.5)

par(mar=c(1,1 , 1, 1))
plot(fitFlow.ag$all0.1.72, c(fitFlow.ag$all10.1.72-fitFlow.ag$all0.1.72), xaxt="n", yaxt="n", pch=19, col=fitFlow.ag$col, ann=F, xlim=c(0, 2), ylim=c(-0.2, 1), cex=1.5)
axis(1, labels=FALSE)
axis(2, labels=FALSE)
mtext("Assessed at 72h", side=3, adj=0.01)
abline(lm(c(fitFlow.ag$all10.1.72-fitFlow.ag$all0.1.72)~fitFlow.ag$all0.1.72))

plot(fitFlow.ag$all0.1, fitFlow.sd$all10.1-fitFlow.sd$all0.1, xaxt="n", yaxt="n", pch=19, col=fitFlow.ag$col, ann=F, xlim=c(0, 2), ylim=c(-0.1, 0.45), cex=1.5)
axis(1, labels=FALSE)
axis(2, las=2)
abline(lm(c(fitFlow.sd$all10.1-fitFlow.sd$all0.1)~fitFlow.ag$all0.1))
axis(2, las=2)
axis(1)
txt <- expression(paste(Delta," replicate variation"))
mtext(txt, side=2, outer=FALSE, line=2.5)

plot(fitFlow.ag$all0.1.72, fitFlow.sd$all10.1.72-fitFlow.sd$all0.1.72, xaxt="n", yaxt="n", pch=19, col=fitFlow.ag$col, ann=F, xlim=c(0, 2), ylim=c(-0.1, 0.45), cex=1.5)
axis(1)
axis(2, labels=FALSE)
abline(lm(c(fitFlow.sd$all10.1-fitFlow.sd$all0.1)~fitFlow.ag$all0.1))
mtext("Initial fitness", side=1, outer=TRUE, line=1.5)
dev.off()
system("open manuscript/figures/Figure2-OD1cor-mean_SD.pdf")

####################
#resistance & tolerance
# Figure 3
####################

fitFlow.ag$MIC24[fitFlow.ag$MIC24==0.0625] <- 0.25
fitFlow.ag$MIC24.10[fitFlow.ag$MIC24.10==0.0625] <- 0.25
fitFlow$MIC24[fitFlow$MIC24==0.0625] <- 0.25
fitFlow$MIC24.10[fitFlow$MIC24.10==0.0625] <- 0.25
#subset(fitFlow, MIC24.10 > 64) #1 line in A5, 5 in A12
fitFlow$MIC24.10[fitFlow$MIC24.10 > 64] <- 256

MIC24.ag.less1 <- subset(fitFlow.ag, MIC24 < 1)
MIC24.ag.less1 <- MIC24.ag.less1[order(MIC24.ag.less1$MIC24, MIC24.ag.less1$place72),]
MIC24.ag.1 <- subset(fitFlow.ag, MIC24 == 1)
MIC24.ag.1 <- MIC24.ag.1[order(MIC24.ag.1$place72),]
MIC24.ag.more1 <- subset(fitFlow.ag, MIC24 > 1)
MIC24.ag.more1 <- MIC24.ag.more1[order(MIC24.ag.more1$MIC24),]

MIC24.less1 <- subset(fitFlow, line %in% MIC24.ag.less1$line) #144
MIC24.less1 <- MIC24.less1[order(MIC24.less1$MIC24, MIC24.less1$place72),]
MIC24.1 <- subset(fitFlow, line %in% MIC24.ag.1$line) #36
MIC24.1 <- MIC24.1[order(MIC24.1$place72),]
MIC24.more1 <- subset(fitFlow, line %in% MIC24.ag.more1$line) #60
MIC24.more1 <- MIC24.more1[order(MIC24.more1$MIC24),]

order72_MIC <- c(MIC24.ag.less1$line, MIC24.ag.1$line, MIC24.ag.more1$line)

fitFlow$tabMIC10[fitFlow.plot$MIC24.10 < 1] <- "<1"
fitFlow$tabMIC10[fitFlow.plot$MIC24.10 == 1] <- "=1"
fitFlow$tabMIC10[fitFlow.plot$MIC24.10 > 1] <- ">1"

highMIC <- subset(fitFlow, MIC24.10 >=8)
#table(highMIC$strain)

############################################################
#Figure 3a
############################################################
pdf("manuscript/figures/Figure3a-BMD_resist.pdf", width=7.5, height=4.5)
plot(1:12, log2(MIC24.ag.less1$MIC24), ylim=c(log2(0.25), log2(256)), yaxt="n", xaxt="n", xlab="Strain", ylab="", pch="-", cex=2.25, col=grey(0.3), xlim=c(1, 22))
points(14:16, log2(MIC24.ag.1$MIC24), pch="-", cex=2.25, col=grey(0.3))
points(18:22, log2(MIC24.ag.more1$MIC24), pch="-", cex=2.25, col=grey(0.3))
points(jitter(rep(1:12, each=12)), log2(MIC24.less1$MIC24.10), pch=19, col=MIC24.less1$col72a)
abline(v=13, lty=2)
points(jitter(rep(14:16, each=12)), log2(MIC24.1$MIC24.10), pch=19, col=MIC24.1$col72a)
abline(v=17, lty=2)
points(jitter(rep(18:22, each=12)), log2(MIC24.more1$MIC24.10), pch=19, col=MIC24.more1$col72a)
axis(2, at=c(log2(0.25),log2(1), log2(4), log2(16), log2(64), log2(256)), labels=c("<1", "1", "4", "16", "64", ">64"), las=2)
#abline(h=log(8), lty=2)
mtext(expression(MIC[50]), side=2, line=2)
axis(1, c(1:12, 14:16, 18:22), labels=c(paste0("A", c(MIC24.ag.less1$line, MIC24.ag.1$line, MIC24.ag.more1$line))), cex.axis=0.8)
axis(1, c(2, 4, 6, 8, 10, 12, 15, 19, 21) , labels=c(paste0("A", c(MIC24.ag.less1$line, MIC24.ag.1$line, MIC24.ag.more1$line)[c(2, 4, 6,8, 10, 12, 14, 17, 19)])),  cex.axis=0.8)
dev.off()
system("open manuscript/figures/Figure3a-BMD_resist.pdf")

########################
#Figure S1 - MIC traces
########################
incMIC <- subset(fitFlow, MIC24.10 > MIC24 & MIC24.10 > 1)
incMIC$lr <- paste(incMIC$line, incMIC$rep, sep="_")
#incMIC <- c("2_7", "3_3", "3_5", "3_7", "3_12", "5_3", "5_10", "5_12", "7_7", "8_5", "8_6", "8_8", "17_2", "19_8", "19_9")
MICdec <- subset(fitFlow, change24 < 0)
MICdec$lr <- paste(MICdec$line, MICdec$rep, sep  = "_")
decMIC <- MICdec$lr
all10$lr <- paste(all10$line, all10$rep, sep  = "_")
all10$col <- "grey"
all10$col[all10$lr %in% decMIC] <- "red"
all10$col[all10$lr %in% incMIC$lr] <- "blue"
all10$lwd <- 0.5
all10$lwd[all10$lr %in% incMIC$lr] <- 3
all10$lwd[all10$lr %in% decMIC] <- 3
all10 <- subset(all10, enviro < 512)

all10.72$lr <- paste(all10.72$line, all10.72$rep, sep  = "_")
all10.72$col <- "grey"
all10.72$col[all10.72$lr %in% decMIC] <- "red"
all10.72$col[all10.72$lr %in% incMIC$lr] <- "blue"
all10.72$lwd <- 0.5
all10.72$lwd[all10.72$lr %in% incMIC$lr] <- 3
all10.72$lwd[all10.72$lr %in% decMIC] <- 3
all10.72 <- subset(all10.72, enviro < 512)


pdf("manuscript/figures/FigureS1-MICTraces.pdf", width=7.5, height=6)
j <- 0
par(mfrow=c(5, 4), mar=c(1, 1, 1, 1), oma=c(3, 4, 1, 1))
for(i in order72_MIC){
  j <- j+1
  sub_t0 <- subset(all0, line == i)
  sub_t0_ag <- aggregate(sub_t0["data"], sub_t0[c("enviro")], median)
  sub <- subset(all10, line == i)
  sub$enviro[sub$enviro == 0] <- 0.0625
  sub_ag <-  aggregate(sub["dataM"], sub[c("enviro", "rep", "lr", "col", "lwd")], median)
  sub_ag1 <- subset(sub_ag, rep == 1)
  plot(log(sub_ag1$enviro), sub_ag1$dataM, type="l", yaxt="n", xaxt="n", ylim=c(0, 2), lwd=0.5, col="grey")
  abline(v=log(fitFlow.ag$MIC24[as.numeric(i)]), col="grey", lwd=2)
  for (k in 1:12){
    sub_ag_sub <- subset(sub_ag, rep == k)
    points(log(sub_ag_sub$enviro), sub_ag_sub$dataM, type="l", col= sub_ag_sub$col, lwd = sub_ag_sub$lwd)
    #points(log(sub_ag_sub$enviro), sub_ag_sub$dataM, type="l", col= ifelse(sub_ag_sub$lr[1] %in% outliers, "black", "black"), lwd = ifelse(sub_ag_sub$lr[1] %in% outliers, 3, 0.5))
  }
  points(log(sub_t0_ag$enviro), sub_t0_ag$data, type="l", lwd=3, col="black")
  if(j %% 4 ==1) axis(2, las=2)
  else axis(2, labels=FALSE)
  if (j > 16) axis(1, at = log(c(0.0625, 0.25, 1, 4, 16, 64)), labels=c("0", "0.25", "1", "4", "16", "64"))
  else axis(1, at = log(c(0.0625, 0.25, 1, 4, 16, 64)), labels=FALSE)

  #abline(h=subset(sub_t0_ag, enviro=="0.0625")$data/2, col="blue")
  mtext(paste0("A", sub$line[1]), side=3, adj=0.01)
}
mtext(expression("Fluconazole concentration ("~mu~"g)"), side=1, line=1.5, outer=TRUE)
#mtext("Fluconazole concentration (ug)", side=1, outer=TRUE, line=2)
mtext("Optical density (24 h)", side=2, outer=TRUE, line=2)
dev.off()
system("open manuscript/figures/FigureS1-MICTraces.pdf")



#################################
#Tolerance above 1ug
#################################
fitFlow_MIC1 <- subset(fitFlow, MIC24.10 == 1)
fitFlow_MIC1$clade <- as.factor(fitFlow_MIC1$clade)
fitFlow_MIC1.ag <- aggregate(fitFlow_MIC1[c("change72.SMG.2", "delta1.72", "all10.1.72", "all0.1.72")], fitFlow_MIC1[c("line", "clade", "zygosity", "col72a")], mean)
fitFlow_MIC1.sd <- aggregate(fitFlow_MIC1[c("change72.SMG.2", "delta1.72", "all10.1.72", "all0.1.72")], fitFlow_MIC1[c("line", "clade", "zygosity", "col72a")], sd)

MIC24.less1_1 <- subset(MIC24.less1, MIC24.10 == 1) #110
MIC24.1_1 <- subset(MIC24.1, MIC24.10 == 1) #36
MIC24.more1_1 <- subset(MIC24.more1, MIC24.10 == 1) #16

table(MIC24.less1_1$line)
table(MIC24.1_1$line)
table(MIC24.more1_1$line)


############################################################
# Figure S2 - Tolerance at different time points
############################################################
#place72 = coloursVa
#[1] 17 11  3  9  4 14 13  8 15  2 19 16  6  7  10 20  5 12 18  1
 #    1  2   3  4  5  6  7  8 9  10 11 12  13 14 15 16  17 18 19 20
#place72_MIC = coloursVa[c(1:5, 7, 8, 10, 12, 13, 16, 5, 8, 10)]
#[1] 17 11 3 9, 14, 13, 15, 19, 16 10, 4, 8, 2, 6, 7, 1, 5, 18, 20, 12
#      1, 2, 3,4, 6, 7, 9, 11, 12, 15, 5, 8, 10, 13, 14, 20, 17, 19, 16, 18


pdf("manuscript/figures/FigureS2-tolerance-aboveMIC.pdf", width=7, height=6.5)
par(mfrow=c(3, 1), mar=c(1, 1, 1, 1), oma=c(3, 4, 1, 1), mgp=c(1.5, 0.75, 0))
beeswarm(SMG24.10~place72_MIC, data = fitFlow, ylim=c(0, 1), yaxt="n", xaxt="n", xlab="strain", ylab="", pch=19, cex=1.2, col=coloursVa[c(1, 2, 3,4, 6, 7, 9, 11, 12, 15, 5, 8, 10, 13, 14, 20, 17, 19, 16, 18)], corral="wrap", xlim=c(1, 20))
#points(1:20, fitFlow.ag$SMG24.2[place72], pch="-", cex=2.25, col=grey(0.2))
points(1:20, SMG24.up[order72_MIC], pch="-", cex=2.25, col=grey(0.2))
points(1:20, SMG24.down[order72_MIC], pch="-", cex=2.25, col=grey(0.2))
arrows(1:20, SMG24.up[order72_MIC], 1:20, SMG24.down[order72_MIC], col=grey(0.2), length=0)
axis(2, las=2)
axis(1, 1:20, labels=FALSE, cex.axis=0.5)
mtext("a 24 h" , side=3, adj=0.01, font=2, cex=1)

beeswarm(SMG48.10~place72_MIC, data = fitFlow, ylim=c(0, 1), yaxt="n", xaxt="n", xlab="strain", ylab="", pch=19, cex=1.2, col=coloursVa[c(1, 2, 3,4, 6, 7, 9, 11, 12, 15, 5, 8, 10, 13, 14, 20, 17, 19, 16, 18)], corral="wrap")
#points(1:20, fitFlow.ag$SMG48.2[place72], pch="-", cex=2.25, col=grey(0.2))
points(1:20, SMG48.up[order72_MIC], pch="-", cex=2.25, col=grey(0.2))
points(1:20, SMG48.down[order72_MIC], pch="-", cex=2.25, col=grey(0.2))
arrows(1:20, SMG48.up[order72_MIC], 1:20, SMG48.down[order72_MIC], col=grey(0.2), length=0)
axis(2, las=2)
axis(1, 1:20, labels=FALSE, cex.axis=0.5)
mtext("b 48 h" , side=3, adj=0.01, font=2, cex=1)

beeswarm(SMG72.10~place72_MIC, data = fitFlow, ylim=c(0, 1), yaxt="n", xaxt="n", xlab="strain", ylab="", pch=19, cex=1.2, col=coloursVa[c(1, 2, 3,4, 6, 7, 9, 11, 12, 15, 5, 8, 10, 13, 14, 20, 17, 19, 16, 18)], corral="wrap")
#points(1:20, fitFlow.ag$SMG72.2[place72], pch="-", cex=2.25, col=grey(0.2))
points(1:20, SMG72.up[order72_MIC], pch="-", cex=2.25, col=grey(0.2))
points(1:20, SMG72.down[order72_MIC], pch="-", cex=2.25, col=grey(0.2))
arrows(1:20, SMG72.up[order72_MIC], 1:20, SMG72.down[order72_MIC], col=grey(0.2), length=0)

axis(2, las=2)
axis(1, 1:20, labels=FALSE, cex.axis=0.5)
mtext("c 72 h" , side=3, adj=0.01, font=2, cex=1)
mtext("Evolved tolerance" , side=2, outer=TRUE, line=2)
axis(1, seq(1, 20, 2), paste0("A", all0.ag.1$line[as.numeric(order72_MIC)][seq(1, 20, 2)]), cex.axis=0.9)
axis(1, seq(2, 20, 2), paste0("A", all0.ag.1$line[as.numeric(order72_MIC)][seq(2, 20, 2)]), cex.axis=0.9)
#text(1:20, -0.15, paste0("A", fitFlow.plot.ag$line[place]), srt=-45, adj=0.1, xpd=NA)
dev.off()
system("open manuscript/figures/FigureS2-tolerance-aboveMIC.pdf")

#################################
#Tolerance traces S3
#################################
pdf("manuscript/figures/FigureS3-MICTraces_72h.pdf", width=7.5, height=6)
j <- 0
par(mfrow=c(5, 4), mar=c(1, 1, 1, 1), oma=c(3, 4, 1, 1))
for(i in order72_MIC){
  j <- j+1
  sub_t0 <- subset(all0.72, line == i)
  sub_t0_ag <- aggregate(sub_t0["data"], sub_t0[c("enviro")], median)
  sub <- subset(all10.72, line == i)
  sub$enviro[sub$enviro == 0] <- 0.0625
  sub_ag <-  aggregate(sub["dataM"], sub[c("enviro", "rep", "lr", "col", "lwd")], median)
  sub_ag1 <- subset(sub_ag, rep == 1)
  plot(log(sub_ag1$enviro), sub_ag1$dataM, type="l", yaxt="n", xaxt="n", ylim=c(0, 2), lwd=0.5, col="grey")
  abline(v=log(fitFlow.ag$MIC24[as.numeric(i)]), col="grey", lwd=2)
  for (k in 1:12){
    sub_ag_sub <- subset(sub_ag, rep == k)
    points(log(sub_ag_sub$enviro), sub_ag_sub$dataM, type="l", col= sub_ag_sub$col, lwd = sub_ag_sub$lwd)
    #points(log(sub_ag_sub$enviro), sub_ag_sub$dataM, type="l", col= ifelse(sub_ag_sub$lr[1] %in% outliers, "black", "black"), lwd = ifelse(sub_ag_sub$lr[1] %in% outliers, 3, 0.5))
  }
  points(log(sub_t0_ag$enviro), sub_t0_ag$data, type="l", lwd=3, col="black")
  if(j %% 4 ==1) axis(2, las=2)
  else axis(2, labels=FALSE)
  if (j > 16) axis(1, at = log(c(0.0625, 0.25, 1, 4, 16, 64)), labels=c("0", "0.25", "1", "4", "16", "64"))
  else axis(1, at = log(c(0.0625, 0.25, 1, 4, 16, 64)), labels=FALSE)

  #abline(h=subset(sub_t0_ag, enviro=="0.0625")$data/2, col="blue")
  mtext(paste0("A", sub$line[1]), side=3, adj=0.01)
}
mtext(expression("Fluconazole concentration ("~mu~"g)"), side=1, line=1.5, outer=TRUE)
#mtext("Fluconazole concentration (ug)", side=1, outer=TRUE, line=2)
mtext("Optical density (24 h)", side=2, outer=TRUE, line=2)
dev.off()
system("open manuscript/figures/FigureS3-MICTraces_72h.pdf")

############################################################
#Figure 3b
#tolerance by replicate
############################################################
table(MIC24.less1$strain)
table(MIC24.1$strain)
table(MIC24.more1$strain)
pdf("manuscript/figures/Figure3b-SMG72-aboveMIC-flip.pdf", width=7.5, height=4.5)
plot(rep(1:12, each = 12), MIC24.less1$SMG72.10, yaxt="n", xaxt="n", xlab="Strain", ylab="", cex=1.2, pch=19, xlim=c(1, 22), col = MIC24.less1$col72a, ylim=c(0,1))
points(MIC24.1$place72_MIC+2, MIC24.1$SMG72.10, cex=1.2, col=MIC24.1$col72a, pch=19)
points(MIC24.more1$place72_MIC+2, MIC24.more1$SMG72.10, cex=1.2, col=MIC24.more1$col72a, pch=19)
abline(v=13, lty=2)
abline(v=17, lty=2)

axis(1, c(1:12, 14:16, 18:22), labels=c(paste0("A", c(MIC24.ag.less1$line, MIC24.ag.1$line, MIC24.ag.more1$line))), cex.axis=0.8)
axis(1, c(2, 4, 6, 8, 10, 12, 15, 19, 21) , labels=c(paste0("A", c(MIC24.ag.less1$line, MIC24.ag.1$line, MIC24.ag.more1$line)[c(2, 4, 6,8, 10, 12, 14, 17, 19)])),  cex.axis=0.8)

points(1:12, MIC24.ag.less1$SMG72.down, pch="-", col="black", cex=2)
points(1:12, MIC24.ag.less1$SMG72.up, pch="-", col="black", cex=2)
arrows(1:12, MIC24.ag.less1$SMG72.up, 1:12, MIC24.ag.less1$SMG72.down, length=0, col="black")
points(14:16, MIC24.ag.1$SMG72.up, pch="-", col="black", cex=2)
points(14:16, MIC24.ag.1$SMG72.down, pch="-", col="black", cex=2)
arrows(14:16, MIC24.ag.1$SMG72.up, 14:16, MIC24.ag.1$SMG72.down, length=0, col="black")
points(18:22, MIC24.ag.more1$SMG72.up, pch="-", col="black", cex=2)
points(18:22, MIC24.ag.more1$SMG72.down, pch="-", col="black", cex=2)
arrows(18:22, MIC24.ag.more1$SMG72.up, 18:22, MIC24.ag.more1$SMG72.down, length=0, col="black")
axis(2, las=2)
mtext("Tolerance above MIC (72h)", side=2, line=3)
dev.off()
system("open manuscript/figures/Figure3b-SMG72-aboveMIC-flip.pdf")

######################################################################
#Figure S4: Change in evolved tolerance versus change in evolved fitness
######################################################################
pdf("manuscript/figures/FigureS4-DeltaEvolFit_72h-DeltaEvolSMG.pdf", width=7, height=5)
k<-0
par(mar=c(1, 1, 1, 1), mfrow=c(5, 4), mar=c(1, 1, 1, 1), oma=c(4, 4, 1, 1))
for(i in c(order72_MIC[1:19])){
  k <- k+1
  plot(subset(fitFlow, line==i)$delta1.72, subset(fitFlow, line==i)$change72.SMG, xlim=c(-0.1, 1.5), ylim=c(-1, 1), xaxt="n", yaxt="n", pch=21)
  t <- cor.test(subset(fitFlow, line==i)$delta1.72, subset(fitFlow, line==i)$change72.SMG, method="spearman")
  if(!is.na(t$p.value)){
    if(t$p.value < 0.05) abline(lm(subset(fitFlow, line==i)$change72.SMG~subset(fitFlow, line==i)$delta1.72), col="black")
    if(t$p.value > 0.05) abline(lm(subset(fitFlow, line==i)$change72.SMG~subset(fitFlow, line==i)$delta1.72), col="black", lty=2)
  }
  if(k > 16) axis(1)
  else axis(1, labels=FALSE)
  if(k %% 4 == 1) axis(2, las=2)
  else axis(2, labels=FALSE)
  text(1, 0.8, paste0("A", i), pos=4, font=2, cex=1.1)
}
mtext("Change in fitness (OD at 72 h)", side=1, outer=TRUE, line=2)
mtext("Change in tolerance (72 h)", side=2, line=2, outer = TRUE)
#A12 plotted separately because breaks the loop (reps have MIC at max)
plot(subset(fitFlow, line==12)$delta1.72, subset(fitFlow, line==12)$change72.SMG, xlim=c(-0.1, 1.5), ylim=c(-1, 1), xaxt="n", yaxt="n", pch=21, bg=ifelse(subset(fitFlow, line ==12)$MIC24.10 ==1 , "blue", "white"))
axis(1)
axis(2, labels=FALSE)
text(1, 0.8, paste0("A", 12), pos=4, font=2, cex=1.1)
dev.off()
system("open manuscript/figures/FigureS4-DeltaEvolFit_72h-DeltaEvolSMG.pdf")

################
#Ploidy
################
pdf("manuscript/figures/Figure4ab-genomeSize.pdf", width=3.25, height=7)
par(mfrow=c(2, 1), mar=c(1, 1, 1, 1), oma=c(3, 4, 1, 1))
plot(c(rep(1, 20), rep(2, 20)), c(flow.ag$t0.G1.mu, flow.ag$t10.G1.1), xlim=c(0.75, 2.25), xaxt="n", yaxt="n", ylim=c(140, 260), col=fitFlow.ag$col72a, pch=19)
axis(2, las=2)
for(i in 1:20){
  if (i %in% c(4, 9, 12))  lines(c(1, 2), c(flow.ag$t0.G1.mu[i], flow.ag$t10.G1.1[i]), lty=2)
  else lines(c(1, 2), c(flow.ag$t0.G1.mu[i], flow.ag$t10.G1.1[i]))
}
points(c(rep(1, 20), rep(2, 20)), c(flow.ag$t0.G1.mu, flow.ag$t10.G1.1), col=fitFlow.ag$col72a, pch=19)

axis(1, c(1, 2), labels=FALSE)
mtext("median genome size \n (FITC intensity)", side=2, line=2.75)
mtext("a" , side=3, adj=0.01, font=2, cex=1)
mtext("   Genome size", side=3, adj=0.1)


plot(c(rep(1, 20), rep(2, 20)), c(flow.ag.cv$t0.G1.mu, flow.ag.cv$t10.G1.1), xlim=c(0.75, 2.25), xaxt="n", yaxt="n", ylim=c(0, 0.4), col=fitFlow.ag$col72a, pch=19)
axis(2, las=2)
for(i in 1:20){
  if (i %in% c(1, 5, 12, 18))  lines(c(1, 2), c(flow.ag.cv$t0.G1.mu[i], flow.ag.cv$t10.G1.1[i]), lty=2)
  else lines(c(1, 2), c(flow.ag.cv$t0.G1.mu[i], flow.ag.cv$t10.G1.1[i]))
}
points(c(rep(1, 20), rep(2, 20)), c(flow.ag.cv$t0.G1.mu, flow.ag.cv$t10.G1.1), col=fitFlow.ag$col72a, pch=19)
axis(1, c(1, 2), labels=FALSE)
axis(1, c(1, 2), labels=c("ancestral", "evolved"))
mtext("b" , side=3, adj=0.01, font=2, cex=1)
mtext("Coefficient of variation \n (FITC intensity)", side=2, line=2.75)
mtext("   Genome size variation", side=3, adj=0.2)
dev.off()
system("open manuscript/figures/Figure4ab-genomeSize.pdf")

#################################
#(non)correlations with genome size
#################################

pdf("manuscript/figures/FigureS5-evolPloidy.pdf", width=7, height=6)
par(mfrow=c(5, 4))
j <- 0
for(i in order72_MIC){
  j <- j+1
  sub <- subset(fitFlow, strain == i)
  t0range <- range(sub$t0.G1.1, na.rm=TRUE)
  plot(sub$rep, sub$t10.G1.1, ylim= c(125, 475), yaxt="n", xaxt="n", pch=19)
  points(sub$rep, sub$t10.G1.2)
  abline(h=t0range[1], col="grey")
  abline(h=t0range[2], col="grey")
  if(j > 16) axis(1, 1:12, cex.axis=0.8)
  else axis(1, 1:12, labels=FALSE)
  if(j %% 4 == 1) axis(2, las=2, cex.axis=0.8, at=c(150, 250, 350, 450))
  else axis(2, labels=FALSE, at=c(150, 250, 350, 450))
  mtext(paste0("A", i), side=3, adj=0.01)
}
mtext("Replicate", side=1, outer=TRUE, line=2)
mtext("Genome size (FITC intensity)", side=2, outer=TRUE, line=2)
system("open manuscript/figures/FigureS5-evolPloidy.pdf")


#################################
#Figure S6
#################################
pdf("manuscript/figures/FigureS6-deltaFit72-evolPloidy.pdf", width=7, height=6)
k<-0
par(mar=c(1, 1, 1, 1), mfrow=c(5, 4), mar=c(1, 1, 1, 1), oma=c(4, 4, 1, 1))
for(i in c(order72_MIC[1:20])){
  k <- k+1
  plot(subset(fitFlow, line==i)$delta1.72, subset(fitFlow, line==i)$t10.G1.1, ylim=c(100, 500), xlim=c(-0.1, 1.4), xaxt="n", yaxt="n", pch=subset(fitFlow, line==i)$pch, bg=ifelse(subset(fitFlow, line ==i)$MIC24.10 >1, "black", "white"))
  if(k < 16){
  t <- cor.test(subset(fitFlow, line==i)$delta1.72, subset(fitFlow, line==i)$t10.G1.1, method="spearman")
  if(t$p.value < 0.05) abline(lm(subset(fitFlow, line==i)$t10.G1.1~subset(fitFlow, line==i)$delta1.72), col="red")
  if(t$p.value > 0.05) abline(lm(subset(fitFlow, line==i)$t10.G1.1~subset(fitFlow, line==i)$delta1.72), col="red", lty=2)
  }
  if(k > 12) axis(1)
  else axis(1, labels=FALSE)
  if(k %% 4 == 1) axis(2, las=2)
  else axis(2, labels=FALSE)
  text(1, 450, paste0("A", i), pos=4, font=2, cex=1.1)
}
mtext("Change in fitness (OD at 72 h)", side=1, outer=TRUE, line=2)
mtext("Evolved genome size (FITC intensity)", side=2, outer=TRUE, line=2)
dev.off()
system("open manuscript/figures/FigureS6-deltaFit72-evolPloidy.pdf")

#################################
#Figure S7
#################################
pdf("manuscript/figures/FigureS7-deltaSMG-evolPloidy.pdf", width=7, height=6)
k<-0
par(mar=c(1, 1, 1, 1), mfrow=c(5, 4), mar=c(1, 1, 1, 1), oma=c(4, 4, 1, 1))
for(i in c(order72_MIC[1:20])){
  k <- k+1
  plot(subset(fitFlow, line==i)$change72.SMG, subset(fitFlow, line==i)$t10.G1.1, ylim=c(100, 500), xlim=c(-0.8, 0.8), xaxt="n", yaxt="n", pch=subset(fitFlow, line==i)$pch)
  if(k < 16){
    t <- cor.test(subset(fitFlow, line==i)$change72.SMG, subset(fitFlow, line==i)$t10.G1.1, method="spearman")
    if(t$p.value < 0.05) abline(lm(subset(fitFlow, line==i)$t10.G1.1~subset(fitFlow, line==i)$change72.SMG), col="red")
    if(t$p.value > 0.05) abline(lm(subset(fitFlow, line==i)$t10.G1.1~subset(fitFlow, line==i)$change72.SMG), col="red", lty=2)
  }
  if(k > 16) axis(1)
  else axis(1, labels=FALSE)
  if(k %% 4 == 1) axis(2, las=2)
  else axis(2, labels=FALSE)
  text(-0.7, 450, paste0("A", i), pos=4, font=2, cex=1.1)
}
mtext("Change in tolerance (SMG at 72 h)", side=1, outer=TRUE, line=2)
mtext("Evolved genome size (FITC intensity)", side=2, outer=TRUE, line=2)
dev.off()
system("open manuscript/figures/FigureS7-deltaSMG-evolPloidy.pdf")

