library(viridis)

#All disk assay and flow data
se <- function(x, na.rm=TRUE) sqrt(var(x, na.rm=TRUE)/(length(x) - 1))
source("scripts/02_MIC-data.R")
source("scripts/03_flow-data.R")

#load in ancestral strain information
lines <- read.csv("data_in/general/LDE-strainTable.csv")
lines$init_mic <- c(4, 1, 0.0125, 0.5, 4, 1, 1, 0.5, 0.0125, 0.0125, 0.0125, 32, 0.0125, 0.0125, 0.0125, 0.0125, 0.0125, 4, 0.0125, 4)

#Showsthe experiments where data was collected (table S3)
BMDtab.0 <- data.frame(exper=1:6, "YPD" = rep("x",6), "0.5" = c("", "", "x", "", "", ""), "1" =c("x", "x", "x", "x", "", ""), "4" = c("", "", "x", "x", "", ""), "8" = c("x", "", "", "x", "x", ""), "32" = c("x", "", "", "x", "x", "x"), "128" = c("", "x", "", "", "x", "x"), "512" = c("", "", "", "", "", "x"))

BMDtab.10 <- data.frame(exper=c(1, 3, 4, 5), "YPD" = rep("x",4), "0.5" = c("", "x", "", ""), "1" =c("x", "x", "", ""), "4" = c("",  "x", "", ""), "8" = c("x","", "x", ""), "32" = c("x", "", "x", "x"), "128" = c( "", "", "x", "x"), "512" = c("", "", "", "x"))

#colours based on initial MIC
#colours
coloursM <- magma(5, direction=-1)
coloursV <- viridis(5, direction=-1)
lines$colours[lines$init_mic == 0.0125] <- coloursM[2]
lines$colours[lines$init_mic == 0.5] <- coloursV[2]
lines$colours[lines$init_mic == 1] <- coloursV[4]
lines$colours[lines$init_mic == 4] <- coloursM[3]
lines$colours[lines$init_mic == 32] <- coloursV[5]

plot(lines$strain, log(lines$init_mic), col=lines$colours, pch=21, cex=2)
#################################################
#order the strains based first on fitness in 1ug FLC
#################################################
#Single drug
all0.1 <- filter(all0.rep, enviro==1)
all0.ag.1 <- aggregate(all0.1["data"], all0.1["line"], mean)
all0.ag.1$sd <- aggregate(all0.1["data"], all0.1["line"], sd)$data
all10.1 <- filter(all10.rep, enviro==1)
all10.ag.1 <- aggregate(all10.1["data"], all10.1["line"], mean)
all10.ag.1$sd <- aggregate(all10.1["data"], all10.1["line"], sd)$data

all0.1.72 <- filter(all0.rep.72, enviro==1)
all0.ag.1.72 <- aggregate(all0.1.72["data"], all0.1.72["line"], mean)
all0.ag.1$sd <- aggregate(all0.1.72["data"], all0.1.72["line"], sd)$data
all10.1.72 <- filter(all10.rep.72, enviro==1)
all10.ag.1.72 <- aggregate(all10.1.72["data"], all10.1.72["line"], mean)
all10.ag.1$sd <- aggregate(all10.1.72["data"], all10.1.72["line"], sd)$data

place <- order(all0.ag.1$data)

all0.0 <- filter(all0.rep, enviro==0)
all0.ag.0 <- aggregate(all0.0["data"], all0.0["line"], mean)
all0.ag.0$sd <- aggregate(all0.0["data"], all0.0["line"], sd)$data
all10.0 <- filter(all10.rep, enviro==0)
all10.ag.0 <- aggregate(all10.0["data"], all10.0["line"], mean)
all10.ag.0$sd <- aggregate(all10.0["data"], all10.0["line"], sd)$data
all0.0.72 <- filter(all0.rep.72, enviro==0)
all0.ag.0.72 <- aggregate(all0.0.72["data"], all0.0.72["line"], mean)
all0.ag.0.72$sd <- aggregate(all0.0.72["data"], all0.0.72["line"], sd)$data
all10.0.72 <- filter(all10.rep.72, enviro==0)
all10.ag.0.72 <- aggregate(all10.0.72["data"], all10.0.72["line"], mean)
all10.ag.0.72$sd <- aggregate(all10.0.72["data"], all10.0.72["line"], sd)$data

newdf3 <- all10.1[1,]
newdf3$place <- 0
newdf3$col <- "#999999"
for(i in 1:20){
  sub <- subset(all10.1, line==i)
  sub$place <- which(place==i)
  sub$col <- rep(subset(lines, strain == i)$colours, 12)
  newdf3 <- bind_rows(newdf3, sub)
}
all10.1 <- newdf3[-1,]
all0.1 <- cbind(all0.1, place=all10.1$place, col=all10.1$col)
all10.1.72 <-  cbind(all10.1.72, place=all10.1$place, col=all10.1$col)
all0.0 <-  cbind(all0.0, place=all10.1$place, col=all10.1$col)
all10.0 <-  cbind(all10.0, place=all10.1$place, col=all10.1$col)
all0.0.72 <-  cbind(all0.0.72, place=all10.1$place, col=all10.1$col)
all10.0.72 <-  cbind(all10.0.72, place=all10.1$place, col=all10.1$col)

##########################
#RESISTANCE TRAITS
##########################
#add place into the MICall dataframe
newdf <- MICall[1,]
newdf$place <- 0
newdf$col <- "#000000"
for(i in 1:20){
  sub <- subset(MICall, line==i)
  sub$place <- which(place==i)
  sub$col <- rep(subset(lines, strain == i)$colours, 12)
  newdf <- rbind(newdf, sub)
}
MICall <- newdf[-1,]
MICall$lr <-  paste(MICall$line, MICall$rep, sep=".")

##################################
#FLOW
##################################
flow$t0.G1.mu <- apply(flow[,2:3], 1, mean)
ddn.flow <- split(flow, flow$line)

#Find the median of all A12 reps in each timepoine
mu0 <- median(ddn.flow[[12]]$t0.G1.mu)
mu3 <- median(ddn.flow[[12]]$t3.G1.1)
mu5 <- median(ddn.flow[[12]]$t5.G1.1)
mu8 <- median(ddn.flow[[12]]$t8.G1.1)
mu10 <- median(ddn.flow[[12]]$t10.G1.1)
mu10b <- median(ddn.flow[[12]]$t10b.G1.1, na.rm=TRUE)

#Use for correlation
# cor0 <- mu0/mu0b
cor3 <- mu0/mu3
cor5 <- mu0/mu5
cor8 <- mu0/mu8
cor10 <- mu0/mu10
cor10b <- mu0/mu10b

flow.ag <- aggregate(flow[c("t0.G1.mu", "t3.G1.1", "t5.G1.1", "t8.G1.1", "t10.G1.1", "t10b.G1.1")], flow["line"], median, na.rm=TRUE)
flow.ag.sd <- aggregate(flow[c("t0.G1.mu", "t3.G1.1", "t5.G1.1", "t8.G1.1", "t10.G1.1", "t10b.G1.1")], flow["line"], sd, na.rm=TRUE)
#flow.ag.cv <- aggregate(flow[c("t0.G1.mu", "t3.G1.1", "t5.G1.1", "t8.G1.1", "t10.G1.1", "t10b.G1.1")], flow["line"], cv)

flow.ag$t3.G1.1 <- flow.ag$t3.G1.1*cor3
flow.ag$t5.G1.1 <- flow.ag$t5.G1.1*cor5
flow.ag$t8.G1.1 <- flow.ag$t8.G1.1*cor8
flow.ag$t10.G1.1 <- flow.ag$t10.G1.1*cor10
flow.ag$t10b.G1.1 <- flow.ag$t10b.G1.1*cor10b

############################################################################
#Join all into one dataframe
#can not find the clade of FH1 to save my life so put in a dummy clade for now
#can not add disk assay data here because there is missing data
############################################################################
fitFlow <- data.frame(strain = rep(lines$strain, each=12), clade = rep(lines$clade, each=12), MTL = rep(lines$MTL, each=12), zygosity = rep(lines$zygosity, each=12), all0.1, all0.1.72[,4:6], all10.1[, 4:6], all10.1.72[,4:6], "t0.G1.1" = flow$t0.G1.mu,"t10.G1.1" =  flow$t10.G1.1*cor10, "t10.G1.2" = flow$t10.G1.2*cor10, "t10.G1.3" = flow$t10.G1.3*cor10,  MICall[, c(2:19)])

fitFlow.ddn <- split(fitFlow, fitFlow$strain)
SMG24.down <- unlist(lapply(fitFlow.ddn, function(x) range(x$SMG24.2, na.rm=TRUE)[1]))
SMG24.up <- unlist(lapply(fitFlow.ddn, function(x) range(x$SMG24.2, na.rm=TRUE)[2]))
SMG48.down <- unlist(lapply(fitFlow.ddn, function(x) range(x$SMG48.2, na.rm=TRUE)[1]))
SMG48.up <- unlist(lapply(fitFlow.ddn, function(x) range(x$SMG48.2, na.rm=TRUE)[2]))
SMG72.down <- unlist(lapply(fitFlow.ddn, function(x) range(x$SMG72.2, na.rm=TRUE)[1]))
SMG72.up <- unlist(lapply(fitFlow.ddn, function(x) range(x$SMG72.2, na.rm=TRUE)[2]))

names(fitFlow)[8:10] <- c("all0.ag.1", "all0.sd.1", "all0.se.1")
names(fitFlow)[13:15] <- c("all0.ag.1.72", "all0.sd.1.72", "all0.se.1.72")
names(fitFlow)[16:18] <- c("all10.ag.1", "all10.sd.1", "all10.se.1")
names(fitFlow)[19:21] <- c("all10.ag.1.72", "all10.sd.1.72", "all10.se.1.72")
fitFlow$lr <- paste(fitFlow$line, fitFlow$rep, sep=".")
fitFlow$pch <- ifelse(is.na(fitFlow$t10.G1.2), 21, 24)
fitFlow$clade <- as.factor(fitFlow$clade)
fitFlow$delta1.24 <- fitFlow$all10.ag.1-fitFlow$all0.ag.1
fitFlow$delta1.72 <- fitFlow$all10.ag.1.72-fitFlow$all0.ag.1.72
fitFlow$MIC24 <- rep(MICall.ag$MIC24, each=12) #this changes t0 to the median
fitFlow$MIC48 <- rep(MICall.ag$MIC48, each=12) #this changes t0 to the median
fitFlow$MIC72 <- rep(MICall.ag$MIC72, each=12) #this changes t0 to the median
fitFlow$SMG24 <- rep(MICall.ag$SMG24, each=12) #this changes t0 to the median
fitFlow$SMG48 <- rep(MICall.ag$SMG48, each=12) #this changes t0 to the median
fitFlow$SMG72 <- rep(MICall.ag$SMG72, each=12) #this changes t0 to the median
fitFlow$change24 <- fitFlow$MIC24.10-fitFlow$MIC24 #this changes t0 to the median
fitFlow$change48 <- fitFlow$MIC48.10-fitFlow$MIC48 #this changes t0 to the median
fitFlow$change72 <- fitFlow$MIC72.10-fitFlow$MIC72 #this changes t0 to the median
fitFlow$change24.SMG <- fitFlow$SMG24.10 - fitFlow$SMG24
fitFlow$change48.SMG <- fitFlow$SMG48.10 - fitFlow$SMG48
fitFlow$change72.SMG <- fitFlow$SMG72.10 - fitFlow$SMG72
fitFlow$change24.SMG.2 <- fitFlow$SMG24.10.2 - fitFlow$SMG24.2
fitFlow$change48.SMG.2 <- fitFlow$SMG48.10.2 - fitFlow$SMG48.2
fitFlow$change72.SMG.2 <- fitFlow$SMG72.10.2 - fitFlow$SMG72.2

fitFlow$lr <- paste(fitFlow$line, fitFlow$rep, sep="-")

#this removes the four non-variable ploidy lines
fitFlow.variable <- subset(fitFlow, line %in%  c(2:4, 6:11, 13:17, 19, 20))

use.wells <- c(8, 13, 16, 19, 22, 23, 26:43, 45:56)
#aggregates
fitFlow.ag <- aggregate(fitFlow[, use.wells], fitFlow[c("clade", "zygosity", "col", "line")], mean, na.rm=TRUE)
#fitFlow.ag$cv.t10.G1.1 <- flow.ag.cv$t10.G1.1
fitFlow.ag$change24adj <- sign(fitFlow.ag$change24)*log(abs(fitFlow.ag$change24)+1)
fitFlow.ag$change72adj <- sign(fitFlow.ag$change72)*log(abs(fitFlow.ag$change72)+1)
fitFlow.ag$change24ave <- fitFlow.ag$all10.ag.1 - fitFlow.ag$all0.ag.1
fitFlow.ag$change72ave <- fitFlow.ag$all10.ag.1.72 - fitFlow.ag$all0.ag.1.72
fitFlow.ag$SMG24.up <- SMG24.up
fitFlow.ag$SMG24.down <- SMG24.down
fitFlow.ag$SMG48.up <- SMG48.up
fitFlow.ag$SMG48.down <- SMG48.down
fitFlow.ag$SMG72.up <- SMG72.up
fitFlow.ag$SMG72.down <- SMG72.down

fitFlow.sd <- aggregate(fitFlow[, use.wells], fitFlow[c("clade", "zygosity", "col", "line")], sd, na.rm=TRUE)

fitFlow.ag.variable <- subset(fitFlow.ag, line %in% c(2:4, 6:11, 13:17, 19, 20))
fitFlow.sd.variable <- subset(fitFlow.sd, line %in% c(2:4, 6:11, 13:17, 19, 20))

###################################
#change MIC values for plotting
###################################
fitFlow.plot <- fitFlow
#change 0 to 0.001 in the variables where looking at change over evol time so can take log
fitFlow.plot$MIC24.10[fitFlow.plot$MIC24.10 < 1] <- 0.25
fitFlow.plot$MIC24.10[fitFlow.plot$MIC24.10 > 128] <- 256
fitFlow.plot$MIC48.10[fitFlow.plot$MIC48.10 < 1] <- 0.25
fitFlow.plot$MIC48.10[fitFlow.plot$MIC48.10 > 128] <- 256
fitFlow.plot$MIC72.10[fitFlow.plot$MIC72.10 < 1] <- 0.25
fitFlow.plot$MIC72.10[fitFlow.plot$MIC72.10 > 128] <- 256
fitFlow.plot$change24[fitFlow.plot$change24==0] <- 0.001
fitFlow.plot$change48[fitFlow.plot$change48==0] <- 0.001
fitFlow.plot$change72[fitFlow.plot$change72==0] <- 0.001
fitFlow.plot$pch <- ifelse(fitFlow.plot$change24 > 1, 19, 21)

#create new dataframes that will be used for plotting purposes
#change all values < 1 to 0.25 and > 128 to 256 for plotting purposes only
fitFlow.plot.ag <- fitFlow.ag
fitFlow.plot.ag$MIC24[fitFlow.plot.ag$MIC24<1] <- 0.25
fitFlow.plot.ag$MIC24.10[fitFlow.plot.ag$MIC24.10<1] <- 0.25
fitFlow.plot.ag$MIC48[fitFlow.plot.ag$MIC48<1] <- 0.25
fitFlow.plot.ag$MIC48.10[fitFlow.plot.ag$MIC48.10<1] <- 0.25
fitFlow.plot.ag$MIC48[fitFlow.plot.ag$MIC48>128] <- 256
fitFlow.plot.ag$MIC48.10[fitFlow.plot.ag$MIC48.10>128] <- 256
fitFlow.plot.ag$MIC72[fitFlow.plot.ag$MIC72<1] <- 0.25
fitFlow.plot.ag$MIC72.10[fitFlow.plot.ag$MIC72.10<1] <- 0.25
fitFlow.plot.ag$MIC72[fitFlow.plot.ag$MIC72>128] <- 256
fitFlow.plot.ag$MIC72.10[fitFlow.plot.ag$MIC72.10>128] <- 256
