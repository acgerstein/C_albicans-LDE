findMIC2 <- function(x){
  temp <- x$enviro[which.max(x$data <  x$data[1]*0.5)-1]
  ifelse(length(temp) ==1, temp, 512)
}

source("scripts/01_raw-OD_to-rep-ag.R")
#all0, all0.rep, all0.ag, all0.48, all0.rep.48, all0.ag.48 ... all10, all10.rep, all10.ag

#make a new column to combine line and rep
all0$lr <-paste(all0$line, all0$rep, sep="_")
all0.48$lr <-paste(all0.48$line, all0.72$rep, sep="_")
all0.72$lr <-paste(all0.72$line, all0.72$rep, sep="_")
all10$lr <-paste(all10$line, all10$rep, sep="_")
all10.48$lr <-paste(all10.48$line, all10.72$rep, sep="_")
all10.72$lr <-paste(all10.72$line, all10.72$rep, sep="_")

#this is not really necessary for the MIC analysis, but is needed for plotting log2(enviro)
all0$enviro[all0$enviro==0] <- 0.000125
all0.48$enviro[all0.48$enviro==0] <- 0.000125
all0.72$enviro[all0.72$enviro==0] <- 0.000125
all10$enviro[all10$enviro==0] <- 0.000125
all10.48$enviro[all10.48$enviro==0] <- 0.000125
all10.72$enviro[all10.72$enviro==0] <- 0.000125

#remove 512 information
all0 <- filter(all0, enviro!= 512)
all0.48 <- filter(all0.48, enviro!= 512)
all0.72 <- filter(all0.72, enviro!= 512)
all10 <- filter(all10, enviro!= 512)
all10.48 <- filter(all10.48, enviro!= 512)
all10.72 <- filter(all10.72, enviro!= 512)


#split by linerep
ddn.0lr <- split(all0, all0$lr)
ddn.0lr.48 <- split(all0.48, all0.48$lr)
ddn.0lr.72 <- split(all0.72, all0.72$lr)
ddn.10lr <- split(all10, all10$lr)
ddn.10lr.48 <- split(all10.48, all10.48$lr)
ddn.10lr.72 <- split(all10.72, all10.72$lr)

# Rosenberg: Supra-MIC Growth (SMG) was calculated as the average growth per well above the MIC divided by the level of growth in the well without drug.

#find MIC for each linerep at 24, 48, 72h
MICs0 <- lapply(ddn.0lr, function(x){
  ag <- aggregate(x["data"], x["enviro"], median)
  MIC <- findMIC2(ag)
   SMG <- sum(ag$data[ag$enviro>MIC])/length(ag$data[ag$enviro>MIC])
   SMG2 <- SMG/subset(ag, enviro==0.000125)$data
   SMG3 <- sum(ag$data[ag$enviro>1])/length(ag$data[ag$enviro>1])
#   SMG3 <- SMG.b/subset(ag, enviro==1)$data
  c(MIC, SMG, SMG2, SMG3)
  })

  MICs0.48 <- mapply(function(x, y){
    ag.24 <- aggregate(x["data"], x["enviro"], median)
    ag.48 <- aggregate(y["data"], y["enviro"], median)
    MIC24 <- findMIC2(ag.24)
    MIC <- findMIC2(ag.48)
     SMG <- sum(ag.48$data[ag.48$enviro>MIC24])/length(ag.48$data[ag.48$enviro>MIC24])
     SMG2 <- SMG/subset(ag.48, enviro==0.000125)$data
     SMG3 <- sum(ag.48$data[ag.48$enviro>1])/length(ag.48$data[ag.48$enviro>1])
     #SMG3 <- SMG.b/subset(ag.48, enviro==1)$data
    c(MIC, SMG, SMG2, SMG3)
    }, x = ddn.0lr, y=ddn.0lr.48)

  MICs0.72 <- mapply(function(x, y){
      ag.24 <- aggregate(x["data"], x["enviro"], median)
      ag.72 <- aggregate(y["data"], y["enviro"], median)
      MIC24 <- findMIC2(ag.24)
      MIC <- findMIC2(ag.72)
      SMG <- sum(ag.72$data[ag.72$enviro>MIC24])/length(ag.72$data[ag.72$enviro>MIC24])
      SMG2 <- SMG/subset(ag.72, enviro==0.000125)$data
      SMG3 <- sum(ag.72$data[ag.72$enviro>1])/length(ag.72$data[ag.72$enviro>1])
      #SMG3 <- SMG.b/subset(ag.72, enviro==1)$data
      c(MIC, SMG, SMG2, SMG3)
      }, x = ddn.0lr, y=ddn.0lr.72)

  MICs10 <- lapply(ddn.10lr, function(x){
    ag <- aggregate(x["dataM"], x["enviro"], median)
    MIC <- findMIC2(ag)
     SMG <- sum(ag$data[ag$enviro>MIC])/length(ag$data[ag$enviro>MIC])
     SMG2 <- SMG/subset(ag, enviro==0.000125)$data
     SMG3 <- sum(ag$data[ag$enviro>1])/length(ag$data[ag$enviro>1])
#     SMG3 <- SMG.b/subset(ag, enviro==1)$data
    c(MIC, SMG, SMG2, SMG3)
    })

  MICs10.48 <- mapply(function(x, y){
        ag.24 <- aggregate(x["dataM"], x["enviro"], median)
        ag.48 <- aggregate(y["dataM"], y["enviro"], median)
        MIC24 <- findMIC2(ag.24)
        MIC <- findMIC2(ag.48)
        SMG <- sum(ag.48$data[ag.48$enviro>MIC24])/length(ag.48$data[ag.48$enviro>MIC24])
        SMG2 <- SMG/subset(ag.48, enviro==0.000125)$data
        SMG3 <- sum(ag.48$data[ag.48$enviro>1])/length(ag.48$data[ag.48$enviro>1])
 #       SMG3 <- SMG.b/subset(ag.48, enviro==1)$data
        c(MIC, SMG, SMG2, SMG3)
        }, x = ddn.10lr, y=ddn.10lr.48)

  MICs10.72 <- mapply(function(x, y){
          ag.24 <- aggregate(x["dataM"], x["enviro"], median)
          ag.72 <- aggregate(y["dataM"], y["enviro"], median)
          MIC24 <- findMIC2(ag.24)
          MIC <- findMIC2(ag.72)
          SMG <- sum(ag.72$data[ag.72$enviro>MIC24])/length(ag.72$data[ag.72$enviro>MIC24])
          SMG2 <- SMG/subset(ag.72, enviro==0.000125)$data
          SMG3 <- sum(ag.72$data[ag.72$enviro>1])/length(ag.72$data[ag.72$enviro>1])
 #         SMG3 <- SMG.b/subset(ag.72, enviro==1)$data
          c(MIC, SMG, SMG2, SMG3)
          }, x = ddn.10lr, y=ddn.10lr.72)


#find MIC for each linerep at 24, 48, 72h
MICall <- data.frame(name= as.character(names(MICs0)), MIC24 =  matrix(unlist(MICs0),ncol=4,byrow=TRUE)[,1], MIC48 = MICs0.48[1,], MIC72 = MICs0.72[1,], MIC24.10 = matrix(unlist(MICs10),ncol=4,byrow=TRUE)[,1], MIC48.10 = MICs10.48[1,], MIC72.10 = MICs10.72[1,], SMG24 =  matrix(unlist(MICs0),ncol=4,byrow=TRUE)[,2], SMG48 = MICs0.48[2,], SMG72 = MICs0.72[2,], SMG24.10 = matrix(unlist(MICs10),ncol=4,byrow=TRUE)[,2], SMG48.10 = MICs10.48[2,], SMG72.10 = MICs10.72[2,], SMG24.2 =  matrix(unlist(MICs0),ncol=4,byrow=TRUE)[,3], SMG48.2 = MICs0.48[3,], SMG72.2 = MICs0.72[3,], SMG24.10.2 = matrix(unlist(MICs10),ncol=4,byrow=TRUE)[,3], SMG48.10.2 = MICs10.48[3,], SMG72.10.2 = MICs10.72[3,], SMG24.3 =  matrix(unlist(MICs0),ncol=4,byrow=TRUE)[,4], SMG48.3 = MICs0.48[4,], SMG72.3 = MICs0.72[4,], SMG24.10.3 = matrix(unlist(MICs10),ncol=4,byrow=TRUE)[,4], SMG48.10.3 = MICs10.48[4,], SMG72.10.3 = MICs10.72[4,])
MICall$name <- as.character(MICall$name)
MICall$line <- unlist(lapply(MICall$name, function(x) strsplit(x, "_")[[1]][1]))
MICall$rep <- unlist(lapply(MICall$name, function(x) strsplit(x, "_")[[1]][2]))
MICall <- MICall[order(factor(MICall$line, 1:20), factor(MICall$rep, 1:12)),]

MICall$change24 <- MICall$MIC24.10 - MICall$MIC24
MICall$change48 <- MICall$MIC48.10 - MICall$MIC48
MICall$change72 <- MICall$MIC72.10 - MICall$MIC72
MICall.ag <- aggregate(MICall[c("MIC24", "MIC48", "MIC72", "MIC24.10", "MIC48.10", "MIC72.10", "change24", "change48", "change72", "SMG24", "SMG48", "SMG72", "SMG24.10", "SMG48.10", "SMG72.10", "SMG24.2", "SMG48.2", "SMG72.2", "SMG24.10.2", "SMG48.10.2", "SMG72.10.2", "SMG24.3", "SMG48.3", "SMG72.3", "SMG24.10.3", "SMG48.10.3", "SMG72.10.3")], MICall[c("line")], median, na.rm=TRUE)

MICall.ag <- MICall.ag[order(factor(MICall.ag$line, 1:20)),]

write.csv(MICall, "tables_intermediate/MIC_all_rep.csv", row.names=FALSE)
write.csv(MICall.ag, "tables_intermediate/MIC_ag.csv", row.names=FALSE)
