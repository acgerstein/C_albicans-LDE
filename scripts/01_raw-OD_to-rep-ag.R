#########################
#STUFF
#########################
library(dplyr)
library(reshape2)
se <- function(x) sqrt(var(x)/(length(x) - 1))

####################################################################################################
######################################
#Load in the datasets
#gen = t0a and t10 from MA4, gen = t0, t3, t13 from MA5
#######################################
##############################
#LIQUID ASSAY MA4a - evolved in 1ug, assayed 0, 1, 4, 8, 32
##############################
ma4.24 <- read.csv("data_in/MIC/raw/OD_MA4_24h_df.csv", sep=",", as.is=TRUE)
ma4.24.0 <- subset(ma4.24, gen== "t0")
ma4b.24 <- read.csv("data_in/MIC/raw/OD_MA4_24h_df_b.csv", sep=",", as.is=TRUE)
ma4b.24$enviro[ma4b.24$enviro==1.1] <- 1
ma4b.24.0 <- subset(ma4b.24, gen== "t0")

ma4.48 <- read.csv("data_in/MIC/raw/OD_MA4_48h_df.csv", sep=",", as.is=TRUE)
ma4.48.0 <- subset(ma4.48, gen== "t0")
ma4b.48 <- read.csv("data_in/MIC/raw/OD_MA4_48h_df_b.csv", sep=",", as.is=TRUE)
ma4b.48$enviro[ma4b.48$enviro==1.1] <- 1
ma4b.48.0 <- subset(ma4b.48, gen== "t0")

ma4.72 <- read.csv("data_in/MIC/raw/OD_MA4_72h_df.csv", sep=",", as.is=TRUE)
ma4.72.0 <- subset(ma4.72, gen== "t0")
ma4b.72 <- read.csv("data_in/MIC/raw/OD_MA4_72h_df_b.csv", sep=",", as.is=TRUE)
ma4b.72$enviro[ma4b.72$enviro==1.1] <- 1
ma4b.72.0 <- subset(ma4b.72, gen== "t0")

##############################
#LIQUID ASSAY MA5 - evolved in 8ug, assayed 0, 1, 4, 8, 32
##############################
ma5.24 <- read.csv("data_in/MIC/raw/OD_MA5_24h_df.csv", sep=",", as.is=TRUE)
ma5.24.0 <- subset(ma5.24, gen== "t0")
ma5.48 <- read.csv("data_in/MIC/raw/OD_MA5_48h_df.csv", sep=",", as.is=TRUE)
ma5.48.0 <- subset(ma5.48, gen== "t0")
ma5.72 <- read.csv("data_in/MIC/raw/OD_MA5_72h_df.csv", sep=",", as.is=TRUE)
ma5.72.0 <- subset(ma5.72, gen== "t0")

ma5.rep.24 <- summarise(group_by(ma5.24, line, rep, enviro, gen), data=mean(data), sd=sd(data), se=se(data))
ma5.rep.48 <- summarise(group_by(ma5.48, line, rep, enviro, gen), data=mean(data), sd=sd(data), se=se(data))
ma5.rep.72 <- summarise(group_by(ma5.72, line, rep, enviro, gen), data=mean(data), sd=sd(data), se=se(data))

############################################################
#LIQUID ASSAY MA6- evolved in 32ug, assayed 0, 8, 32, 128
############################################################
ma6.24 <- read.csv("data_in/MIC/raw/OD_MA6_24h_df.csv", sep=",", as.is=TRUE)
ma6.24.0 <- subset(ma6.24, gen=="t0c")
ma6.48 <- read.csv("data_in/MIC/raw/OD_MA6_48h_df.csv", sep=",", as.is=TRUE)
ma6.48.0 <- subset(ma6.48, gen=="t0c")
ma6.72 <- read.csv("data_in/MIC/raw/OD_MA6_72h_df.csv", sep=",", as.is=TRUE)
ma6.72.0 <- subset(ma6.72, gen=="t0c")

############################################################
#LIQUID ASSAY MA7- evolved in 128ug, assayed  0, 32, 128, 512
############################################################
ma7.24 <- read.csv("data_in/MIC/raw/OD_MA7_24h_df.csv", sep=",", as.is=TRUE)
ma7.24.0 <- subset(ma7.24, gen=="t0b")
ma7.48 <- read.csv("data_in/MIC/raw/OD_MA7_48h_df.csv", sep=",", as.is=TRUE)
ma7.48.0 <- subset(ma7.48, gen=="t0b")
ma7.72 <- read.csv("data_in/MIC/raw/OD_MA7_72h_df.csv", sep=",", as.is=TRUE)
ma7.72.0 <- subset(ma7.72, gen=="t0b")

##############################
#LIQUID ASSAY MA47 - evolved in 1ug, t3 and MA7 - evolved in 128
##############################
ma47.24 <- read.csv("data_in/MIC/raw/OD_MA47_24h_df.csv", sep=",", as.is=TRUE)
ma47.24 <- filter(ma47.24, enviro!="512")
ma47.24.0 <- subset(ma47.24, gen=="t0d")
ma47.48 <- read.csv("data_in/MIC/raw/OD_MA47_48h_df.csv", sep=",", as.is=TRUE)
ma47.48 <- filter(ma47.48, enviro!="512")
ma47.48.0 <- subset(ma47.48, gen=="t0d")
ma47.72 <- read.csv("data_in/MIC/raw/OD_MA47_72h_df.csv", sep=",", as.is=TRUE)
ma47.72 <- filter(ma47.72, enviro!="512")
ma47.72.0 <- subset(ma47.72, gen=="t0d")

ma4t3 <- subset(ma47.24, gen %in% c("t0d", "t03") & enviro==1)
ma4t3.ag <- summarise(group_by(ma47.24, gen, line),data= mean(data), sd = sd(data), se= se(data))
ma4t3.48 <- subset(ma47.48, gen %in% c("t0d", "t03") & enviro==1)
ma4t3.ag.48 <- summarise(group_by(ma47.48, gen, line),data= mean(data), sd = sd(data), se= se(data))
ma4t3.72 <- subset(ma47.72, gen %in% c("t0d", "t03") & enviro==1)
ma4t3.ag.72 <- summarise(group_by(ma47.72, gen, line),data= mean(data), sd = sd(data), se= se(data))

##############################
#LIQUID ASSAY MA445 - MA4 t0, 3, 10 evolved in 1ug, MA4 t0, t10 evolved in 8, MA5 - evolved in 8
##############################
ma45.all <- read.csv("data_in/MIC/raw/OD_MA4MA5_alltp_df.csv", as.is=TRUE)

ma45t0 <- subset(ma45.all, gen == "MA4_t0")
ma45t0.ag <- summarise(group_by(ma45t0, line, enviro), data24= mean(OD.24), data48= mean(OD.48), sd24 = sd(OD.24), se24= se(OD.24), sd48 = sd(OD.48), se48= se(OD.48))

ma5bt0 <- subset(ma45.all, gen=="MA5_t0")
ma5bt0.ag <- summarise(group_by(ma5bt0, line),data24= mean(OD.24), data48= mean(OD.48), sd24 = sd(OD.24), se24= se(OD.24), sd48 = sd(OD.48), se48= se(OD.48))
ma5bt10 <- subset(ma45.all, gen=="MA5_t10")
ma5bt10.ag <- summarise(group_by(ma5bt10, line),data24= mean(OD.24), data48= mean(OD.48), sd24 = sd(OD.24), se24= se(OD.24), sd48 = sd(OD.48), se48= se(OD.48))


####################################################################################
#Compare t0s
#gen = t0a and t10 from MA4, gen = t0, t3, t13 from MA5

#AGGREGATE AMONG THE SIX DATASETS
####################################################################

#######
#24h
#######
ma4.0.rep <- data.frame(summarise(group_by(ma4.24.0, line, col, enviro), mean(data), exper="4"))
names(ma4.0.rep)[2] <- "rep"
ma4b.0.rep <- data.frame(summarise(group_by(ma4b.24.0, line, col, enviro), mean(data), exper="4b"))
names(ma4b.0.rep)[2] <- "rep"
ma5.0.rep <- data.frame(summarise(group_by(ma5.24.0, line, rep, enviro), mean(data), exper="5"))
ma6.0.rep <- data.frame(summarise(group_by(ma6.24.0, line, rep, enviro), mean(data), exper="6"))
ma7.0.rep <- data.frame(summarise(group_by(ma7.24.0, line, rep, enviro), mean(data), exper="7"))
ma47.0.rep <- data.frame(summarise(group_by(ma47.24.0, line, rep, enviro), mean(data), exper="47"))

#all0 <- rbind(ma4.0.rep, ma4b.0.rep, ma5.0.rep, ma6.0.rep, ma7.0.rep, ma47.0.rep)
all0 <- rbind(ma4.0.rep, ma4b.0.rep, ma5.0.rep, ma6.0.rep, ma47.0.rep)
names(all0)[4] <- "data"
write.csv(all0, "data_in/MIC/aggregated/OD_all0_24.csv",  row.names=FALSE)

#average across all experiments for each of the 12 lines (= then 12 replicates)
all0.rep <- summarise(group_by(all0, line, rep, enviro), data = mean(data), sd = sd(data), se= se(data))
write.csv(all0.rep, "data_in/MIC/aggregated/OD_all0_24_rep.csv",  row.names=FALSE)

ma4.0.ag <- summarise(group_by(ma4.24.0, line, enviro), mean(data), rep="4")
ma4b.0.ag <- summarise(group_by(ma4b.24.0, line, enviro), mean(data), rep="4b")
ma5.0.ag <- summarise(group_by(ma5.24.0, line, enviro), mean(data), rep="5")
ma6.0.ag <- summarise(group_by(ma6.24.0, line, enviro), mean(data), rep="6")
ma7.0.ag <- summarise(group_by(ma7.24.0, line, enviro), mean(data), rep="7")
ma47.0.ag <- summarise(group_by(ma47.24.0, line, enviro), mean(data), rep="47")

#all0ag <- data.frame(rbind(ma4.0.ag, ma4b.0.ag, ma5.0.ag, ma6.0.ag, ma7.0.ag, ma47.0.ag))
all0ag <- data.frame(rbind(ma4.0.ag, ma4b.0.ag, ma5.0.ag, ma6.0.ag, ma47.0.ag))
names(all0ag)[3] <- "data"
all0.ag <- summarise(group_by(all0ag, line, enviro), data = mean(data), sd = sd(data), se= se(data))
write.csv(all0.ag, "data_in/MIC/aggregated/OD_all0_24_ag.csv",  row.names=FALSE)

#######
#48h
#######
ma4.0.48.rep <- data.frame(summarise(group_by(ma4.48.0, line, col, enviro), mean(data), exper="4"))
names(ma4.0.48.rep)[2] <- "rep"
ma4b.0.48.rep <- data.frame(summarise(group_by(ma4b.48.0, line, col, enviro), mean(data), exper="4b"))
names(ma4b.0.48.rep)[2] <- "rep"
ma5.0.48.rep <- data.frame(summarise(group_by(ma5.48.0, line, rep, enviro), mean(data), exper="5"))
ma6.0.48.rep <- data.frame(summarise(group_by(ma6.48.0, line, rep, enviro), mean(data), exper="6"))
ma7.0.48.rep <- data.frame(summarise(group_by(ma7.48.0, line, rep, enviro), mean(data), exper="7"))
ma47.0.48.rep <- data.frame(summarise(group_by(ma47.48.0, line, rep, enviro), mean(data), exper="47"))

#all0.48 <- data.frame(rbind(ma4.0.48.rep, ma4b.0.48.rep, ma5.0.48.rep, ma6.0.48.rep, ma7.0.48.rep, ma47.0.48.rep))
all0.48 <- data.frame(rbind(ma4.0.48.rep, ma4b.0.48.rep, ma5.0.48.rep, ma6.0.48.rep,ma47.0.48.rep))
names(all0.48)[4] <- "data"
write.csv(all0.48, "data_in/MIC/aggregated/OD_all0_48.csv",  row.names=FALSE)

all0.rep.48 <- summarise(group_by(all0.48, line, rep, enviro), data = mean(data), sd = sd(data), se = se(data))
write.csv(all0.rep.48, "data_in/MIC/aggregated/OD_all0_48_rep.csv",  row.names=FALSE)

ma4.0.48.ag <- summarise(group_by(ma4.48.0, line, enviro), mean(data), rep="4")
ma4b.0.48.ag <- summarise(group_by(ma4b.48.0, line, enviro), mean(data), rep="4b")
ma5.0.48.ag <- summarise(group_by(ma5.48.0, line, enviro), mean(data), rep="5")
ma6.0.48.ag <- summarise(group_by(ma6.48.0, line, enviro), mean(data), rep="6")
ma7.0.48.ag <- summarise(group_by(ma7.48.0, line, enviro), mean(data), rep="7")
ma47.0.48.ag <- summarise(group_by(ma47.48.0, line, enviro), mean(data), rep="47")

#all0.48ag <- data.frame(rbind(ma4.0.48.ag, ma4b.0.48.ag, ma5.0.48.ag, ma6.0.48.ag, ma7.0.48.ag, ma47.0.48.ag))
all0.48ag <- data.frame(rbind(ma4.0.48.ag, ma4b.0.48.ag, ma5.0.48.ag, ma6.0.48.ag,  ma47.0.48.ag))
names(all0.48ag)[3] <- "data"

all0.ag.48 <- summarise(group_by(all0.48ag, line, enviro), data = mean(data), sd = sd(data), se = se(data))
write.csv(all0.ag.48, "data_in/MIC/aggregated/OD_all0_48_ag.csv",  row.names=FALSE)

#######
#72h
#######
ma4.0.72.rep <- data.frame(summarise(group_by(ma4.72.0, line, col, enviro), mean(data), exper="4"))
names(ma4.0.72.rep)[2] <- "rep"
ma4b.0.72.rep <- data.frame(summarise(group_by(ma4b.72.0, line, col, enviro), mean(data), exper="4b"))
names(ma4b.0.72.rep)[2] <- "rep"
ma5.0.72.rep <- data.frame(summarise(group_by(ma5.72.0, line, rep, enviro), mean(data), exper="5"))
ma6.0.72.rep <- data.frame(summarise(group_by(ma6.72.0, line, rep, enviro), mean(data), exper="6"))
ma7.0.72.rep <- data.frame(summarise(group_by(ma7.72.0, line, rep, enviro), mean(data), exper="7"))
ma47.0.72.rep <- data.frame(summarise(group_by(ma47.72.0, line, rep, enviro), mean(data), exper="47"))

#all0.72 <- data.frame(rbind(ma4.0.72.rep, ma4b.0.72.rep, ma5.0.72.rep, ma6.0.72.rep, ma7.0.72.rep, ma47.0.72.rep))
all0.72 <- data.frame(rbind(ma4.0.72.rep, ma4b.0.72.rep, ma5.0.72.rep, ma6.0.72.rep,  ma47.0.72.rep))
names(all0.72)[4] <- "data"
write.csv(all0.72, "data_in/MIC/aggregated/OD_all0_72.csv",  row.names=FALSE)

all0.rep.72 <- summarise(group_by(all0.72, line, rep, enviro), data = mean(data), sd = sd(data), se = se(data))
write.csv(all0.rep.72, "data_in/MIC/aggregated/OD_all0_72_rep.csv",  row.names=FALSE)

ma4.0.72.ag <- summarise(group_by(ma4.72.0, line, enviro), mean(data), rep="4")
ma4b.0.72.ag <- summarise(group_by(ma4b.72.0, line, enviro), mean(data), rep="4b")
ma5.0.72.ag <- summarise(group_by(ma5.72.0, line, enviro), mean(data), rep="5")
ma6.0.72.ag <- summarise(group_by(ma6.72.0, line, enviro), mean(data), rep="6")
ma7.0.72.ag <- summarise(group_by(ma7.72.0, line, enviro), mean(data), rep="7")
ma47.0.72.ag <- summarise(group_by(ma47.72.0, line, enviro), mean(data), rep="47")


#all0.72ag <- data.frame(rbind(ma4.0.72.ag, ma4b.0.72.ag, ma5.0.72.ag, ma6.0.72.ag, ma7.0.72.ag, ma47.0.72.ag))
all0.72ag <- data.frame(rbind(ma4.0.72.ag, ma4b.0.72.ag, ma5.0.72.ag, ma6.0.72.ag, ma47.0.72.ag))
names(all0.72ag)[3] <- "data"
all0.ag.72 <- summarise(group_by(all0.72ag, line, enviro), data = mean(data), sd = sd(data), se = se(data))
write.csv(all0.ag.72, "data_in/MIC/aggregated/OD_all0_72_ag.csv",  row.names=FALSE)


#####################################################################################
#Compare t10s
#gen = t10 MA4, MA6, MA7
#AGGREGATE AMONG THE FIVE DATASETS
######################################

#######
#24h
#######
ma4.10.rep <- ma4.24 %>%
	filter(gen=="t10") %>%
		group_by(line, col, enviro) %>%
			summarise(dataM = mean(data))
			ma4.10.rep$exper = "4"
names(ma4.10.rep)[2] <-"rep"
ma4.10.rep <- data.frame(ma4.10.rep)

ma4b.10.rep <-ma4b.24 %>%
	filter(gen=="t10") %>%
		group_by(line, col, enviro) %>%
			summarise(dataM = mean(data))
ma4b.10.rep$exper = "4b"
names(ma4b.10.rep)[2] <-"rep"
ma4b.10.rep <- data.frame(ma4b.10.rep)

ma6.10.rep <- ma6.24 %>%
	filter(gen=="t10") %>%
		group_by(line, rep, enviro) %>%
			summarise(dataM = mean(data))
ma6.10.rep$exper = "6"
ma6.10.rep <- data.frame(ma6.10.rep)

ma7.10.rep <- ma7.24 %>%
	filter(gen=="t10") %>%
		group_by(line, rep, enviro) %>%
			summarise(dataM = mean(data))
ma7.10.rep$exper = "7"
ma7.10.rep <- data.frame(ma7.10.rep)

#all10 <- rbind(ma4.10.rep, ma4b.10.rep, ma6.10.rep, ma7.10.rep)
all10 <- rbind(ma4.10.rep, ma4b.10.rep, ma6.10.rep)
write.csv(all10, "data_in/MIC/aggregated/OD_all10_24.csv",  row.names=FALSE)
all10.rep <- summarise(group_by(all10, line, rep, enviro), data = mean(dataM), sd = sd(dataM), se = se(dataM))
write.csv(all10.rep, "data_in/MIC/aggregated/OD_all10_24_rep.csv",  row.names=FALSE)

ma4.10.ag <- ma4.24 %>%
	filter(gen=="t10") %>%
		group_by(line, enviro) %>%
			summarise(dataM = mean(data), sd = sd(data), se=se(data))
ma4.10.ag$exper <- "4"

ma4b.10.ag <- ma4b.24 %>%
	filter(gen=="t10") %>%
		group_by(line, enviro) %>%
			summarise(dataM = mean(data), sd = sd(data), se=se(data))
ma4b.10.ag$exper <- "4b"

ma6.10.ag <- ma6.24 %>%
	filter(gen=="t10") %>%
		group_by(line, enviro) %>%
			summarise(dataM = mean(data), sd = sd(data), se=se(data))
ma6.10.ag$exper <- "6"

ma7.10.ag <- ma7.24 %>%
	filter(gen=="t10") %>%
		group_by(line, enviro) %>%
			summarise(dataM = mean(data), sd = sd(data), se=se(data))
ma7.10.ag$exper <- "7"

#this is aggregated by experiment among lines
#all10ag <- data.frame(rbind(ma4.10.ag, ma4b.10.ag, ma6.10.ag, ma7.10.ag))
all10ag <- data.frame(rbind(ma4.10.ag, ma4b.10.ag, ma6.10.ag))
#this is aggregated completely
all10.ag <- summarise(group_by(all10ag, line, enviro), data = mean(dataM), sd = sd(dataM), se = se(dataM))
write.csv(all10.ag, "data_in/MIC/aggregated/OD_all10_24_ag.csv",  row.names=FALSE)

#######
#48h
#######
names(ma4.48)[3] <- "rep"
ma4.10.rep.48 <- ma4.48 %>%
	filter(gen=="t10") %>%
		group_by(line, rep, enviro) %>%
			summarise(dataM = mean(data))
ma4.10.rep.48$exper = "4"

names(ma4b.48)[3] <- "rep"
ma4b.10.rep.48 <- ma4b.48 %>%
	filter(gen=="t10") %>%
		group_by(line, rep, enviro) %>%
			summarise(dataM = mean(data))
ma4b.10.rep.48$exper = "4b"

ma6.10.rep.48 <- ma6.48 %>%
	filter(gen=="t10") %>%
		group_by(line, rep, enviro) %>%
			summarise(dataM = mean(data))
ma6.10.rep.48$exper = "6"
names(ma6.10.rep.48)[2] <- "rep"

ma7.10.rep.48 <- ma7.48 %>%
	filter(gen=="t10") %>%
		group_by(line, rep, enviro) %>%
			summarise(dataM = mean(data))
ma7.10.rep.48$exper = "7"
names(ma7.10.rep.48)[2] <- "rep"

#all10.48 <- data.frame(rbind(ma4.10.rep.48, ma4b.10.rep.48, ma6.10.rep.48, ma7.10.rep.48))
all10.48 <- data.frame(rbind(ma4.10.rep.48, ma4b.10.rep.48, ma6.10.rep.48))
write.csv(all10.48, "data_in/MIC/aggregated/OD_all10_48.csv",  row.names=FALSE)

all10.rep.48 <- summarise(group_by(all10.48, line, rep, enviro), data = mean(dataM), sd = sd(dataM), se = se(dataM))
write.csv(all10.rep.48, "data_in/MIC/aggregated/OD_all10_48_rep.csv",  row.names=FALSE)

ma4.10.ag.48 <- ma4.48 %>%
	filter(gen=="t10") %>%
		group_by(line, enviro) %>%
			summarise(dataM = mean(data), sd = sd(data), se=se(data), exper="4")
ma4.10.ag.48$exper = "4"
ma4b.10.ag.48 <- ma4b.48 %>%
	filter(gen=="t10") %>%
		group_by(line, enviro) %>%
			summarise(dataM = mean(data), sd = sd(data), se=se(data), exper="4b")
ma4b.10.ag.48$exper = "4b"

ma6.10.ag.48 <- ma6.48 %>%
	filter(gen=="t10") %>%
		group_by(line, enviro) %>%
			summarise(dataM = mean(data), sd = sd(data), se=se(data), exper="6")

ma7.10.ag.48 <- ma7.48 %>%
	filter(gen=="t10") %>%
		group_by(line, enviro) %>%
			summarise(dataM = mean(data), sd = sd(data), se=se(data), exper="7")

#all10.48ag <- data.frame(rbind(ma4.10.ag.48, ma4b.10.ag.48, ma6.10.ag.48, ma7.10.ag.48))
all10.48ag <- data.frame(rbind(ma4.10.ag.48, ma4b.10.ag.48, ma6.10.ag.48))
all10.ag.48 <- summarise(group_by(all10.48ag, line, enviro), data = mean(dataM), sd = sd(dataM), se = se(dataM))
write.csv(all10.ag.48, "data_in/MIC/aggregated/OD_all10_48_ag.csv",  row.names=FALSE)

#######
#72h
#######
names(ma4.72)[3] <- "rep"
ma4.10.rep.72 <- ma4.72 %>%
	filter(gen=="t10") %>%
		group_by(line, rep, enviro) %>%
			summarise(dataM = mean(data), exper="4")

names(ma4b.72)[3] <- "rep"
ma4b.10.rep.72 <- ma4b.72 %>%
	filter(gen=="t10") %>%
		group_by(line, rep, enviro) %>%
			summarise(dataM = mean(data), exper="4b")

ma6.10.rep.72 <- ma6.72 %>%
	filter(gen=="t10") %>%
		group_by(line, rep, enviro) %>%
			summarise(dataM = mean(data), exper="6")

ma7.10.rep.72 <- ma7.72 %>%
	filter(gen=="t10") %>%
		group_by(line, rep, enviro) %>%
			summarise(dataM = mean(data), exper="7")

#all10.72 <- data.frame(rbind(ma4.10.rep.72, ma4b.10.rep.72, ma6.10.rep.72, ma7.10.rep.72))
all10.72 <- data.frame(rbind(ma4.10.rep.72, ma4b.10.rep.72, ma6.10.rep.72))
write.csv(all10.72, "data_in/MIC/aggregated/OD_all10_72.csv",  row.names=FALSE)

all10.rep.72 <- summarise(group_by(all10.72, line, rep, enviro), data = mean(dataM), sd = sd(dataM), se = se(dataM))
write.csv(all10.rep.72, "data_in/MIC/aggregated/OD_all10_72_rep.csv",  row.names=FALSE)

ma4.10.ag.72 <- ma4.72 %>%
	filter(gen=="t10") %>%
		group_by(line, enviro) %>%
			summarise(dataM = mean(data), sd = sd(data), se=se(data), exper="4")
ma4b.10.ag.72 <- ma4b.72 %>%
	filter(gen=="t10") %>%
		group_by(line, enviro) %>%
			summarise(dataM = mean(data), sd = sd(data), se=se(data), exper="4b")
ma6.10.ag.72 <- ma6.72 %>%
	filter(gen=="t10") %>%
		group_by(line, enviro) %>%
			summarise(dataM = mean(data), sd = sd(data), se=se(data), exper="6")
ma7.10.ag.72 <- ma7.72 %>%
	filter(gen=="t10") %>%
		group_by(line, enviro) %>%
			summarise(dataM = mean(data), sd = sd(data), se=se(data), exper="7")

#all10.72ag <- data.frame(rbind(ma4.10.ag.72, ma4b.10.ag.72, ma6.10.ag.72, ma7.10.ag.72))
all10.72ag <- data.frame(rbind(ma4.10.ag.72, ma4b.10.ag.72, ma6.10.ag.72))
all10.ag.72 <- summarise(group_by(all10.72ag, line, enviro), data = mean(dataM), sd = sd(dataM), se = se(dataM))
write.csv(all10.ag.72, "data_in/MIC/aggregated/OD_all10_72_ag.csv",  row.names=FALSE)
