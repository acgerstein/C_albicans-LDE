cv <-function(x, ...) sd(x,na.rm=TRUE)/mean(x,na.rm=TRUE)

#ANC & LDE
flow <- read.csv("data_in/flow/MA4_flow.csv", sep=",", as.is=TRUE)
names(flow)[1] <- "name"

#factor out line and well information
flow$line <- unlist(lapply(as.character(flow$name), function(x) strsplit(x, "-")[[1]][1]))
flow$well <- unlist(lapply(as.character(flow$name), function(x) strsplit(x, "-")[[1]][2]))

flow$line[flow$line=="SC"] <- "A17"
flow$line[flow$line=="FH"] <- "A18"
flow$line[flow$line=="DS"] <- "A19"
flow$line[flow$line=="T1"] <- "A20"

flow <- flow[order(flow$line),]
flow$line <- as.numeric(substring(flow$line, 2, 3))

flow.ag <- aggregate(flow[c("t0.G1.1", "t3.G1.1", "t5.G1.1", "t8.G1.1", "t10.G1.1", "t10b.G1.1")], flow["line"], median, na.rm=TRUE)
flow.ag.sd <- aggregate(flow[c("t0.G1.1", "t3.G1.1", "t5.G1.1", "t8.G1.1", "t10.G1.1", "t10b.G1.1")], flow["line"], sd, na.rm=TRUE)
flow.ag.cv <- aggregate(flow[c("t0.G1.1", "t3.G1.1", "t5.G1.1", "t8.G1.1", "t10.G1.1", "t10b.G1.1")], flow["line"], cv, na.rm=TRUE)
