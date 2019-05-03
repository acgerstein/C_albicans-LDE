#################################
#Load data and libraries
#################################
source("scripts/04_combine-MIC-flow.R")

#################################
#Fitness in FLC1
#################################
# The majority of evolved replicates grew to a higher optical density (OD) in 1ug/mL than the ancestral strains, when measured at either 24 h or 72 h (Figure 1, paired t-test results in Table S1).
comp1 <- list()
for(i in 1:20){
  sub0 <- subset(all0.1, line==i)
  sub10 <- subset(all10.1, line==i)
  test <- t.test(sub0$all0.1, sub10$all10.1)
  comp1[[i]] <- c(test$statistic, test$parameter, test$p.value)
}
FLC1.ttest <- do.call("rbind", comp1)
FLC1.ttest <- data.frame(line =1:20, FLC1.ttest)
names(FLC1.ttest)[4] <- "p.value"
FLC1.ttest$sig <- ifelse(FLC1.ttest$p.value < 0.05, "*", "")
FLC1.ttest$t <- round(FLC1.ttest$t, 2)
FLC1.ttest$df <- round(FLC1.ttest$df, 2)
FLC1.ttest$p.value <- round(FLC1.ttest$p.value, 4)
write.csv(FLC1.ttest, "manuscript/tables/TablsS1a-FLC1ttest_24h.csv", row.names=FALSE)

comp1.72 <- list()
for(i in 1:20){
  sub0 <- subset(all0.1.72, line==i)
  sub10 <- subset(all10.1.72, line==i)
  test <- t.test(sub0$ all0.1.72, sub10$all10.1.72)
  comp1.72[[i]] <- c(test$statistic, test$parameter, test$p.value)
}
FLC1.ttest.72 <- do.call("rbind", comp1.72)
FLC1.ttest.72 <- data.frame(line =1:20, FLC1.ttest.72)
names(FLC1.ttest.72)[4] <- "p.value"
FLC1.ttest.72$sig <- ifelse(FLC1.ttest.72$p.value < 0.05, "*", "")
print(FLC1.ttest.72)
FLC1.ttest.72$t <- round(FLC1.ttest.72$t, 2)
FLC1.ttest.72$df <- round(FLC1.ttest.72$df, 2)
FLC1.ttest.72$p.value <- round(FLC1.ttest.72$p.value, 4)
write.csv(FLC1.ttest.72, "manuscript/tables/TableS1b-FLC1ttest_72h.csv", row.names=FALSE)

##########################
#https://paperpile.com/view/adc154f8-62c6-08e2-9d84-eb0bbc7c7222
#The initial ratio was computed as log(wi_final/wmax_start) whereas the ratio of fitness gains was calculated as log(wi_final/wi_start) – log(wmax_final/wmax_start).

wmax_start <- max(fitFlow.ag$all0.1)
wmax_final <- max(fitFlow.ag$all10.1)
init <- 1+log10(fitFlow.ag$all0.1/wmax_start)/100
fitnessGain <- 1+(log10(fitFlow.ag$all10.1/fitFlow.ag$all0.1) - log10(wmax_final/wmax_start))/100

plot(init, fitnessGain)
summary(lm(fitnessGain~init))
#r-squared 0.7191

wmax_start72 <- max(fitFlow.ag$all0.1.72)
wmax_final72 <- max(fitFlow.ag$all10.1.72)
init72 <- 1+log(fitFlow.ag$all0.1.72/wmax_start72)/100
fitnessGain72 <- 1+(log(fitFlow.ag$all10.1.72/fitFlow.ag$all0.1.72) - log(wmax_final72/wmax_start72))/100
summary(lm(fitnessGain72~init72))
#r-squared 0.955
plot(init72, fitnessGain72)

##################
#whole strain set
##################
#The degree of fitness improvement was strongly influenced by initial strain fitness, while neither mating type nor clade had a significant effect (Figure 2, ANOVA tests with change in fitness values transformed using the log-modululs transformation to meet parametric assumptions; 24 h---initial fitness: F1, 12 = 15.55, p = 0.002, clade: F4, 12 = 1.48, p = 0.27, MAT zygosity: F1, 12 = 1.70, p = 0.22; 72 h---initial fitness: F1, 12 = 7.20, p = 0.020, clade: F4, 12 = 12.37, p = 0.16, MAT zygosity: F1, 11 = 0.04, p = 0.84).
#https://blogs.sas.com/content/iml/2014/07/14/log-transformation-of-pos-neg.html

delta1.aov <- aov(delta1.24~ all0.1 +zygosity + as.factor(clade) , data=fitFlow.ag)
# Df Sum Sq Mean Sq F value Pr(>F)
# all0.ag.1         1 0.5410  0.5410   7.994 0.0143 *
# zygosity          1 0.0184  0.0184   0.271 0.6112
# as.factor(clade)  4 0.1734  0.0433   0.640 0.6430
# Residuals        13 0.8798  0.0677


delta1.72.aov <- aov(delta1.72~ all0.ag.1.72 +zygosity + clade , data=fitFlow.ag)
# Df Sum Sq Mean Sq F value       Pr(>F)
# all0.ag.1.72  1 1.3723  1.3723 158.730 0.0000000116 ***
# zygosity      1 0.0021  0.0021   0.248        0.626
# clade         4 0.0681  0.0170   1.970        0.159
# Residuals    13 0.1124  0.0086

delta1.variance24.aov <- aov(c(fitFlow.sd$all10.1-fitFlow.sd$all0.1) ~ fitFlow.ag$all0.1 + fitFlow.sd$zygosity + as.factor(fitFlow.sd$clade))
# Df  Sum Sq Mean Sq F value    Pr(>F)
# fitFlow.ag$all0.1            1 0.24434 0.24434  40.941 0.0000235 ***
#   fitFlow.sd$zygosity          1 0.01325 0.01325   2.220     0.160
# as.factor(fitFlow.sd$clade)  4 0.01950 0.00488   0.817     0.537
# Residuals                   13 0.07759 0.00597

delta1.variance72.aov <- aov(c(fitFlow.sd$all10.1.72-fitFlow.sd$all0.1.72) ~ fitFlow.ag$all0.1.72 +zygosity + clade , data=fitFlow.ag)
# Df  Sum Sq Mean Sq F value Pr(>F)
# fitFlow.ag$all0.1.72  1 0.03437 0.03437   6.352 0.0256 *
#   zygosity              1 0.00614 0.00614   1.136 0.3060
# clade                 4 0.02195 0.00549   1.014 0.4356
# Residuals            13 0.07033 0.00541


# There is also a significant negative correlation between the degree of fitness change in low drug and the variance in evolved fitness among replicates – replicate lines were more variable in strain backgrounds that increased in fitness the most.
cor.test(fitFlow.ag$delta1.24, fitFlow.sd$all10.ag.1, paired=TRUE, method="spearman") #S = 694, p-value = 0.03451, rho = 0.48
cor.test(fitFlow.ag$delta1.72, fitFlow.sd$all10.ag.1.72, paired=TRUE, method="spearman") #S = 476, p-value = 0.00286, rho = 0.64

#################################
#Fitness in YPD
#################################
# When fitness was assessed in a permissive environment (standard lab YPD, the base medium for the evolutionary environment), there was a significant reduction in fitness for evolved strains compared to ancestral most strain backgrounds at 24 h

comp0 <- list()
stat <- c()
for(i in 1:20){
  sub0 <- subset(all0.0, line==i)
  s0 <- shapiro.test(sub0$data)
  sub10 <- subset(all10.0, line==i)
  s10 <- shapiro.test(sub10$data)
  if(s0$p.value > 0.05 & s10$p.value > 0.05){
      v <- var.test(sub0$data, sub10$data)
      if(v$p.value > 0.05){
         test <- t.test(sub0$data, sub10$data, var.equal=TRUE)
         stat <- append(stat, "t-test")
       }
      else{
         test <- t.test(sub0$data, sub10$data, var.equal=FALSE)
         stat <- append(stat, "t-test-welch")
       }
    }
  else{
    test <- wilcox.test(sub0$data, sub10$data)
    test$parameter <- NA
    stat <- append(stat, "wilcox")
  }
  comp0[[i]] <- c(statistic = test$statistic, df = test$parameter, "p.value"= test$p.value, "imp" = ifelse(mean(sub10$data) > mean(sub0$data), "evol", "anc"))
}
FLC0.ttest <- do.call("rbind", comp0)
FLC0.ttest <- data.frame(line =1:20, "test" = stat, FLC0.ttest)
names(FLC0.ttest) <- c("line", "test", "statistic", "df", "p.value", "imp")
FLC0.ttest$sig <- ifelse(FLC0.ttest$p.value < 0.05, "*", "")
FLC0.ttest$statistic<- round(as.numeric(FLC0.ttest$statistic), 2)
FLC0.ttest$df <- round(as.numeric(FLC0.ttest$df), 2)
FLC0.ttest$p.value <- round(as.numeric(FLC0.ttest$p.value), 4)
write.csv(FLC0.ttest, "/Users/acgerstein/Documents/Postdoc/Papers/MutAccum/tables/170802FLC0ttest.csv", row.names=FALSE)


comp0.72 <- list()
stat <- c()
for(i in 1:20){
  sub0 <- subset(all0.0.72, line==i)
  s0 <- shapiro.test(sub0$data)
  sub10 <- subset(all10.0.72, line==i)
  s10 <- shapiro.test(sub10$data)
  if(s0$p.value > 0.05 & s10$p.value > 0.05){
      v <- var.test(sub0$data, sub10$data)
      if(v$p.value > 0.05){
         test <- t.test(sub0$data, sub10$data, var.equal=TRUE)
         stat <- append(stat, "t-test")
       }
      else{
         test <- t.test(sub0$data, sub10$data, var.equal=FALSE)
         stat <- append(stat, "t-test-welch")
       }
    }
  else{
    test <- wilcox.test(sub0$data, sub10$data)
    test$parameter <- NA
    stat <- append(stat, "wilcox")
  }
  comp0.72[[i]] <- c(statistic = test$statistic, df = test$parameter, "p.value"= test$p.value, "imp" = ifelse(mean(sub10$data) > mean(sub0$data), "evol", "anc"))
}
FLC0.ttest.72 <- do.call("rbind", comp0.72)
FLC0.ttest.72 <- data.frame(line =1:20, "test" = stat, FLC0.ttest.72)
names(FLC0.ttest.72) <- c("line", "test", "statistic", "df", "p.value", "imp")
FLC0.ttest.72$sig <- ifelse(FLC0.ttest.72$p.value < 0.05, "*", "")
FLC0.ttest.72$statistic<- round(as.numeric(FLC0.ttest.72$statistic), 2)
FLC0.ttest.72$df <- round(as.numeric(FLC0.ttest.72$df), 2)
FLC0.ttest.72$p.value <- round(as.numeric(FLC0.ttest.72$p.value), 4)
write.csv(FLC0.ttest.72, "/Users/acgerstein/Documents/Postdoc/Papers/MutAccum/tables/170802FLC0ttest-72h.csv", row.names=FALSE)

test0 <-data.frame(line =1:20, YPD=FLC0.ttest$sig, YPD.imp = FLC0.ttest$imp, YPD.72=FLC0.ttest.72$sig, YPD.imp.72 = FLC0.ttest.72$imp)
test0

cor.test(fitFlow.ag$delta1.24, fitFlow.sd$all10.ag.1, paired=TRUE, method="spearman") #S = 694, p-value = 0.03451, rho = 0.48

# For our low drug evolved strains, however, we recovered a significant positive relationship among strain backgrounds when we compare mean fitness assessed at 24 h at low drug and the permissive environment (spearman’s rank correlation test; S = 396, p = 0.0008, rho = 0.702; there is not enough variation in the permissive environment at 72 h for correlative purposes)
cor.test(all10.ag.0$data, fitFlow.ag$all10.ag.1, method="spearman")

parameter <- c()
p.value <- c()
rho <- c()
for(i in 1:20){
  test <- cor.test(subset(all10.0, line==i)$data, subset(all10.1, line==i)$data, method="spearman")
  parameter[i] <- test$statistic
  p.value[i] <- test$p.value
  rho[i] <- test$estimate
}

ypdCor <- data.frame(line = 1:20, s =parameter, p.value, rho)
write.csv(ypdCor, "/Users/acgerstein/Documents/Postdoc/Papers/MutAccum/tables/171030YPD-FLC1_cor-24h.csv", row.names=FALSE)

#################################
#Higher drug
#################################
# After ~100 generations of evolution at a low level of fluconazole, only 18 replicates from seven strain backgrounds increased in MIC beyond 1ug/mL measured through a broth microdilution assay (BMD). In most of these replicates, MIC increased only slightly, to 4 ug/mL (five replicates from strain A3, four replicates from A8, and one replicate from A4 and A7). Clinical resistance (MIC >= 8 at 24 hours, Fothergill 2014) was acquired in only eight replicates, six from strain A5 (which had an initial MIC of 4), and one replicate each from strains A17 and A19, two backgrounds with a very low initial MIC (Figure 1A).

MICto1 <- subset(fitFlow, MIC24 < 1 & MIC24.10 == 1) #97 of 132 reps increased to 1, 12 increased above 1
MICinc1 <- subset(fitFlow, MIC24 <= 1 & MIC24.10 > 1)
MICinc8 <- subset(fitFlow, MIC24 <= 8 & MIC24.10 >= 8)

nrow(subset(fitFlow, MIC24 == 1 & MIC24.10 < 1)) #4
nrow(subset(fitFlow, MIC24 > 1 & MIC24.10 == 1)) #40

 # Interestingly, the MIC of some replicates from initially-high MIC strains actually decreased to 1ug/mL  (two replicates from A5, eight replicates from A12, all 12 replicates from A20).

MICdec <- subset(fitFlow, change24 < 0)

#SMG
SMG24up <- subset(fitFlow, SMG24.10.2 > rep(SMG24.up, each=12)) #72 (13 strains)
SMG24down <- subset(fitFlow, SMG24.10.2 < rep(SMG24.down, each=12)) #79 (17 strains)

table(SMG24up$strain)
# 2  3  4  6  8  9 11 12 13 14 15 17 19
# 3 12  2 11  4  5  7  8  3  3  1  7  6
table(SMG24down$strain)
# 1  4  5  7  8  9 10 11 12 13 14 15 16 17 18 19 20
# 6  2  1  2  2  3  7  3  1  3  5 10 10  2 12  2  8

SMG72up <- subset(fitFlow, SMG72.10.2 > rep(SMG72.up, each=12)) #73 (13 strains)
SMG72down <- subset(fitFlow, SMG72.10.2 < rep(SMG72.down, each=12)) #91 (17 strains)

table(SMG72up$strain)
# 1  3  4  5  6  7  8  9 10 11 12 13 15 17 19 20
# 3 12  3  1  2  4  2  6  2 11  8  3  2  2  7  5
table(SMG72down$strain)
# 1  2  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
# 1  8  8  1  4  1  5  4  6  1  1  6  9  8  9  7  6  2  4



# Mean strain improvements in resistance and trailing growth were not influenced by ancestral strain behaviour, mating type, nor clade (ANOVA tests; change in MIC---initial MIC: F1, 12 = 2.32, p = 0.15; clade = F3, 12 = 0.54, p = 0.76; zygosity: F1, 12 = 0.035, p = 0.854; change in trailing growth---initial trailing: F1, 12 = 2.94, p = 0.11; clade = F3, 12 = 0.38, p = 0.80; zygosity: F1, 12 =0.75, p = 0.4; change in RAD20---initial RAD: F1, 12 = 0.042, p = 0.84; clade = F3, 12 = 2.98, p = 0.07; zygosity: F1, 12 =0.48, p = 0.5). By contrast, improvements in FoG were negatively correlated with initial FoG (but were not influenced by clade or zygosity), similar to the pattern seen with initial fitness in the low drug environment (change in FoG---initial FoG: F1, 12 = 6.45, p = 0.03; clade = F3, 12 = 0.53, p = 0.67; zygosity: F1, 12 = 1.99, p = 0.18).

change24.SMG.aov <- lme(change24.SMG.2 ~ all0.ag.1+zygosity+as.factor(clade), random =list(strain =~1), data=fitFlow, na.action=na.omit)
anova(change24.SMG.aov)
#   numDF denDF   F-value p-value
# (Intercept)          1   217 0.3145514  0.5755
# all0.ag.1            1   217 0.5146739  0.4739
# zygosity             1    14 2.5293538  0.1341
# as.factor(clade)     4    14 1.7892436  0.1871

change72.SMG.aov <- lme(change72.SMG.2~ all0.ag.1+zygosity + as.factor(clade) ,random =list(strain =~1), data=fitFlow, na.action=na.omit)
anova(change72.SMG.aov)
# numDF denDF   F-value p-value
# (Intercept)          1   217 0.3145514  0.5755
# all0.ag.1            1   217 0.5146739  0.4739
# zygosity             1    14 2.5293538  0.1341
# as.factor(clade)     4    14 1.7892436  0.1871

cor.test(fitFlow.sd$SMG24.10.2, fitFlow.sd$all10.ag.1, method="spearman") #S = 960, p-value = 0.2341

cor.test(fitFlow.sd$SMG72.10.2, fitFlow.sd$all10.ag.1, method="spearman") #S = 588, p-value = 0.01183,  cor = 0.56

cor.test(fitFlow.sd$SMG72.10.2, fitFlow.sd$all10.ag.1, method="spearman") #S = 588, p-value = 0.01183, rho = 0.56
cor.test(fitFlow.sd$SMG72.10.2, fitFlow.sd$all10.ag.1.72, method="spearman") #S = 762, p-value = 0.06173, rho = 0.43

#comment here when capacity
shapiro.test(fitFlow.ag.variable$all0.ag.1) #W = 0.94881, p-value = 0.4711
shapiro.test(fitFlow.ag.variable$all0.ag.1.72) #W = 0.93587, p-value = 0.3015
shapiro.test(fitFlow.ag.variable$t10.G1.1) #W = 0.9694, p-value = 0.8288
shapiro.test(fitFlow.ag.variable$change24.SMG.2) #W = 0.86816, p-value = 0.02549
shapiro.test(fitFlow.ag.variable$change72.SMG.2) #W = 0.84432, p-value = 0.01125

qqnorm(fitFlow.ag.variable$change24)
qqnorm(fitFlow.ag.variable$t10.G1.1)

#is there really a point at looking at the aggreagate among strain data?
fit0_flow.variable <- cor.test(fitFlow.ag.variable$all0.ag.1, fitFlow.ag.variable$t10.G1.1) #t = -0.33692, df = 14, p-value = 0.7412
fit10.72_flow.variable <- cor.test(fitFlow.ag.variable$all10.ag.1.72, fitFlow.ag.variable$t10.G1.1) #t = -0.32112, df = 14, p-value = 0.7529

fitD_flow.variable24 <- cor.test(fitFlow.ag.variable$change24.SMG.2, fitFlow.ag.variable$t10.G1.1, method="spearman") #S = 778.07, p-value = 0.5941
fitD_flow.variable72 <- cor.test(fitFlow.ag.variable$change72.SMG.2, fitFlow.ag.variable$t10.G1.1, method="spearman") #S = 1027.3, p-value = 0.04325, rho = -0.51
plot(fitFlow.ag.variable$change72.SMG.2, fitFlow.ag.variable$t10.G1.1)


# Interestingly, there was a significant correlation between the variance of evolved low drug fitness (measured at 24 h) and variance of evolved genome size (Pearson's correlation, t14 = 5.65, p-value = 0.00002, cor = 0.80, this remains significant when the four non-variable lines are removed: t14 = 2.31, p-value = 0.037, cor = 0.51).

shapiro.test(fitFlow.ag$t10.G1.1) #W = 0.95807, p-value = 0.5061
qqnorm(fitFlow.ag$t10.G1.1)
shapiro.test(fitFlow.sd$all10.ag.1) #W = 0.93818, p-value = 0.2214
shapiro.test(fitFlow.sd.variable$all10.ag.1) #W = 0.97123, p-value = 0.8581
shapiro.test(fitFlow.sd.sub$FoG20.10)
qqnorm(fitFlow.sd.sub$FoG20.10)

fit10_flowCV <- cor.test(fitFlow.sd$all10.ag.1, fitFlow.sd$t10.G1.1, method="spearman") #t = 5.6545, df = 18, p-value = 0.00002304
fit10_flowCV.variable <- cor.test(fitFlow.sd.variable$all10.ag.1, fitFlow.sd.variable$t10.G1.1) #t = 2.3115, df = 14, p-value = 0.03655

flowCV_SMG24 <- cor.test(fitFlow.sd.variable$t10.G1.1, fitFlow.sd.variable$SMG24.10.2, method="spearman")
fit10_SMG24 <- cor.test(fitFlow.sd.variable$all10.ag.1, fitFlow.sd.variable$SMG24.10.2, method="spearman")

flowCV_SMG72 <- cor.test(fitFlow.sd.variable$t10.G1.1, fitFlow.sd.variable$SMG72.10.2, method="spearman")
fit10_SMG72<- cor.test(fitFlow.sd.variable$all10.ag.1, fitFlow.sd.variable$SMG72.10.2, method="spearman")

 plot(fitFlow.sd.variable$all10.ag.1, fitFlow.sd.variable$SMG72.10)
