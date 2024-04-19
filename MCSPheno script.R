library(haven)
library(data.table)
library(mice)

pheno1 <-read_sav('C:/Users/lucas/Documents/documents_20230428/GDAC_2022_18_LUCIANO_mcs_basic_demographics_v0003_shareable_20220215.sav')
pheno2 <-read_sav('C:/Users/lucas/Documents/documents_20230428/GDAC_2022_18_LUCIANO_mcs_cm_structure_pheno_data_2023-04-26_17-21-10.sav', encoding="latin1")
pheno3 <-read_sav('C:/Users/lucas/Documents/documents_20230428/GDAC_2022_18_LUCIANO_mcs_family_structure_pheno_data_2023-04-26_17-21-10.sav')
pheno4 <-read_sav('C:/Users/lucas/Documents/documents_20230428/GDAC_2022_18_LUCIANO_mcs_hhgrid_structure_pheno_data_2023-04-26_17-21-10.sav')
pheno5 <-read_sav('C:/Users/lucas/Documents/documents_20230428/GDAC_2022_18_LUCIANO_mcs_parent_cm_structure_pheno_data_2023-04-26_17-21-10.sav')
pheno6 <-read_sav('C:/Users/lucas/Documents/documents_20230428/GDAC_2022_18_LUCIANO_mcs_parent_structure_pheno_data_2023-04-26_17-21-10.sav', encoding="latin1")
pheno7 <-read_sav('C:/Users/lucas/Documents/documents_20230428/GDAC_2022_18_LUCIANO_mcs_data_struct_fam_2023-08-02_16-46-43.sav')
pheno8 <-read_sav('C:/Users/lucas/Documents/documents_20230428/GDAC_2022_18_LUCIANO_3_mcs_cm_structure_pheno_data_2023-07-10_14-52-34.sav')
pheno9 <-read_sav('C:/Users/lucas/Documents/documents_20230428/GDAC_2022_18_LUCIANO_4_mcs_data_struct_cm_long_2023-12-18_15-56-15.sav', encoding="latin1")
pheno10 <-read_sav('C:/Users/lucas/Documents/documents_20230428/GDAC_2022_18_LUCIANO_4_mcs_data_struct_par_long_2023-12-18_15-56-15.sav', encoding="latin1")

CogModelW5 <-read.csv('C:/Users/lucas/Documents/documents_20230428/mcs5_lucianoIDs.csv')
CogModelW6 <-read.csv('C:/Users/lucas/Documents/documents_20230428/mcs6_lucianoIDs.csv')
CogModelW5 <- CogModelW5[,-3:-5]
CogModelW6 <- CogModelW6[,-3:-5]
colnames(CogModelW5) <- toupper(colnames(CogModelW5))
colnames(CogModelW6) <- toupper(colnames(CogModelW6))

pheno3n <- cbind(pheno3[,1:2] ,(sapply(pheno3[,3:50],as.numeric)))
pheno3ns <- pheno3n[,-8]
mice_imputed_P1 <- data.frame(
  original = pheno3ns$ADHTYS00,
  imputed_cart = complete(mice(pheno3ns, method = "cart"))$ADHTYS00
)
mice_imputed_P2 <- data.frame(
  original = pheno3n$BDHTYS00,
  imputed_cart = complete(mice(pheno3n, method = "cart"))$BDHTYS00
)
mice_imputed_P3 <- data.frame(
  original = pheno3n$CDHTYS00,
  imputed_cart = complete(mice(pheno3n, method = "cart"))$CDHTYS00
)
mice_imputed_P4 <- data.frame(
  original = pheno3n$DDHTYS00,
  imputed_cart = complete(mice(pheno3n, method = "cart"))$DDHTYS00
)
mice_imputed_P5 <- data.frame(
  original = pheno3n$EHTYS00,
  imputed_cart = complete(mice(pheno3n, method = "cart"))$EHTYS00
)
mice_imputed_P6 <- data.frame(
  original = pheno3n$FDHTYS00,
  imputed_cart = complete(mice(pheno3n, method = "cart"))$FDHTYS00
)

load("~/ImputedData/mice_imputed_P1.RData")
load("~/ImputedData/mice_imputed_P2.RData")
load("~/ImputedData/mice_imputed_P3.RData")
load("~/ImputedData/mice_imputed_P4.RData")
load("~/ImputedData/mice_imputed_P5.RData")
load("~/ImputedData/mice_imputed_P6.RData")

pheno3sub <- as.data.frame(cbind(pheno3[,1],pheno3[,9],pheno3[,4],pheno3[,16],pheno3[,23],pheno3[,30],pheno3[,41]))
pheno3sub[,2:7] <- sapply(pheno3sub[,2:7],as.numeric)
# pheno3sub[,2:7] <- sapply(pheno3sub[,2:7],function(x){x*(-1) + 3})

pheno3sub$ADHTYS00 <-mice_imputed_P1$imputed_cart
pheno3sub$BDHTYS00 <-mice_imputed_P2$imputed_cart
pheno3sub$CDHTYS00 <-mice_imputed_P3$imputed_cart
pheno3sub$DDHTYS00 <-mice_imputed_P4$imputed_cart
pheno3sub$EHTYS00 <-mice_imputed_P5$imputed_cart
pheno3sub$FDHTYS00 <-mice_imputed_P6$imputed_cart

pheno3sub$SingleParentScore <- rowSums(pheno3sub[,2:7], na.rm = TRUE)

CANTAB <- cbind(pheno2[,1:2],pheno2[,1538:1648])
CANTABsub <- CANTAB[38:44]
CANTABsub[is.na(CANTABsub)] <- 0
CANTAB[38:44] <- CANTABsub
CANTAB$Distraction<-paste(rowSums(CANTAB[38:42]))
CANTAB<-transform(CANTAB, Distraction= as.numeric(Distraction))
CANTAB$Interruption<-paste(rowSums(CANTAB[43:44]))
CANTAB<-transform(CANTAB, Interruption= as.numeric(Interruption))

CANTAB2 <- cbind(pheno2[,1:2],pheno2[,1645],pheno2[,1686:1694],pheno2[,1718])
CANTAB2sub <- CANTAB2[5:11]
CANTAB2sub[is.na(CANTAB2sub)] <- 0
CANTAB2[5:11] <- CANTAB2sub
CANTAB2$Distraction<-paste(rowSums(CANTAB2[5:9]))
CANTAB2<-transform(CANTAB2, Distraction= as.numeric(Distraction))
CANTAB2$Interruption<-paste(rowSums(CANTAB2[10:11]))
CANTAB2<-transform(CANTAB2, Interruption= as.numeric(Interruption))
CANTAB2[,16] <- (pheno2[,1638])

pheno7$IMD_Decile <- pheno7$FIMDSCOE
pheno7$IMD_Decile[!is.na(pheno7$FIMDSCON)] = pheno7$FIMDSCON[!is.na(pheno7$FIMDSCON)]
pheno7$IMD_Decile[!is.na(pheno7$FISIMDSC)] = pheno7$FISIMDSC[!is.na(pheno7$FISIMDSC)]
pheno7$IMD_Decile[!is.na(pheno7$FIWIMDSC)] = pheno7$FIWIMDSC[!is.na(pheno7$FIWIMDSC)]

pheno1sub <- pheno1[,1:4]
colnames(pheno1sub) <- toupper(colnames(pheno1sub))
pheno10sub <- pheno10
pheno10sub <- merge(pheno1sub, pheno10sub)
pheno10sub <- pheno10sub[pheno10sub$MFC == 'M',]
pheno10sub2 <- pheno10sub
pheno10sub2[is.na(pheno10sub2)] <- 0
pheno10sub2$APCIPR00[pheno10sub2$APCIPR00 > 0] <- 1
pheno10sub2$APCICH00[pheno10sub2$APCICH00 > 0] <- 1
pheno10sub3 <- cbind(pheno10sub2$APCIPR00,pheno10sub2$APCICH00)
pheno10sub2$AnySmoking<-paste(as.numeric(rowSums(pheno10sub3[,1:2])))
pheno10sub2$AnySmoking[pheno10sub2$AnySmoking > 1] <- 1

pheno10sub4 <- pheno10sub[pheno10sub$MFC == 'M',]
pheno10sub4$AnySmoking <- pheno10sub2$AnySmoking
pheno10sub5 <- pheno10sub4[pheno10sub4$AnySmoking == '1',]
pheno10sub5$APSMCH00[is.na(pheno10sub5$APSMCH00)] <- 1
pheno10sub6 <- pheno10sub5[pheno10sub5$APSMCH00 %in% c("2", "3"),]
pheno10sub62 <- pheno10sub5[pheno10sub5$APSMCH00 == 1,]
pheno10sub62$APWHCH00[is.na(pheno10sub62$APWHCH00)] <- 10
pheno10sub62$APCICH00[is.na(pheno10sub62$APCICH00)] <- 0
pheno10sub62$APCIPR00[is.na(pheno10sub62$APCIPR00)] <- 0

pheno10sub6$SMOKINGSEVERITY <- 0
pheno10sub6$SMOKINGSEVERITY[pheno10sub6$APCIPR00 < 11] <- 5
pheno10sub6$SMOKINGSEVERITY[pheno10sub6$APCIPR00 > 19] <- 7
pheno10sub6$SMOKINGSEVERITY[pheno10sub6$SMOKINGSEVERITY == 0] <- 6

pheno10sub7 <- pheno10sub62[pheno10sub62$APWHCH00 %in% c("1", "2", "3"),]
pheno10sub8 <- pheno10sub7[pheno10sub7$APCICH00 == 0,]
pheno10sub8$SMOKINGSEVERITY <- 0
pheno10sub8$SMOKINGSEVERITY[pheno10sub8$APCIPR00 < 11] <- 2
pheno10sub8$SMOKINGSEVERITY[pheno10sub8$APCIPR00 > 19] <- 4
pheno10sub8$SMOKINGSEVERITY[pheno10sub8$SMOKINGSEVERITY == 0] <- 3

pheno10sub9 <- pheno10sub7[pheno10sub7$APCICH00 > 0,]
pheno10sub9$APCICH00[pheno10sub9$APCICH00 > 95] <- 1
pheno10sub9$APCIPR00[is.na(pheno10sub9$APCIPR00)] <- 0

pheno10sub9$SMOKINGSEVERITY <- 0
i <- 1
for (i in  1:929) {
if (pheno10sub9$APCIPR00[i] > pheno10sub9$APCICH00[i]) {pheno10sub9$SMOKINGSEVERITY[i] <- pheno10sub9$APCIPR00[i]}
else {pheno10sub9$SMOKINGSEVERITY[i] <- pheno10sub9$APCICH00[i]}
}

pheno10sub9$SMOKINGSEVERITY[pheno10sub9$SMOKINGSEVERITY < 11] <- 2
pheno10sub9$SMOKINGSEVERITY[pheno10sub9$SMOKINGSEVERITY > 19] <- 4
pheno10sub9$SMOKINGSEVERITY[pheno10sub9$SMOKINGSEVERITY > 10] <- 3

pheno10sub10 <- pheno10sub62[pheno10sub62$APWHCH00 > 3,]
pheno10sub10$APCICH00[pheno10sub10$APCICH00 > 95] <- 1
pheno10sub10$APCICH00[is.na(pheno10sub10$APCICH00)] <- 0
pheno10sub10$APCIPR00[is.na(pheno10sub10$APCIPR00)] <- 0

pheno10sub10$SMOKINGSEVERITY <- 0
i <- 1
for (i in  1:215) {
  if (pheno10sub10$APCIPR00[i] > pheno10sub10$APCICH00[i]) {pheno10sub10$SMOKINGSEVERITY[i] <- pheno10sub10$APCIPR00[i]}
  else {pheno10sub10$SMOKINGSEVERITY[i] <- pheno10sub10$APCICH00[i]}
}

pheno10sub10$SMOKINGSEVERITY[pheno10sub10$SMOKINGSEVERITY < 11] <- 5
pheno10sub10$SMOKINGSEVERITY[pheno10sub10$SMOKINGSEVERITY > 19] <- 7
pheno10sub10$SMOKINGSEVERITY[pheno10sub10$SMOKINGSEVERITY > 10] <- 6

pheno10sub11 <- pheno10sub4[pheno10sub4$AnySmoking == '0',]
pheno10sub11$SMOKINGSEVERITY <- 1

pheno10full <- rbind(pheno10sub6, pheno10sub8, pheno10sub9, pheno10sub10, pheno10sub11)

pheno9sub <- pheno9[,-12:-26]
pheno9sub <- na.omit(pheno9sub)
pheno9sub <- pheno9sub[!is.na(pheno9sub$BCENVI00),]

pc <- prcomp(pheno9sub[,5:20])
print(pc)
summary(pc)

phenoHOME <- pheno9sub[,5:20]
pc <- prcomp(phenoHOME)
print(pc)
summary(pc)

library(corrplot)
library(psych)
library(ggplot2)
library(EFAtools)


datamatrix <- cor(phenoHOME)
corrplot(datamatrix, method="number")
KMO(cor(phenoHOME))


parallel <- fa.parallel(phenoHOME)

fa.none <- fa(r=phenoHOME, 
              nfactors = 3, 
              # covar = FALSE, SMC = TRUE,
              fm="pa", # type of factor analysis we want to use ("pa" is principal axis factoring)
              max.iter=100, # (50 is the default, but we have changed it to 100
              rotate="none") # none rotation
print(fa.none)
fa.diagram(fa.none)

pc <- pca(phenoHOME, nfactors = 3, rotate = "varimax", scores = TRUE)
pc
summary(pc)

pheno9full <- cbind(pheno9sub, pc[["scores"]])

pheno7e <- pheno7[,-2:-40]
pheno8e <- pheno8[,-3:-17]
pheno8e <- pheno8e[,-4:-9]
pheno8e <- pheno8e[!is.na(pheno8e$SWMTOTER),]
pheno9e <- pheno9full[,-3:-20]
pheno10e <- pheno10full[,-3:-8]
pheno11e <- as.data.frame(cbind(pheno3sub$LUCIANO2_FID,pheno3sub$SingleParentScore,pheno3sub$ADHTYS00))
colnames(pheno11e) <- c("LUCIANO2_FID","SingleParentScore","Wave1SingleParent")

CANTABsub <- CANTAB[,-3:-4]
CANTABsub <- CANTABsub[,-5:-42]
CANTABsub <- CANTABsub[,-6:-60]
CANTABsub$ECTECH0B[is.na(CANTABsub$ECTECH0B)] <- 0
CANTABsub$EITIRC00[is.na(CANTABsub$EITIRC00)] <- 3
CANTABsub$EMCS5AGE[is.na(CANTABsub$EMCS5AGE)] <- mean(CANTABsub$EMCS5AGE, na.rm = TRUE)

pheno8e <- merge(pheno8e, CANTABsub)
pheno8e$LUCIANO2_SID[pheno8e$LUCIANO2_SID == ""] <- NA
pheno8e <- pheno8e[!is.na(pheno8e$LUCIANO2_SID),]
pheno8e <- pheno8e[!is.na(pheno8e$SWMTOTER),]
Q <- quantile(pheno8e$SWMTOTER, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(pheno8e$SWMTOTER)
pheno8e<- subset(pheno8e, pheno8e$SWMTOTER > (Q[1] - 1.5*iqr) & pheno8e$SWMTOTER < (Q[2]+1.5*iqr))

COG.lm <- lm(SWMTOTER ~ EMCS5AGE + FCCSEX00 + Distraction + Interruption + ECTECH0B + EITIRC00, data=pheno8e)
anova(COG.lm)
summary(COG.lm)
pheno8e$SWMTotErrAdj <- resid(COG.lm)

Q <- quantile(CogModelW5$BETA, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(CogModelW5$BETA)
CogModelW5<- subset(CogModelW5, CogModelW5$BETA > (Q[1] - 1.5*iqr) & CogModelW5$BETA < (Q[2]+1.5*iqr))

CogModelW5 <- merge(CogModelW5, CANTABsub)

COG2.lm <- lm(BETA ~ EMCS5AGE + FCCSEX00 + Distraction + Interruption + ECTECH0B + EITIRC00, data=CogModelW5)
anova(COG2.lm)
summary(COG2.lm)
CogModelW5$BETAdjst <- resid(COG2.lm)

merge1 <- merge(pheno7e, pheno8e, all=TRUE)
merge2 <- merge(pheno9e, pheno8e, all=TRUE)
merge3 <- merge(pheno10e, pheno8e, by="LUCIANO2_FID", all=TRUE)
merge4 <- merge(pheno8e, CogModelW5)
merge5 <- merge(pheno8e, CogModelW6)
merge6 <- merge(pheno11e, pheno8e, by="LUCIANO2_FID", all=TRUE)

merge1$LUCIANO2_SID[merge1$LUCIANO2_SID == ""] <- NA
merge1 <- merge1[!is.na(merge1$LUCIANO2_SID),]
merge1 <- merge1[!is.na(merge1$IMD_Decile),]
merge2$LUCIANO2_SID[merge2$LUCIANO2_SID == ""] <- NA
merge2 <- merge2[!is.na(merge2$LUCIANO2_SID),]
merge3 <- merge3[,-2]
names(merge3)[names(merge3) == 'LUCIANO2_SID.y'] <- 'LUCIANO2_SID'
merge3$LUCIANO2_SID[merge3$LUCIANO2_SID == ""] <- NA
merge3 <- merge3[!is.na(merge3$LUCIANO2_SID),]
merge3$AnySmoking <- as.numeric(merge3$AnySmoking)
merge3$SWMTotErrAdj <- as.numeric(merge3$SWMTotErrAdj)
merge6$LUCIANO2_SID[merge6$LUCIANO2_SID == ""] <- NA
merge6 <- merge6[!is.na(merge6$LUCIANO2_SID),]
merge6$SingleParentScore <- as.numeric(merge6$SingleParentScore)
merge6$Wave1SingleParent <- as.numeric(merge6$Wave1SingleParent)

cor.test(merge1$IMD_Decile, merge1$SWMTotErrAdj)
cor.test(merge2$RC1, merge2$SWMTotErrAdj)
cor.test(merge2$RC2, merge2$SWMTotErrAdj)
cor.test(merge2$RC3, merge2$SWMTotErrAdj)
cor.test(merge3$AnySmoking, merge3$SWMTotErrAdj)
cor.test(merge3$SMOKINGSEVERITY, merge3$SWMTotErrAdj)
cor.test(merge4$BETAdjst, merge4$SWMTotErrAdj)
cor.test(merge5$BETA, merge5$SWMTotErrAdj)
cor.test(merge6$SingleParentScore, merge6$SWMTotErrAdj)
cor.test(merge6$Wave1SingleParent, merge6$SWMTotErrAdj)

pheno6e <- as.data.frame(cbind(pheno6$LUCIANO2_FID, pheno6$PNUM, pheno6$EACAQ00, pheno6$FDNVQ00, pheno6$FDACAQ00))
colnames(pheno6e) <- c("LUCIANO2_FID", "PNUM", "Wave5EA", "Wave6EA", "AllWavesEA")

pheno6e$Wave5EA[pheno6e$Wave5EA == 96] <- -1
pheno6e$Wave5EA[pheno6e$Wave5EA == 95] <- -2
pheno6e$Wave5EA[is.na(pheno6e$Wave5EA)] <- -2
pheno6e$Wave6EA[pheno6e$Wave6EA == 96] <- -1
pheno6e$Wave6EA[pheno6e$Wave6EA == 95] <- -2
pheno6e$Wave6EA[is.na(pheno6e$Wave6EA)] <- -2

pheno6e$AllWavesEA2 <- 0
i <- 1
for (i in  1:19452) {
  if (pheno6e$Wave5EA[i] > pheno6e$Wave6EA[i]) {pheno6e$AllWavesEA2[i] <- pheno6e$Wave5EA[i]}
  else {pheno6e$AllWavesEA2[i] <- pheno6e$Wave6EA[i]}
}

pheno6e$AllWavesEA2[pheno6e$AllWavesEA2 == -2] <- NA
pheno6e$AllWavesEA2[pheno6e$AllWavesEA2 == -1] <- 1

pheno6e <- merge(pheno1sub,pheno6e)
pheno6e <- pheno6e[pheno6e$MFC == "M",]
pheno6e <- pheno6e[!is.na(pheno6e$AllWavesEA2),]
merge7 <- merge(pheno6e, pheno8e, by="LUCIANO2_FID", all=TRUE)
merge7 <- merge7[,-3]
names(merge7)[names(merge7) == 'LUCIANO2_SID.y'] <- 'LUCIANO2_SID'
merge7$LUCIANO2_SID[merge7$LUCIANO2_SID == ""] <- NA
merge7 <- merge7[!is.na(merge7$LUCIANO2_SID),]
merge7$AllWavesEA2 <- as.numeric(merge7$AllWavesEA2)

cor.test(merge7$AllWavesEA2, merge7$SWMTotErrAdj)

library(ggpubr)
ggdensity(pheno8e$SWMTotErrAdj, 
          main = "SWM Task",
          xlab = "Total errors")

ggdensity(CogModelW5$BETA, 
          main = "CG Task",
          xlab = "Beta")

ggdensity(CogModelW6$BETA, 
          main = "CG Task",
          xlab = "Beta")

cor.test(pheno8e$CGTDELAY, pheno8e$SWMTotErrAdj)

merge8 <- merge(pheno7e, CogModelW5)
merge9 <- merge(pheno9e, CogModelW5)
merge10 <- merge(pheno10e, CogModelW5, by="LUCIANO2_FID")
merge11 <- merge(pheno11e, CogModelW5, by="LUCIANO2_FID")

merge8$LUCIANO2_SID[merge8$LUCIANO2_SID == ""] <- NA
merge8 <- merge8[!is.na(merge8$LUCIANO2_SID),]
merge8 <- merge8[!is.na(merge8$IMD_Decile),]
merge9$LUCIANO2_SID[merge9$LUCIANO2_SID == ""] <- NA
merge9 <- merge9[!is.na(merge9$LUCIANO2_SID),]
merge10 <- merge10[,-2]
names(merge10)[names(merge10) == 'LUCIANO2_SID.y'] <- 'LUCIANO2_SID'
merge10$LUCIANO2_SID[merge10$LUCIANO2_SID == ""] <- NA
merge10 <- merge10[!is.na(merge10$LUCIANO2_SID),]
merge10$AnySmoking <- as.numeric(merge10$AnySmoking)
merge10$BETAdjst <- as.numeric(merge10$BETAdjst)
merge11$LUCIANO2_SID[merge11$LUCIANO2_SID == ""] <- NA
merge11 <- merge11[!is.na(merge11$LUCIANO2_SID),]
merge11$SingleParentScore <- as.numeric(merge11$SingleParentScore)
merge11$Wave1SingleParent <- as.numeric(merge11$Wave1SingleParent)
merge12 <- merge(pheno6e, CogModelW5, by="LUCIANO2_FID")
merge12 <- merge12[,-3]
names(merge12)[names(merge12) == 'LUCIANO2_SID.y'] <- 'LUCIANO2_SID'
merge12$LUCIANO2_SID[merge12$LUCIANO2_SID == ""] <- NA
merge12 <- merge12[!is.na(merge12$LUCIANO2_SID),]
merge12$AllWavesEA2 <- as.numeric(merge12$AllWavesEA2)

cor.test(merge8$IMD_Decile, merge8$BETAdjst)
cor.test(merge9$RC1, merge9$BETAdjst)
cor.test(merge9$RC2, merge9$BETAdjst)
cor.test(merge9$RC3, merge9$BETAdjst)
cor.test(merge10$AnySmoking, merge10$BETAdjst)
cor.test(merge10$SMOKINGSEVERITY, merge10$BETAdjst)
cor.test(merge11$SingleParentScore, merge11$BETAdjst)
cor.test(merge11$Wave1SingleParent, merge11$BETAdjst)
cor.test(merge12$AllWavesEA2, merge12$BETAdjst)

Q <- quantile(CogModelW6$BETA, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(CogModelW6$BETA)
CogModelW6<- subset(CogModelW6, CogModelW6$BETA > (Q[1] - 1.5*iqr) & CogModelW6$BETA < (Q[2]+1.5*iqr))

CANTAB2sub <- CANTAB2[,-5:-11]
CANTAB2sub$FCCATECH[is.na(CANTAB2sub$FCCATECH)] <- 0
CANTAB2sub$FCTIRC00[is.na(CANTAB2sub$FCTIRC00)] <- 3

CogModelW6 <- merge(CogModelW6, CANTAB2sub)
CogModelW6 <- CogModelW6[!is.na(CogModelW6$FCMCS6AG),]

COG3.lm <- lm(BETA ~ FCMCS6AG + FCCSEX00 + Distraction + Interruption + FCCATECH + FCTIRC00, data=CogModelW6)
anova(COG3.lm)
summary(COG3.lm)
CogModelW6$BETAdjst <- resid(COG3.lm)

CANTAB2 <- CANTAB2[!is.na(CANTAB2$CGTDELAY),]
Q <- quantile(CANTAB2$CGTDELAY, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(CANTAB2$CGTDELAY)
CANTAB2<- subset(CANTAB2, CANTAB2$CGTDELAY > (Q[1] - 1.5*iqr) & CANTAB2$CGTDELAY < (Q[2]+1.5*iqr))

CANTAB2sub <- CANTAB2[,-5:-11]
CANTAB2sub$FCCATECH[is.na(CANTAB2sub$FCCATECH)] <- 0
CANTAB2sub$FCTIRC00[is.na(CANTAB2sub$FCTIRC00)] <- 3

CogModelW6 <- merge(CogModelW6, CANTAB2sub)
CogModelW6 <- CogModelW6[!is.na(CogModelW6$FCMCS6AG),]

COG3.lm <- lm(BETA ~ FCMCS6AG + FCCSEX00 + Distraction + Interruption + FCCATECH + FCTIRC00, data=CogModelW6)
anova(COG3.lm)
summary(COG3.lm)
CogModelW6$BETAdjst <- resid(COG3.lm)

merge13 <- merge(pheno7e, CogModelW6)
merge14 <- merge(pheno9e, CogModelW6)
merge15 <- merge(pheno10e, CogModelW6, by="LUCIANO2_FID")
merge16 <- merge(pheno11e, CogModelW6, by="LUCIANO2_FID")

merge13$LUCIANO2_SID[merge13$LUCIANO2_SID == ""] <- NA
merge13 <- merge13[!is.na(merge13$LUCIANO2_SID),]
merge13 <- merge13[!is.na(merge13$IMD_Decile),]
merge14$LUCIANO2_SID[merge14$LUCIANO2_SID == ""] <- NA
merge14 <- merge14[!is.na(merge14$LUCIANO2_SID),]
merge15 <- merge15[,-2]
names(merge15)[names(merge15) == 'LUCIANO2_SID.y'] <- 'LUCIANO2_SID'
merge15$LUCIANO2_SID[merge15$LUCIANO2_SID == ""] <- NA
merge15 <- merge15[!is.na(merge15$LUCIANO2_SID),]
merge15$AnySmoking <- as.numeric(merge15$AnySmoking)
merge15$BETAdjst <- as.numeric(merge15$BETAdjst)
merge16$LUCIANO2_SID[merge16$LUCIANO2_SID == ""] <- NA
merge16 <- merge16[!is.na(merge16$LUCIANO2_SID),]
merge16$SingleParentScore <- as.numeric(merge16$SingleParentScore)
merge16$Wave1SingleParent <- as.numeric(merge16$Wave1SingleParent)
merge17 <- merge(pheno6e, CogModelW6, by="LUCIANO2_FID")
merge17 <- merge17[,-3]
names(merge17)[names(merge17) == 'LUCIANO2_SID.y'] <- 'LUCIANO2_SID'
merge17$LUCIANO2_SID[merge17$LUCIANO2_SID == ""] <- NA
merge17 <- merge17[!is.na(merge17$LUCIANO2_SID),]
merge17$AllWavesEA2 <- as.numeric(merge17$AllWavesEA2)

cor.test(merge13$IMD_Decile, merge13$BETAdjst)
cor.test(merge14$RC1, merge14$BETAdjst)
cor.test(merge14$RC2, merge14$BETAdjst)
cor.test(merge14$RC3, merge14$BETAdjst)
cor.test(merge15$AnySmoking, merge15$BETAdjst)
cor.test(merge15$SMOKINGSEVERITY, merge15$BETAdjst)
cor.test(merge16$SingleParentScore, merge16$BETAdjst)
cor.test(merge16$Wave1SingleParent, merge16$BETAdjst)
cor.test(merge17$AllWavesEA2, merge17$BETAdjst)

genolabels <- fread('C:/Users/lucas/Documents/documents_20230428/GDAC_2022_18_LUCIANO_mcs_genoID_projectID_link.txt')
colnames(genolabels) <- c("geno_fid","geno_sid","LUCIANO2_FID","LUCIANO2_SID")

merge1 <- merge(merge1,genolabels)
merge2 <- merge(merge2,genolabels)
merge3 <- merge(merge3,genolabels)
merge6 <- merge(merge6,genolabels)
merge7 <- merge(merge7,genolabels)
merge8 <- merge(merge8,genolabels)
merge9 <- merge(merge9,genolabels)
merge10 <- merge(merge10,genolabels)
merge11 <- merge(merge11,genolabels)
merge12 <- merge(merge12,genolabels)
merge13 <- merge(merge13,genolabels)
merge14 <- merge(merge14,genolabels)
merge15 <- merge(merge15,genolabels)
merge16 <- merge(merge16,genolabels)
merge17 <- merge(merge17,genolabels)
covar <- merge(CANTABsub,genolabels)

merge1 <- merge1[, c(24, 25, 4, 3)]
merge2 <- merge2[, c(26, 27, 6, 3, 4, 5)]
merge3 <- merge3[, c(26, 27, 6, 4, 5)]
merge6 <- merge6[, c(25, 26, 5, 3, 4)]
merge7 <- merge7[, c(29, 30, 9, 8)]

merge8 <- merge8[, c(28, 29, 14, 3)]
merge9 <- merge9[, c(30, 31, 29, 3, 4, 5)]
merge10 <- merge10[, c(29, 30, 28, 3, 4)]
merge11 <- merge11[, c(29, 30, 28, 3, 4)]
merge12 <- merge12[, c(32, 33, 18, 7)]
merge13 <- merge13[, c(16, 17, 15, 3)]
merge14 <- merge14[, c(18, 19, 17, 3, 4, 5)]
merge15 <- merge15[, c(17, 18, 16, 3, 4)]
merge16 <- merge16[, c(17, 18, 16, 3, 4)]
merge17 <- merge17[, c(21, 22, 14, 7)]
covar <- covar[,c(21,22,19,20,14,15,5)]
covar$EMCS5AGE[is.na(covar$EMCS5AGE)] <- mean(covar$EMCS5AGE, na.rm = TRUE)


PCA <- fread('C:/Users/lucas/Documents/GRM/Wave5PCA.txt')
colnames(PCA) <- c("geno_fid","geno_sid","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
covar <- merge(covar,PCA)

write.table(merge1, file = 'SWMIMD.txt', row.names=FALSE, col.names=FALSE, quote = FALSE, sep = " ")
write.table(merge2, file = 'SWMHOME.txt', row.names=FALSE, col.names=FALSE, quote = FALSE, sep = " ")
write.table(merge3, file = 'SWMSMOKE.txt', row.names=FALSE, col.names=FALSE, quote = FALSE, sep = " ")
write.table(merge6, file = 'SWMSINGLE.txt', row.names=FALSE, col.names=FALSE, quote = FALSE, sep = " ")
write.table(merge7, file = 'SWMEA.txt', row.names=FALSE, col.names=FALSE, quote = FALSE, sep = " ")

write.table(merge8, file = 'CGT5IMD.txt', row.names=FALSE, col.names=FALSE, quote = FALSE, sep = " ")
write.table(merge9, file = 'CGT5HOME.txt', row.names=FALSE, col.names=FALSE, quote = FALSE, sep = " ")
write.table(merge10, file = 'CGT5SMOKE.txt', row.names=FALSE, col.names=FALSE, quote = FALSE, sep = " ")
write.table(merge11, file = 'CGT5SINGLE.txt', row.names=FALSE, col.names=FALSE, quote = FALSE, sep = " ")
write.table(merge12, file = 'CGT5EA.txt', row.names=FALSE, col.names=FALSE, quote = FALSE, sep = " ")

write.table(merge13, file = 'CGT6IMD.txt', row.names=FALSE, col.names=FALSE, quote = FALSE, sep = " ")
write.table(merge14, file = 'CGT6HOME.txt', row.names=FALSE, col.names=FALSE, quote = FALSE, sep = " ")
write.table(merge15, file = 'CGT6SMOKE.txt', row.names=FALSE, col.names=FALSE, quote = FALSE, sep = " ")
write.table(merge16, file = 'CGT6SINGLE.txt', row.names=FALSE, col.names=FALSE, quote = FALSE, sep = " ")
write.table(merge17, file = 'CGT6EA.txt', row.names=FALSE, col.names=FALSE, quote = FALSE, sep = " ")
write.table(covar, file = 'COVAR.txt', row.names=FALSE, col.names=FALSE, quote = FALSE, sep = " ")

merge0 <- merge1[,-4]
write.table(merge0, file = 'SWM.txt', row.names=FALSE, col.names=FALSE, quote = FALSE)


Wave5IDsSWM <- pheno8e[,1:2]
Wave5IDsCGT <- CogModelW5[,1:2]
Wave5IDs <- rbind(Wave5IDsSWM, Wave5IDsCGT )
Wave5IDs <- unique(Wave5IDs)
Wave5IDs <- merge(Wave5IDs,genolabels)
Wave5IDs <- Wave5IDs[,3:4]

Wave6IDsCGT <- CogModelW6[,1:2]
Wave6IDs <- merge(Wave6IDsCGT,genolabels)
Wave6IDs <- Wave6IDs[,3:4]

write.table(Wave5IDs, file = 'Wave5IDs.txt', row.names=FALSE, col.names=FALSE, quote = FALSE, sep = " ")
write.table(Wave6IDs, file = 'Wave6IDs.txt', row.names=FALSE, col.names=FALSE, quote = FALSE, sep = " ")

pheno1$samples_2_use <- as.character(pheno1$samples_2_use)

motha <- pheno1[pheno1$mfc == "M",]
motha <- motha[motha$samples_2_use == "1",]
motha <- motha[,1:2]
colnames(motha) <- toupper(colnames(motha))
motha <- merge(motha, genolabels)
motha <- motha[,3:4]

fatha <- pheno1[pheno1$mfc == "F",]
fatha <- fatha[fatha$samples_2_use == "1",]
fatha <- fatha[,1:2]
colnames(fatha) <- toupper(colnames(fatha))
fatha <- merge(fatha, genolabels)
fatha <- fatha[,3:4]

baba <- pheno1[pheno1$mfc == "C",]
baba <- baba[baba$samples_2_use == "1",]
baba <- baba[,1:2]
colnames(baba) <- toupper(colnames(baba))
baba <- merge(baba, genolabels)
baba <- baba[,3:4]

hereIamat <- merge(baba,fatha,by="geno_fid",all=TRUE)
colnames(hereIamat) <- c("geno_fid","individualID","fatherID")
campgrenada <- merge(hereIamat,motha,by="geno_fid",all=TRUE)
colnames(campgrenada) <- c("geno_fid","geno_sid","fatherID","motherID")
campgrenada <- merge(campgrenada,genolabels,all=TRUE)
campisvery <- pheno1[,1:5]
campisvery <- campisvery[,-3:-4]
colnames(campisvery) <- toupper(colnames(campisvery))
entertaining <- merge(campgrenada,campisvery)
entertaining <- entertaining[,-1:-2]
entertaining$andtheysaywellhavesomefunifitstopsraining <- -9
entertaining <- entertaining[!is.na(entertaining$geno_sid),]
entertaining$fatherID[is.na(entertaining$fatherID)] <- 0
entertaining$motherID[is.na(entertaining$motherID)] <- 0

fam <- fread('C:/Users/lucas/Documents/documents_20230428/MCS_topmed.fam')
fam <- fam[,1:2]
colnames(fam) <- c("geno_fid","geno_sid")
famfinal <- merge(fam,entertaining)

write.table(famfinal, file = 'MCS_topmed.fam', row.names=FALSE, col.names=FALSE, quote = FALSE, sep = " ")

rstudioapi::writeRStudioPreference("data_viewer_max_columns", 1000L)
