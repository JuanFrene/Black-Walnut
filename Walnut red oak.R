library(vegan)
library(agricolae)
library(RColorBrewer)
library(gplots)
library(reshape)
library(stringr)
library(ggplot2)
library(reshape2)
library(Hmisc)
library(stats)
library(ggfortify)
library(plyr)
library(cluster)
library(ggthemes)
library(ggcorrplot)
library(lavaan)
library(semPlot)
library(lavaanPlot)
library(factoextra)
library(ggpubr)

Fame <- read.table("clipboard", header = TRUE)
TabCom <- read.table("clipboard", header = TRUE)
TabPro <- read.table("clipboard", header = TRUE)
TabNor <- read.table("clipboard", header = TRUE)

str(TabCom)
head(TabCom)
summary(TabCom)
ncol(Fame)


TabCom$Plant = factor(TabCom$Plant, levels=c("BW", "RO","BS"))
TabCom0<- TabCom[TabCom$Depth=="0-10",]
TabCom10<- TabCom[TabCom$Depth=="10-20",]
TabPro$Plant = factor(TabPro$Plant, levels=c("BW","BS","RO"))
TabPro0<- TabPro[TabPro$Depth=="0-10",]
TabPro10<- TabPro[TabPro$Depth=="10-20",]
TabPro0$Sampling = factor(TabPro0$Sampling, levels=c("Autumn","Winter","Spring", "Summer"))
TabPro10$Sampling = factor(TabPro10$Sampling, levels=c("Autumn","Winter","Spring", "Summer"))
TabNorBW<- TabNor[TabNor$Plant=="BW",]
TabNorRO<- TabNor[TabNor$Plant=="RO",]



TabPro$Treatment = factor(TabPro$Treatment, levels=c("BW", "RO","DZ_BW", "DZ_RO"))
FAME <- as.matrix(Fame[,5:75])
mdf2$Treatment=as.factor(mdf2$Treatment)
df2$Treatment = factor(df2$Treatment, levels=c("BW", "RO","DZ_BW", "DZ_RO"))

#TOC
ggplot(data=TabPro0, aes(x=Sampling, y=TOC, fill=Plant)) +
  geom_bar(stat="identity", position='dodge')+
  theme_few()+ scale_fill_manual(values=c('Red','blue','green'))+ ggtitle("Total organic carbon")+ theme(
    plot.title = element_text(color="Black", size=14, face="bold.italic"))+ 
  xlab("Season") + ylab("%") 
#+ geom_errorbar(aes(ymin=TOC-DS_TOC, ymax=TOC+DS_TOC)) 
  

ggplot(data=TabPro10, aes(x=Sampling, y=X., fill=Plant)) +
  geom_bar(stat="identity", position='dodge')+
  theme_minimal()+ scale_fill_manual(values=c('Red','blue','green'))+ ggtitle("Soil moisture")+ theme(
    plot.title = element_text(color="Black", size=14, face="bold.italic"))+ 
  xlab("Season") + ylab("%") + ylim(0,25)

aovTOC <- aov(TOC~Plant*Season, data=TabCom10)
summary(aovTOC)
lsdaovTOC <- LSD.test(aovTOC, c("Plant", "Season"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaovTOC
aovTOC10 <- aov(TOC~Plant, data=TabCom10[37:48,])
summary(aovTOC10)
lsdaaovpH <- LSD.test(aovTOC10, c("Plant"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovpH
aovpH <- aov(pH~Plant, data=TabCom0[36:48,])
summary(aovpH)
lsdaaovpH <- LSD.test(aovpH, c("Plant"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovpH

#TC
ggplot(data=TabPro0, aes(x=Pant, y=TOC, fill=Depth)) +
  geom_bar(stat="identity", position='stack')+
  theme_minimal()+ scale_fill_manual(values=c('Red','blue','green'))+ ggtitle("Total organic carbon")+ theme(
    plot.title = element_text(color="Black", size=14, face="bold.italic"))+ 
  xlab("Season") + ylab("%")

ggplot(data=TabPro10, aes(x=Sampling, y=X., fill=Plant)) +
  geom_bar(stat="identity", position='dodge')+
  theme_minimal()+ scale_fill_manual(values=c('Red','blue','green'))+ ggtitle("Soil moisture")+ theme(
    plot.title = element_text(color="Black", size=14, face="bold.italic"))+ 
  xlab("Season") + ylab("%") + ylim(0,25)

aovTC <- aov(TC~Plant*Season, data=TabCom10)
summary(aovTC)
lsdaovTC <- LSD.test(aovTOC, c("Plant", "Season"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaovTC
aovTC <- aov(TC~Plant, data=TabCom10[1:12,])
summary(aovTC)
lsdovTC <- LSD.test(aovTC, c("Plant"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdovTC

#TN
ggplot(data=TabPro0, aes(x=Pant, y=TN, fill=Depth)) +
  geom_bar(stat="identity", position='stack')+
  theme_minimal()+ scale_fill_manual(values=c('Red','blue','green'))+ ggtitle("Total Nitrogen")+ theme(
    plot.title = element_text(color="Black", size=14, face="bold.italic"))+ 
  xlab("Season") + ylab("%")

ggplot(data=TabPro10, aes(x=Sampling, y=TN, fill=Plant)) +
  geom_bar(stat="identity", position='dodge')+
  theme_minimal()+ scale_fill_manual(values=c('Red','blue','green'))+ ggtitle("Total Nitrogen")+ theme(
    plot.title = element_text(color="Black", size=14, face="bold.italic"))+ 
  xlab("Season") + ylab("%") + ylim(0,25)

aovTN <- aov(TN~Plant*Season, data=TabCom10)
summary(aovTN)
lsdaovTN <- LSD.test(aovTN, c("Plant", "Season"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaovTN
aovTN <- aov(TN~Plant, data=TabCom10[37:48,])
summary(aovTN)
lsdovTN <- LSD.test(aovTN, c("Plant"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdovTN

#C.N
ggplot(data=TabPro0, aes(x=Pant, y=TN, fill=Depth)) +
  geom_bar(stat="identity", position='stack')+
  theme_minimal()+ scale_fill_manual(values=c('Red','blue','green'))+ ggtitle("Total Nitrogen")+ theme(
    plot.title = element_text(color="Black", size=14, face="bold.italic"))+ 
  xlab("Season") + ylab("%")

ggplot(data=TabPro10, aes(x=Sampling, y=TN, fill=Plant)) +
  geom_bar(stat="identity", position='dodge')+
  theme_minimal()+ scale_fill_manual(values=c('Red','blue','green'))+ ggtitle("Total Nitrogen")+ theme(
    plot.title = element_text(color="Black", size=14, face="bold.italic"))+ 
  xlab("Season") + ylab("%") + ylim(0,25)

aovC.N <- aov(C.N~Plant*Season, data=TabCom10)
summary(aovC.N)
lsdaovC.N <- LSD.test(aovTN, c("Plant", "Season"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaovC.N
aovC.N <- aov(C.N~Plant, data=TabCom10[37:48,])
summary(aovC.N)
lsdaovC.N <- LSD.test(aovC.N, c("Plant"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaovC.N


#Soil moisture
ggplot(data=TabPro0, aes(x=Pant, y=X., fill=Depth)) +
  geom_bar(stat="identity", position='stack')+
  theme_minimal()+ scale_fill_manual(values=c('Red','blue','green'))+ ggtitle("Soil moisture")+ theme(
    plot.title = element_text(color="Black", size=14, face="bold.italic"))+ 
  xlab("Season") + ylab("%")

ggplot(data=TabPro10, aes(x=Sampling, y=X., fill=Plant)) +
  geom_bar(stat="identity", position='dodge')+
  theme_minimal()+ scale_fill_manual(values=c('Red','blue','green'))+ ggtitle("Soil moisture")+ theme(
    plot.title = element_text(color="Black", size=14, face="bold.italic"))+ 
  xlab("Season") + ylab("%") + ylim(0,25)

aovpH <- aov(Hum~Treatment*Season, data=TabCom)
summary(aovpH)
lsdaaovpH <- LSD.test(aovpH, c("Treatment", "Season"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovpH
aovpH <- aov(Hum~Plant, data=TabCom0[36:48,])
summary(aovpH)
lsdaaovpH <- LSD.test(aovpH, c("Plant"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovpH
aovpH <- aov(pH~Plant, data=TabCom0[36:48,])
summary(aovpH)
lsdaaovpH <- LSD.test(aovpH, c("Plant"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovpH

#pH
ggplot(data=TabPro0, aes(x=Sampling, y=pH, fill=Plant)) +
  geom_bar(stat="identity", position='dodge')+
  theme_minimal()+ scale_fill_manual(values=c('Red','blue','green'))+ ggtitle("pH")+ theme(
  plot.title = element_text(color="Black", size=14, face="bold.italic"))+ 
  xlab("Season")

ggplot(data=TabPro10, aes(x=Sampling, y=pH, fill=Plant)) +
  geom_bar(stat="identity", position='dodge')+
  theme_minimal()+ scale_fill_manual(values=c('Red','blue','green'))+ ggtitle("pH")+ theme(
    plot.title = element_text(color="Black", size=14, face="bold.italic"))+ 
  xlab("Season")


aovpH <- aov(pH~Treatment*Season, data=TabCom0)
summary(aovpH)
lsdaaovpH <- LSD.test(aovpH, c("Treatment", "Season"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovpH
aovpH <- aov(pH~Treatment, data=TabCom0[36:48,])
summary(aovpH)
lsdaaovpH <- LSD.test(aovpH, c("Treatment"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovpH
aovpH <- aov(pH~Plant, data=TabCom0[36:48,])
summary(aovpH)
lsdaaovpH <- LSD.test(aovpH, c("Plant"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovpH

aovpH <- aov(pH~Treatment*Season, data=TabCom10)
summary(aovpH)
lsdaaovpH <- LSD.test(aovpH, c("Treatment", "Season"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovpH

aovpH <- aov(pH~Treatment, data=TabCom10[36:48,])
summary(aovpH)
lsdaaovpH <- LSD.test(aovpH, c("Treatment"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovpH
aovpH <- aov(pH~Plant, data=TabCom10[36:48,])
summary(aovpH)
lsdaaovpH <- LSD.test(aovpH, c("Plant"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovpH


#BG
ggplot(data=TabPro0, aes(x=Sampling, y=BG, fill=Plant)) +
  geom_bar(stat="identity", position='dodge', color='black')+
  theme_few()+ scale_fill_manual(values=c('black','gray90','gray'))+ ggtitle("Glucosidase")+ theme(
    plot.title = element_text(color="Black", size=14, face="bold.italic"))+ 
  xlab("Season")+ ylab("mg p-nitrofenol / kg soil h")+
  geom_errorbar(position = "dodge", aes(ymin=BG-0, ymax=BG+DS_BG), width = 0.9, size=0.1)

ggplot(data=TabPro10, aes(x=Sampling, y=BG, fill=Plant)) +
  geom_bar(stat="identity", position='dodge', color='black')+
  theme_few()+ scale_fill_manual(values=c('black','gray90','gray'))+ ggtitle("Glucosidase")+ theme(
    plot.title = element_text(color="Black", size=14, face="bold.italic"))+ 
  xlab("Season")+ ylab("mg p-nitrofenol / kg soil h")+
  scale_linetype_manual(values = 'black')+
  geom_errorbar(position = "dodge", aes(ymin=BG-0, ymax=BG+DS_BG), width = 0.9, size=0.1)


aovBG0 <- aov(BG~Treatment*Season, data=TabCom0)
summary(aovBG0)
lsdaaovBG0 <- LSD.test(aovBG0, c("Treatment", "Season"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovBG0
aovBG0 <- aov(BG~Treatment, data=TabCom0[36:48,])
summary(aovBG0)
lsdaaovBG0 <- LSD.test(aovBG0, c("Treatment"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovBG0
aovBG0 <- aov(BG~Plant, data=TabCom0[36:48,])
summary(aovBG0)
lsdaaovBG0 <- LSD.test(aovBG0, c("Plant"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovBG0

aovBG10 <- aov(BG~Treatment*Season, data=TabCom10)
summary(aovBG10)
aovBG10 <- aov(BG~Plant*Season, data=TabCom10)
summary(aovBG10)
aovBG10 <- aov(BG~Treatment, data=TabCom10[36:48,])
summary(aovBG10)
lsdaaovBG10 <- LSD.test(aovBG10, c("Treatment"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovBG10
aovBG10 <- aov(BG~Plant, data=TabCom10[36:48,])
summary(aovBG10)
lsdaaovBG10 <- LSD.test(aovBG10, c("Plant"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovBG10

#NAG
 ggplot(data=TabPro0, aes(x=Sampling, y=NAG, fill=Plant)) +
  geom_bar(stat="identity", position='dodge', color='black')+
  theme_few()+ scale_fill_manual(values=c('black','gray90','gray'))+ ggtitle("Glucosaminidase")+ theme(
    plot.title = element_text(color="Black", size=14, face="bold.italic"))+ 
  xlab("Season")+ ylab("mg p-nitrofenol / kg soil h")+
   geom_errorbar(position = "dodge", aes(ymin=NAG-0, ymax=NAG+DS_NAG), width = 0.9, size=0.1)
 

ggplot(data=TabPro10, aes(x=Sampling, y=NAG, fill=Plant)) +
  geom_bar(stat="identity", position='dodge', color='black')+
  theme_few()+ scale_fill_manual(values=c('black','gray90','gray'))+ ggtitle("Glucosamidase")+ theme(
    plot.title = element_text(color="Black", size=14, face="bold.italic"))+ 
  xlab("Season")+ ylab("mg p-nitrofenol / kg soil h")+
  geom_errorbar(position = "dodge", aes(ymin=NAG-0, ymax=NAG+10), width = 0.9, size=0.1)


aovNAG0 <- aov(NAG~Treatment*Season, data=TabCom0)
summary(aovNAG0)
lsdaovNAG0 <- LSD.test(aovNAG0, c("Treatment", "Season"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaovNAG0
aovNAG0 <- aov(NAG~Treatment, data=TabCom0[36:48,])
summary(aovNAG0)
lsdaovNAG0 <- LSD.test(aovNAG0, c("Treatment"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaovNAG0
aovNAG0 <- aov(NAG~Plant, data=TabCom0[36:48,])
summary(aovNAG0)
lsdaovNAG0 <- LSD.test(aovNAG0, c("Plant"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaovNAG0

aovNAG10 <- aov(NAG~Treatment*Season, data=TabCom10)
summary(aovNAG10)
aovNAG10 <- aov(NAG~Treatment, data=TabCom10[36:48,])
summary(aovNAG10)
lsdaovNAG10 <- LSD.test(aovNAG10, c("Treatment"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaovNAG10
aovNAG10 <- aov(NAG~Plant, data=TabCom10[36:48,])
summary(aovNAG10)
lsdaovNAG10 <- LSD.test(aovNAG10, c("Plant"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaovNAG10

#PME
ggplot(data=TabPro0, aes(x=Sampling, y=PME, fill=Plant)) +
  geom_bar(stat="identity", position='dodge', color='black')+
  theme_few()+ scale_fill_manual(values=c('black','gray90','gray'))+ ggtitle("Acid Phosphatase")+ theme(
    plot.title = element_text(color="Black", size=14, face="bold.italic"))+ 
  xlab("Season")+ ylab("mg p-nitrofenol / kg soil h")+
  geom_errorbar(position = "dodge", aes(ymin=PME-0, ymax=PME+DS_PME), width = 0.9, size=0.1)


ggplot(data=TabPro10, aes(x=Sampling, y=PME, fill=Plant)) +
  geom_bar(stat="identity", position='dodge', color='black')+
  theme_few()+ scale_fill_manual(values=c('black','gray90','gray'))+ ggtitle("Acid Phosphatase")+ theme(
    plot.title = element_text(color="Black", size=14, face="bold.italic"))+ 
  xlab("Season")+ ylab("mg p-nitrofenol / kg soil h")+
  geom_errorbar(position = "dodge", aes(ymin=PME-0, ymax=PME+DS_PME), width = 0.9, size=0.1)


aovPME0 <- aov(PME~Treatment*Season, data=TabCom0)
summary(aovPME0)
aovPME0 <- aov(PME~Treatment, data=TabCom0[36:48,])
summary(aovPME0)
lsdaovPME0 <- LSD.test(aovPME0, c("Treatment"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaovPME0
aovPME0 <- aov(PME~Plant, data=TabCom0[36:48,])
summary(aovPME0)
lsdaovPME0 <- LSD.test(aovPME0, c("Plant"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaovPME0

aovPME10 <- aov(PME~Treatment*Season, data=TabCom10)
summary(aovPME10)
aovPME10 <- aov(PME~Treatment, data=TabCom10[36:48,])
summary(aovPME10)
lsdaovPME10 <- LSD.test(aovPME10, c("Treatment"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaovPME10
aovPME10 <- aov(PME~Plant, data=TabCom10[36:48,])
summary(aovPME10)
lsdaovPME10 <- LSD.test(aovPME10, c("Plant"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaovPME10
aovPME

#AS
ggplot(data=TabPro0, aes(x=Sampling, y=AS, fill=Plant)) +
geom_bar(stat="identity", position='dodge', color='black')+
  theme_few()+ scale_fill_manual(values=c('black','gray90','gray'))+ ggtitle("Arysulphatase")+ theme(
    plot.title = element_text(color="Black", size=14, face="bold.italic"))+ 
  xlab("Season")+ ylab("mg p-nitrofenol / kg soil h")+
  geom_errorbar(position = "dodge", aes(ymin=AS-0, ymax=AS+DS_AS), width = 0.9, size=0.1)


ggplot(data=TabPro10, aes(x=Sampling, y=AS, fill=Plant)) +
  geom_bar(stat="identity", position='dodge', color='black')+
  theme_few()+ scale_fill_manual(values=c('black','gray90','gray'))+ ggtitle("Arysulphatase")+ theme(
    plot.title = element_text(color="Black", size=14, face="bold.italic"))+ 
  xlab("Season")+ ylab("mg p-nitrofenol / kg soil h")+
  geom_errorbar(position = "dodge", aes(ymin=AS-0, ymax=AS+DS_AS), width = 0.9, size=0.1)

aovAS0 <- aov(AS~Treatment*Season, data=TabCom0)
summary(aovAS0)
aovAS0 <- aov(AS~Treatment, data=TabCom0[36:48,])
summary(aovAS0)
lsdaovAS0 <- LSD.test(aovAS0, c("Treatment"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaovAS0
aovAS0 <- aov(AS~Plant, data=TabCom0[36:48,])
summary(aovAS0)
lsdaovAS0 <- LSD.test(aovAS0, c("Plant"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaovAS0

aovAS10 <- aov(AS~Treatment*Season, data=TabCom10)
summary(aovAS10)
aovAS10 <- aov(AS~Treatment, data=TabCom10[36:48,])
summary(aovAS10)
lsdaovAS10 <- LSD.test(aovAS10, c("Treatment"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaovAS10
aovAS10 <- aov(AS~Plant, data=TabCom10[36:48,])
summary(aovAS10)
lsdaovAS10 <- LSD.test(aovAS10, c("Plant"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaovAS10

#Combinar graficos
ggarrange(e1, e2, e3, e4, e5, e6, e7,e8, ncol=4, nrow=2, 
          labels=c('a','b','c','d','e','f','g','g'))


#Total FAME 
ggplot(data=TabPro0, aes(x=Sampling, y=Total_PLFA, fill=Plant)) +
  geom_bar(stat="identity", position='dodge')+
  theme_minimal()+ scale_fill_manual(values=c('Red','blue','green'))+ ggtitle("Total FAME")+ theme(
    plot.title = element_text(color="Black", size=14, face="bold.italic"))+ 
  xlab("Season")+ ylab("nmol / g soil")

ggplot(mapping = aes(x=Plant, y=Total_FAME)) +
  geom_bar(data = TabPro0[1:12,], width = 0.8, stat = "identity"+ scale_fill_manual(values=c('Red','blue','green')))+ 
  geom_bar(data = TabPro10[1:12,], width = 0.78, stat = "identity", fill = c('Red','blue','green')) +
  theme_few() 



ggplot(data=TabPro10, aes(x=Sampling, y=Total_FAME, fill=Plant)) +
  geom_bar(stat="identity", position='dodge', color='black')+
  theme_minimal()+ scale_fill_manual(values=c('Red','blue','green'))+ ggtitle("Total FAME")+ theme(
    plot.title = element_text(color="Black", size=14, face="bold.italic"))+ 
  xlab("Season")+ ylab("nmol / g soil")

aovTotal_PLFA0 <- aov(Total_PLFA~Treatment*Season, data=TabCom0)
summary(aovTotal_PLFA0)
lsdaaovTotal_PLFA0 <- LSD.test(aovTotal_PLFA0, c("Treatment", "Season"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovTotal_PLFA0
aovTotal_PLFA0 <- aov(Total_PLFA~Treatment, data=TabCom0[36:48,])
summary(aovTotal_PLFA0)
lsdaaovTotal_PLFA0 <- LSD.test(aovTotal_PLFA0, c("Treatment"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovTotal_PLFA0
aovTotal_PLFA0 <- aov(Total_PLFA~Plant, data=TabCom0[36:48,])
summary(aovTotal_PLFA0)
lsdaaovTotal_PLFA0 <- LSD.test(aovTotal_PLFA0, c("Plant"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovTotal_PLFA0

aovTotal_PLFA10 <- aov(Total_PLFA~Treatment*Season, data=TabCom10)
summary(aovTotal_PLFA10)
aovTotal_PLFA10 <- aov(Total_PLFA~Treatment, data=TabCom10[36:48,])
summary(aovTotal_PLFA10)
lsdaaovTotal_PLFA10 <- LSD.test(aovTotal_PLFA10, c("Treatment"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovTotal_PLFA10
aovTotal_PLFA10 <- aov(Total_PLFA~Plant, data=TabCom10[36:48,])
summary(aovTotal_PLFA10)
lsdaaovTotal_PLFA10 <- LSD.test(aovTotal_PLFA10, c("Plant"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovTotal_PLFA10

#Bacteria
ggplot(data=TabPro0, aes(x=Sampling, y=Bacteria, fill=Plant)) +
  geom_bar(stat="identity", position='dodge')+
  theme_minimal()+ scale_fill_manual(values=c('Red','blue','green'))+ ggtitle("Total Bacteria")+ theme(
    plot.title = element_text(color="Black", size=14, face="bold.italic"))+ 
  xlab("Season")+ ylab("nmol / g soil")

ggplot(data=TabPro10, aes(x=Sampling, y=Bacteria, fill=Plant)) +
  geom_bar(stat="identity", position='dodge')+
  theme_minimal()+ scale_fill_manual(values=c('Red','blue','green'))+ ggtitle("Total Bacteria")+ theme(
    plot.title = element_text(color="Black", size=14, face="bold.italic"))+ 
  xlab("Season")+ ylab("nmol / g soil")

aovBacteria0 <- aov(Bacteria~Treatment*Season, data=TabCom0)
summary(aovBacteria0)
lsdaaovBacteria0 <- LSD.test(aovBacteria0, c("Treatment", "Season"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovBacteria0
aovBacteria0 <- aov(Bacteria~Treatment, data=TabCom0[36:48,])
summary(aovBacteria0)
lsdaaovBacteria0 <- LSD.test(aovBacteria0, c("Treatment"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovBacteria0
aovBacteria0 <- aov(Bacteria~Plant, data=TabCom0[36:48,])
summary(aovBacteria0)
lsdaaovBacteria0 <- LSD.test(aovBacteria0, c("Plant"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovBacteria0

aovBacteria10 <- aov(Bacteria~Treatment*Season, data=TabCom10)
summary(aovBacteria10)
aovBacteria10 <- aov(Bacteria~Treatment, data=TabCom10[36:48,])
summary(aovBacteria10)
lsdaaovBacteria10 <- LSD.test(aovBacteria10, c("Treatment"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovBacteria10
aovBacteria10 <- aov(Bacteria~Plant, data=TabCom10[36:48,])
summary(aovBacteria10)
lsdaaovBacteria10 <- LSD.test(aovBacteria10, c("Plant"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovBacteria10

#Fungi
ggplot(data=TabPro0, aes(x=Sampling, y=Fungy, fill=Plant)) +
  geom_bar(stat="identity", position='dodge')+
  theme_minimal()+ scale_fill_manual(values=c('Red','blue','green'))+ ggtitle("Total Fungy")+ theme(
    plot.title = element_text(color="Black", size=14, face="bold.italic"))+ 
  xlab("Season")+ ylab("nmol / g soil")

ggplot(data=TabPro10, aes(x=Sampling, y=Fungy, fill=Plant)) +
  geom_bar(stat="identity", position='dodge')+
  theme_minimal()+ scale_fill_manual(values=c('Red','blue','green'))+ ggtitle("Total Fungy")+ theme(
    plot.title = element_text(color="Black", size=14, face="bold.italic"))+ 
  xlab("Season")+ ylab("nmol / g soil")

aovFungy0 <- aov(Fungy~Treatment*Season, data=TabCom0)
summary(aovFungy0)
lsdaaovFungy0 <- LSD.test(aovFungy0, c("Treatment", "Season"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovFungy0
aovFungy0 <- aov(Fungy~Treatment, data=TabCom0[36:48,])
summary(aovFungy0)
lsdaaovFungy0 <- LSD.test(aovFungy0, c("Treatment"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovFungy0
aovFungy0 <- aov(Fungy~Plant, data=TabCom0[36:48,])
summary(aovFungy0)
lsdaaovFungy0 <- LSD.test(aovFungy0, c("Plant"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovFungy0

aovFungy10 <- aov(Fungy~Treatment*Season, data=TabCom10)
summary(aovFungy10)
aovFungy10 <- aov(Fungy~Treatment, data=TabCom10[36:48,])
summary(aovFungy10)
lsdaaovFungy10 <- LSD.test(aovFungy10, c("Treatment"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovFungy10
aovFungy10 <- aov(Fungy~Plant, data=TabCom10[36:48,])
summary(aovFungy10)
lsdaaovFungy10 <- LSD.test(aovFungy10, c("Plant"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovFungy10

#Gram+
ggplot(data=TabPro0, aes(x=Sampling, y=GramPos, fill=Plant)) +
  geom_bar(stat="identity", position='dodge')+
  theme_minimal()+ scale_fill_manual(values=c('Red','blue','green'))+ ggtitle("Gram Positive")+ theme(
    plot.title = element_text(color="Black", size=14, face="bold.italic"))+ 
  xlab("Season")+ ylab("nmol / g soil")

ggplot(data=TabPro10, aes(x=Sampling, y=GramPos, fill=Plant)) +
  geom_bar(stat="identity", position='dodge')+
  theme_minimal()+ scale_fill_manual(values=c('Red','blue','green'))+ ggtitle("Gram Positive")+ theme(
    plot.title = element_text(color="Black", size=14, face="bold.italic"))+ 
  xlab("Season")+ ylab("nmol / g soil")

aovGramPos0 <- aov(GramPos~Treatment*Season, data=TabCom0)
summary(aovGramPos0)
lsdaaovGramPos0 <- LSD.test(aovGramPos0, c("Treatment", "Season"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovGramPos0
aovGramPos0 <- aov(GramPos~Treatment, data=TabCom0[36:48,])
summary(aovGramPos0)
lsdaaovGramPos0 <- LSD.test(aovGramPos0, c("Treatment"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovGramPos0
aovGramPos0 <- aov(GramPos~Plant, data=TabCom0[36:48,])
summary(aovGramPos0)
lsdaaovGramPos0 <- LSD.test(aovGramPos0, c("Plant"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovGramPos0

aovGramPos10 <- aov(GramPos~Treatment*Season, data=TabCom10)
summary(aovGramPos10)
aovGramPos10 <- aov(GramPos~Treatment, data=TabCom10[36:48,])
summary(aovGramPos10)
lsdaaovGramPos10 <- LSD.test(aovGramPos10, c("Treatment"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovGramPos10
aovGramPos10 <- aov(GramPos~Plant, data=TabCom10[36:48,])
summary(aovGramPos10)
lsdaaovGramPos10 <- LSD.test(aovGramPos10, c("Plant"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovGramPos10

#Gram-
ggplot(data=TabPro0, aes(x=Sampling, y=GramNeg, fill=Plant)) +
  geom_bar(stat="identity", position='dodge')+
  theme_minimal()+ scale_fill_manual(values=c('Red','blue','green'))+ ggtitle("Gram Negative")+ theme(
    plot.title = element_text(color="Black", size=14, face="bold.italic"))+ 
  xlab("Season")+ ylab("nmol / g soil")

ggplot(data=TabPro10, aes(x=Sampling, y=GramNeg, fill=Plant)) +
  geom_bar(stat="identity", position='dodge')+
  theme_minimal()+ scale_fill_manual(values=c('Red','blue','green'))+ ggtitle("Gram Negative")+ theme(
    plot.title = element_text(color="Black", size=14, face="bold.italic"))+ 
  xlab("Season")+ ylab("nmol / g soil")

aovGramNeg0 <- aov(GramNeg~Treatment*Season, data=TabCom0)
summary(aovGramNeg0)
lsdaaovGramNeg0 <- LSD.test(aovGramNeg0, c("Treatment", "Season"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovGramNeg0
aovGramNeg0 <- aov(GramNeg~Treatment, data=TabCom0[36:48,])
summary(aovGramNeg0)
lsdaaovGramNeg0 <- LSD.test(aovGramNeg0, c("Treatment"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovGramNeg0
aovGramNeg0 <- aov(GramNeg~Plant, data=TabCom0[36:48,])
summary(aovGramNeg0)
lsdaaovGramNeg0 <- LSD.test(aovGramNeg0, c("Plant"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovGramNeg0

aovGramNeg10 <- aov(GramNeg~Treatment*Season, data=TabCom10)
summary(aovGramNeg10)
aovGramNeg10 <- aov(GramNeg~Treatment, data=TabCom10[36:48,])
summary(aovGramNeg10)
lsdaaovGramNeg10 <- LSD.test(aovGramNeg10, c("Treatment"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovGramNeg10
aovGramNeg10 <- aov(GramNeg~Plant, data=TabCom10[36:48,])
summary(aovGramNeg10)
lsdaaovGramNeg10 <- LSD.test(aovGramNeg10, c("Plant"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovGramNeg10

#Actinomycetes
ggplot(data=TabPro0, aes(x=Sampling, y=Actinomycetes, fill=Plant)) +
  geom_bar(stat="identity", position='dodge')+
  theme_minimal()+ scale_fill_manual(values=c('Red','blue','green'))+ ggtitle("Actinomycetes")+ theme(
    plot.title = element_text(color="Black", size=14, face="bold.italic"))+ 
  xlab("Season")+ ylab("nmol / g soil")

ggplot(data=TabPro10, aes(x=Sampling, y=Actinomycetes, fill=Plant)) +
  geom_bar(stat="identity", position='dodge')+
  theme_minimal()+ scale_fill_manual(values=c('Red','blue','green'))+ ggtitle("Actinomycetes")+ theme(
    plot.title = element_text(color="Black", size=14, face="bold.italic"))+ 
  xlab("Season")+ ylab("nmol / g soil")

aovActinomycetes0 <- aov(Actinomycetes~Treatment*Season, data=TabCom0)
summary(aovActinomycetes0)
aovActinomycetes0 <- aov(Actinomycetes~Treatment, data=TabCom0[36:48,])
summary(aovActinomycetes0)
lsdaaovActinomycetes0 <- LSD.test(aovActinomycetes0, c("Treatment"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovActinomycetes0
aovActinomycetes0 <- aov(Actinomycetes~Plant, data=TabCom0[36:48,])
summary(aovActinomycetes0)
lsdaaovActinomycetes0 <- LSD.test(aovActinomycetes0, c("Plant"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovActinomycetes0

aovActinomycetes10 <- aov(Actinomycetes~Treatment*Season, data=TabCom10)
summary(aovActinomycetes10)
aovActinomycetes10 <- aov(Actinomycetes~Treatment, data=TabCom10[36:48,])
summary(aovActinomycetes10)
lsdaaovActinomycetes10 <- LSD.test(aovActinomycetes10, c("Treatment"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovActinomycetes10
aovActinomycetes10 <- aov(Actinomycetes~Plant, data=TabCom10[36:48,])
summary(aovActinomycetes10)
lsdaaovActinomycetes10 <- LSD.test(aovActinomycetes10, c("Plant"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovActinomycetes10

#AMF
ggplot(data=TabPro0, aes(x=Sampling, y=AMF, fill=Plant)) +
  geom_bar(stat="identity", position='dodge')+
  theme_minimal()+ scale_fill_manual(values=c('Red','blue','green'))+ ggtitle("AMF")+ theme(
    plot.title = element_text(color="Black", size=14, face="bold.italic"))+ 
  xlab("Season")+ ylab("nmol / g soil")

ggplot(data=TabPro10, aes(x=Sampling, y=AMF, fill=Plant)) +
  geom_bar(stat="identity", position='dodge')+
  theme_minimal()+ scale_fill_manual(values=c('Red','blue','green'))+ ggtitle("AMF")+ theme(
    plot.title = element_text(color="Black", size=14, face="bold.italic"))+ 
  xlab("Season")+ ylab("nmol / g soil")

aovAMF0 <- aov(AMF~Treatment*Season, data=TabCom0)
summary(aovAMF0)
aovAMF0 <- aov(AMF~Treatment, data=TabCom0[36:48,])
summary(aovAMF0)
lsdaaovAMF0 <- LSD.test(aovAMF0, c("Treatment"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovAMF0
aovAMF0 <- aov(AMF~Plant, data=TabCom0[36:48,])
summary(aovAMF0)
lsdaaovAMF0 <- LSD.test(aovAMF0, c("Plant"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovAMF0

aovAMF10 <- aov(AMF~Treatment*Season, data=TabCom10)
summary(aovAMF10)
aovAMF10 <- aov(AMF~Treatment, data=TabCom10[36:48,])
summary(aovAMF10)
lsdaaovAMF10 <- LSD.test(aovAMF10, c("Treatment"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovAMF10
aovAMF10 <- aov(AMF~Plant, data=TabCom10[36:48,])
summary(aovAMF10)
lsdaaovAMF10 <- LSD.test(aovAMF10, c("Plant"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovAMF10

#Saprophic_fungy
ggplot(data=TabPro0, aes(x=Sampling, y=Saprophic_fungy, fill=Plant)) +
  geom_bar(stat="identity", position='dodge')+
  theme_minimal()+ scale_fill_manual(values=c('Red','blue','green'))+ ggtitle("Saprophytic fungi")+ theme(
    plot.title = element_text(color="Black", size=14, face="bold.italic"))+ 
  xlab("Season")+ ylab("nmol / g soil")

ggplot(data=TabPro10, aes(x=Sampling, y=Saprophic_fungy, fill=Plant)) +
  geom_bar(stat="identity", position='dodge')+
  theme_minimal()+ scale_fill_manual(values=c('Red','blue','green'))+ ggtitle("Saprophytic fungi")+ theme(
    plot.title = element_text(color="Black", size=14, face="bold.italic"))+ 
  xlab("Season")+ ylab("nmol / g soil")

aovSaprophic_fungy0 <- aov(Saprophic_fungy~Treatment*Season, data=TabCom0)
summary(aovSaprophic_fungy0)
aovSaprophic_fungy0 <- aov(Saprophic_fungy~Treatment, data=TabCom0[36:48,])
summary(aovSaprophic_fungy0)
lsdaaovSaprophic_fungy0 <- LSD.test(aovSaprophic_fungy0, c("Treatment"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovSaprophic_fungy0
aovSaprophic_fungy0 <- aov(Saprophic_fungy~Plant, data=TabCom0[36:48,])
summary(aovSaprophic_fungy0)
lsdaaovSaprophic_fungy0 <- LSD.test(aovSaprophic_fungy0, c("Plant"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovSaprophic_fungy0

aovSaprophic_fungy10 <- aov(Saprophic_fungy~Treatment*Season, data=TabCom10)
summary(aovSaprophic_fungy10)
aovSaprophic_fungy10 <- aov(Saprophic_fungy~Treatment, data=TabCom10[36:48,])
summary(aovSaprophic_fungy10)
lsdaaovSaprophic_fungy10 <- LSD.test(aovSaprophic_fungy10, c("Treatment"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovSaprophic_fungy10
aovSaprophic_fungy10 <- aov(Saprophic_fungy~Plant, data=TabCom10[36:48,])
summary(aovSaprophic_fungy10)
lsdaaovSaprophic_fungy10 <- LSD.test(aovSaprophic_fungy10, c("Plant"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovSaprophic_fungy10

#Mesofauna
ggplot(data=TabPro0, aes(x=Sampling, y=Mesofauna, fill=Plant)) +
  geom_bar(stat="identity", position='dodge')+
  theme_minimal()+ scale_fill_manual(values=c('Red','blue','green'))+ ggtitle("Mesofauna")+ theme(
    plot.title = element_text(color="Black", size=14, face="bold.italic"))+ 
  xlab("Season")+ ylab("nmol / g soil")

ggplot(mapping = aes(x=Sampling, y=Mesofauna, fill=Plant)) +
  geom_bar(data = TabPro0, width = 0.8, stat = "identity", position='dodge') +
  geom_bar(data = TabPro10, width = 0.8, stat = "identity", fill = "white") +
  theme_classic() + scale_y_continuous(expand = c(0, 0))+ 
  scale_fill_manual(values=c('Red','blue','green'))

ggplot(data=TabPro10, aes(x=Sampling, y=Mesofauna, fill=Plant)) +
  geom_bar(stat="identity", position='dodge')+
  theme_minimal()+ scale_fill_manual(values=c('Red','blue','green'))+ ggtitle("Mesofauna")+ theme(
    plot.title = element_text(color="Black", size=14, face="bold.italic"))+ 
  xlab("Season")+ ylab("nmol / g soil")

aovMesofauna0 <- aov(Mesofauna~Treatment*Season, data=TabCom0)
summary(aovMesofauna0)
aovMesofauna0 <- aov(Mesofauna~Treatment, data=TabCom0[36:48,])
summary(aovMesofauna0)
lsdaaovMesofauna0 <- LSD.test(aovMesofauna0, c("Treatment"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovMesofauna0
aovMesofauna0 <- aov(Mesofauna~Plant, data=TabCom0[36:48,])
summary(aovMesofauna0)
lsdaaovMesofauna0 <- LSD.test(aovMesofauna0, c("Plant"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovMesofauna0

aovMesofauna10 <- aov(Mesofauna~Treatment*Season, data=TabCom10)
summary(aovMesofauna10)
aovMesofauna10 <- aov(Mesofauna~Treatment, data=TabCom10[36:48,])
summary(aovMesofauna10)
lsdaaovMesofauna10 <- LSD.test(aovMesofauna10, c("Treatment"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovMesofauna10
aovMesofauna10 <- aov(Mesofauna~Plant, data=TabCom10[36:48,])
summary(aovMesofauna10)
lsdaaovMesofauna10 <- LSD.test(aovMesofauna10, c("Plant"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovMesofauna10
head(TabPro)
#F/B
ggplot(data=TabPro0, aes(x=Sampling, y=F.B, fill=Plant)) +
  geom_bar(stat="identity", position='dodge')+
  theme_minimal()+ scale_fill_manual(values=c('Red','blue','green'))+ ggtitle("F.B")+ theme(
    plot.title = element_text(color="Black", size=14, face="bold.italic"))+ 
  xlab("Season")+ ylab("nmol / g soil")

ggplot(data=TabPro10, aes(x=Sampling, y=F.B, fill=Plant)) +
  geom_bar(stat="identity", position='dodge')+
  theme_minimal()+ scale_fill_manual(values=c('Red','blue','green'))+ ggtitle("F.B")+ theme(
    plot.title = element_text(color="Black", size=14, face="bold.italic"))+ 
  xlab("Season")+ ylab("nmol / g soil")

aovF.B0 <- aov(F.B~Treatment*Season, data=TabCom0)
summary(aovF.B0)
aovF.B0 <- aov(F.B~Treatment, data=TabCom0[36:48,])
summary(aovF.B0)
lsdaaovF.B0 <- LSD.test(aovF.B0, c("Treatment"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovF.B0
aovF.B0 <- aov(F.B~Plant, data=TabCom0[36:48,])
summary(aovF.B0)
lsdaaovF.B0 <- LSD.test(aovF.B0, c("Plant"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovF.B0

aovF.B10 <- aov(F.B~Treatment*Season, data=TabCom10)
summary(aovF.B10)
aovF.B10 <- aov(F.B~Treatment, data=TabCom10[36:48,])
summary(aovF.B10)
lsdaaovF.B10 <- LSD.test(aovF.B10, c("Treatment"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovF.B10
aovF.B10 <- aov(F.B~Plant, data=TabCom10[36:48,])
summary(aovF.B10)
lsdaaovF.B10 <- LSD.test(aovF.B10, c("Plant"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovF.B10

#Soil Moisture
ggplot(data=TabPro0, aes(x=Sampling, y=Hum, fill=Plant)) +
  geom_bar(stat="identity", position='dodge')+
  theme_minimal()+ scale_fill_manual(values=c('Red','blue','green'))+ ggtitle("F.B")+ theme(
    plot.title = element_text(color="Black", size=14, face="bold.italic"))+ 
  xlab("Season")+ ylab("nmol / g soil")

ggplot(data=TabPro10, aes(x=Sampling, y=F.B, fill=Plant)) +
  geom_bar(stat="identity", position='dodge')+
  theme_minimal()+ scale_fill_manual(values=c('Red','blue','green'))+ ggtitle("F.B")+ theme(
    plot.title = element_text(color="Black", size=14, face="bold.italic"))+ 
  xlab("Season")+ ylab("nmol / g soil")

aovHum0 <- aov(Hum~Treatment*Season, data=TabCom0)
summary(aovHum0)
aovF.B0 <- aov(F.B~Treatment, data=TabCom0[36:48,])
summary(aovF.B0)
lsdaaovF.B0 <- LSD.test(aovF.B0, c("Treatment"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovF.B0
aovHum0 <- aov(Hum~Plant, data=TabCom0[36:48,])
summary(aovHum0)
lsdaaovHum0 <- LSD.test(aovHum0, c("Plant"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovHum0

aovF.B10 <- aov(F.B~Treatment*Season, data=TabCom10)
summary(aovF.B10)
aovF.B10 <- aov(F.B~Treatment, data=TabCom10[36:48,])
summary(aovF.B10)
lsdaaovF.B10 <- LSD.test(aovF.B10, c("Treatment"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovF.B10
aovHum10 <- aov(Hum~Plant, data=TabCom10[36:48,])
summary(aovHum10)
lsdaaovHum10 <- LSD.test(aovHum10, c("Plant"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaaovHum10

#PCA
autoplot(prcomp(TabCom[1:36,6:19]))
prcomp(TabCom[1:12,17:26])
TabCom0$Both = factor(TabCom0$Both)
autoplot(prcomp(TabCom0[,17:26]), data = TabCom0, colour = "Plant", shape='Season', loadings = TRUE, loadings.colour = 'black',
         loadings.label = TRUE, loadings.label.size = 3, loadings.label.colour = "black",
         main='PCA 0-10 cm', frame = FALSE, pval=TRUE, frame.type = 'norm' )+ theme_few()
autoplot(prcomp(TabCom10[,16:24]), data = TabCom10, colour = "Plant", shape='Season',loadings = TRUE, loadings.colour = 'black',
         loadings.label = TRUE, loadings.label.size = 3, loadings.label.colour = "black",
         main='PCA 10-20 cm', frame = FALSE, frame.type= "convex")+ theme_few()
        
autoplot(prcomp(TabCom[,13:16]), data = TabCom, shape='Plant', colour = "Season", loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3,
         main='PCA Soil enzyme activity', frame = TRUE)+ theme_few() 
autoplot(prcomp(TabCom0[25:36,8:23]), data = TabCom, colour = "Plant", loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3,
         main='PCA Spring', frame = TRUE)+ theme_bw() 
autoplot(prcomp(TabCom0[37:48,8:23]), data = TabCom, colour = "Plant", loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3,
         main='PCA Summer', frame = TRUE)+ theme_bw() 

PCA05 = autoplot(prcomp(TabCom0[,13:26]), data = TabCom0, colour = "Plant", shape='Season',size=2.5, loadings = TRUE, loadings.colour = 'black',
         loadings.label = TRUE,, loadings.label.size = 5, loadings.label.colour = "black",
         main='PCA 0-10 cm', frame = FALSE, pval=TRUE, frame.type = 'norm' )+ theme_few() +  scale_y_continuous(limits = c(-0.2, 0.3))
 autoplot(prcomp(TabCom10[,13:26]), data = TabCom10, colour = "Plant", shape='Season',size=2.5,loadings = TRUE, loadings.colour = 'black',
         loadings.label = TRUE, loadings.label.size = 5, loadings.label.colour = "black",
         main='BG  NAG PME  AS', frame = FALSE, frame.type= "convex")+ theme_few()+  scale_y_continuous(limits = c(-0.3, 0.2))

ggarrange(PCA05, PCA510, ncol = 2, nrow=1, 
          labels = c("a", "b"))

#RDA
rda.Enz.FAME <- rda(TabCom[,17:26]~., TabCom[,13:16], scale=T)
rda.Enz.FAME0 <- rda(TabCom0[,17:26]~., TabCom0[,13:16], scale=T)
rda.Enz.FAME10 <- rda(TabCom10[,17:26]~., TabCom10[,13:16], scale=T)
summary(rda.Enz.FAME)
ordiplot(rda.Enz.FAME0, type = 'n', xlab="RDA1 (41.88%)", ylab="RDA2 (3.67%)",
         ylim = c(-3,3), xlim = c(-1,1))
ordiellipse(rda.Enz.FAME0, TabCom0$Both, display="sites", kind = "se", draw = "lines",show.groupslabels = as.character(TabCom0$Both), col = 1:12)
legend ('topright', col = 1:4, pch = 1, legend = levels (TabProPB$Treatment))
ef1 <- envfit (rda.Enz.FAME, TabCom[-5,16:25], choices = c(2,3), permutations = 0)
ef2 <- envfit (rda.Enz.FAME, TabCom[-5,10:15], choices = c(2,3), permutations = 0)
plot (ef1, col='red')
plot (ef2)
anova (rda.Enz.FAME10, by = 'margin', parallel = 4)

rda.qui.FAME <- rda(TabCom[,17:26]~., TabCom[,7:12], scale=T)
rda.qui.FAME0 <- rda(TabCom0[,17:26]~., TabCom0[,7:12], scale=T)
rda.qui.FAME10 <- rda(TabCom10[,17:26]~., TabCom10[,7:12], scale=T)
summary(rda.qui.FAME)
ordiplot(rda.qui.FAME, type = 'n', xlab="RDA1 (64.81%)", ylab="RDA2 (6.10%)",
         ylim = c(-3,3), xlim = c(-1,1))
ordiellipse(rda.Enz.FAME0, TabCom0$Both, display="sites", kind = "se", draw = "lines",show.groupslabels = as.character(TabCom0$Both), col = 1:12)
legend ('topright', col = 1:4, pch = 1, legend = levels (TabProPB$Treatment))
ef1 <- envfit (rda.qui.FAME, TabCom[-5,16:25], choices = c(2,3), permutations = 0)
ef2 <- envfit (rda.qui.FAME, TabCom[-5,6:11], choices = c(2,3), permutations = 0)
plot (ef1, col='red')
plot (ef2)
anova (rda.qui.FAME10, by = 'margin', parallel = 4)

rda.qui.Enz <- rda(TabCom[,13:16]~., TabCom[,7:12], scale=T)
rda.qui.Enz0 <- rda(TabCom0[,13:16]~., TabCom0[,7:12], scale=T)
rda.qui.Enz10 <- rda(TabCom10[,13:16]~., TabCom10[,7:12], scale=T)
summary(rda.qui.Enz)
ordiplot(rda.qui.Enz0, type = 'n', xlab="RDA1 (41.88%)", ylab="RDA2 (3.67%)",
         ylim = c(-3,3), xlim = c(-1,1))
ordiellipse(rda.Enz.FAME0, TabCom0$Both, display="sites", kind = "se", draw = "lines",show.groupslabels = as.character(TabCom0$Both), col = 1:12)
legend ('topright', col = 1:4, pch = 1, legend = levels (TabProPB$Treatment))
ef1 <- envfit (rda.qui.Enz, TabCom[-5,16:25], choices = c(2,3), permutations = 0)
ef2 <- envfit (rda.qui.Enz, TabCom[-5,10:15], choices = c(2,3), permutations = 0)
plot (ef1, col='red')
plot (ef2)
anova (rda.qui.Enz10, by = 'margin', parallel = 4)

rda.todo.FAME <- rda(TabCom[,17:26]~., TabCom[,7:16], scale=T)
rda.todo.FAME0 <- rda(TabCom0[,17:26]~., TabCom0[,7:16], scale=T)
rda.todo.FAME10 <- rda(TabCom10[,17:26]~., TabCom10[,7:16], scale=T)
summary(rda.todo.FAME)
ordiplot(rda.todo.FAME, type = 'n', xlab="RDA1 (64.81%)", ylab="RDA2 (6.10%)",
         ylim = c(-3,3), xlim = c(-1,1))
ordiellipse(rda.Enz.FAME0, TabCom0$Both, display="sites", kind = "se", draw = "lines",show.groupslabels = as.character(TabCom0$Both), col = 1:12)
legend ('topright', col = 1:4, pch = 1, legend = levels (TabProPB$Treatment))
ef1 <- envfit (rda.todo.FAME, TabCom[-5,16:25], choices = c(2,3), permutations = 0)
ef2 <- envfit (rda.todo.FAME, TabCom[-5,6:11], choices = c(2,3), permutations = 0)
plot (ef1, col='red')
plot (ef2)
anova (rda.todo.FAME10, by = 'margin', parallel = 4)


tablaCor<-cor(TabNor[,-(1:6)], method = 'pearson')
p.mat <-cor_pmat(TabNor[,-(1:6)], method = 'pearson')
ggcorrplot(tablaCor, p.mat = p.mat, hc.order = TRUE,
           type = 'lower', insig='blank')
write.table(tablaCor,"corr R BW en txt.txt")
write.table(p.mat, "corr P BW en txt.txt")


#SEM
myModel <- ' 
Bacteria ~ Hum  + pH + TOC+ Plant_Type
Fungi ~ Hum  + pH + TOC+ Plant_Type
PCA_enz ~ Bacteria + Hum + pH + Fungi + TOC + Plant_Type
TOC~Plant_Type + pH + Hum

Hum ~~Hum 
pH ~~ pH
Plant_Type~~Plant_Type
Bacteria~~Fungi

'

fit<- sem(model = myModel, data = TabNor)
summary(fit, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)
lavaanPlot(model=fit, node_options = list(shape = 'box', fontname = 'Helvetica'), 
           edge_options = list(color='black'), 
           stand = TRUE, coefs = TRUE, covs = TRUE, 
           stars = c('regress'))
parameterEstimates(fit)
fitmeasures(fit, c('df', 'rmsea', 'cfi', 'aic', 'bic', 'chisq', 'pvalue'))
head(TabNor)
varTable(TabNor)

data <- read.table("clipboard", header = TRUE)

g <- ggplot(data, aes(Month, TempAver))
g + geom_area() + 
  labs(title="Bar Chart", 
       subtitle="Mean Temperature", 
       caption="Source: Frequency of Manufacturers from 'mpg' dataset") +
  theme_few()
