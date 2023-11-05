packages <- c('ggthemes','dplyr', "ape", "ShortRead", "Biostrings", "phyloseq", "corncob","DESeq2", "microbiome", "DECIPHER", "phangorn", "tibble", "lme4", "lmerTest", "ggplot2", "vegan", "car", "rcompanion", "emmeans", "RVAideMemoire")
sapply(packages, require, character.only = TRUE)              

# Set working directory
setwd("C:/Users/juanp/Google Drive/labs/NCSU/Black Walnut project/Bacteria 16S/")#CC Network/

# Set plotting theme
theme_set(theme_bw())
# Assign variables for imported data
sharedfile = "BwBac.shared"
taxfile = "BwBac.taxonomy"

mothur_data <- import_mothur(mothur_shared_file = sharedfile,
                             mothur_constaxonomy_file = taxfile)

#Agregar la data
sampledata = sample_data(data.frame(
  Plant = c('BW', 'BS', 'BS', 'BS','BW', 'BW', 'BW', 'BW', 'BW', 'RO',
            'RO','BW', 'RO', 'RO','BS', 'BS', 'BS', 'BW','BW','BW','BW',
            'BW','BW','BW','RO','RO','RO','RO','BS', 'BS','RO', 'BS',
            'BW','BW','BW','BW','BW','RO','RO','RO','RO','RO','BS','BS','BS','RO','RO'),
  Subsamples = c(1,1,2,3,1 ,2,3,4,5,1, 2,2,3,4,1, 2,3,3,4,1, 5,2,3,4, 1,2,3,4,1, 2,1,3,1,2, 3,4,5,1,2, 3,4,2,1,2, 3,3,4),
  Season = c('Autumn','Autumn','Autumn','Autumn','Autumn','Autumn','Autumn','Autumn',
             'Autumn','Autumn','Autumn','Autumn','Autumn','Autumn','Autumn','Autumn',
             'Autumn','Autumn','Autumn','Spring','Autumn','Spring','Spring',
             'Spring','Spring','Spring','Spring','Spring','Spring','Spring','Autumn',
             'Spring','Spring','Spring','Spring','Spring','Spring','Spring','Spring',
             'Spring','Spring','Autumn','Spring','Spring','Spring','Autumn','Autumn'),
  row.names=sample_names(mothur_data),
  stringsAsFactors=FALSE,
  Depth = c('0-10','0-10','0-10','0-10','10-20','10-20','10-20','10-20','10-20',
            '10-20','10-20','0-10','10-20','10-20','10-20','10-20','10-20','0-10',
            '0-10','0-10','0-10','0-10','0-10','0-10','0-10','0-10',
            '0-10','0-10','0-10','0-10','0-10','0-10','10-20','10-20','10-20',
            '10-20','10-20','10-20','10-20','10-20','10-20','0-10','10-20',
            '10-20','10-20','0-10','0-10')))


#Phylum resumen
PhylumMean <- read.table("clipboard", header = TRUE)

#Metatabla
Data<- read.table("clipboard", header = TRUE)
nrow(Data)
meta <- Data
nrow(meta)

#Random phylogenetic tree
random_tree = rtree(ntaxa(mothur_data), rooted=TRUE, tip.label=taxa_names(mothur_data))

## Building the Phyloseq Object
physeq1 = merge_phyloseq(mothur_data, sampledata , random_tree)

# Now, let's make sure that the correct information is included in our phyloseq object
# we want to duplicate the phyloseq object, so that we can have 1 original and edit the other

physeq1.summary <- summary(physeq1@tax_table) # Should include asv info
physeq1@tax_table # Should include taxonomic info
physeq1@sam_data # Should reflect the mapping file that we imported
write.table(ps.T3@tax_table, "ps.T3 tax_table en txt.txt")

ps.A2 = subset_taxa(physeq1, Rank1 ==  "Archaea") #Kingdom == "Bacteria" | Kingdom == "Archaea"
ps.B2 = subset_taxa(physeq1, Rank1 ==  "Bacteria") #Kingdom == "Bacteria" | Kingdom == "Archaea"
ps.T2 = subset_taxa(physeq1, Rank1 ==  "Bacteria" |Rank1 ==  "Archaea") #Kingdom == "Bacteria" | Kingdom == "Archaea"
sample_sums(ps.2)
ps.pruned <- prune_samples(sample_sums(ps.T3)>=1000, ps.T3) 
ps.A3 = prune_taxa(taxa_sums(ps.A2) >= 5, ps.A2)
ps.B3 = prune_taxa(taxa_sums(ps.B2) >= 5, ps.B2)
ps.T3 = prune_taxa(taxa_sums(ps.T2) >= 5, ps.T2)
ps.perc <- transform_sample_counts(ps.T3, function(x) x / sum(x) * 100) 
ps.perc.T3 = prune_taxa(taxa_sums(ps.perc) >= 0.01, ps.perc)
ps.perc.T2 = prune_taxa(taxa_sums(ps.perc) >= 2, ps.perc)
ps.T4 = prune_taxa(taxa_sums(ps.T2) >= 1000, ps.T2)

###ANALYSIS###
# 1. NMDS plot####
plot_ordination(ps.A3, ordinate(ps.A3, "NMDS", "bray"), color = "Plant", shape = 'Depth') + theme_few() + guides(colour = guide_legend(override.aes = list(size = 8))) + guides(shape = guide_legend(override.aes = list(size = 8))) + geom_point(size = 4)
plot_ordination(ps.A3, ordinate(ps.A3, "NMDS", "bray"), color = "Plant", shape = 'Season') + theme_few() + guides(colour = guide_legend(override.aes = list(size = 8))) + guides(shape = guide_legend(override.aes = list(size = 8))) + geom_point(size = 4)+ geom_polygon(aes(fill = Plant))

plot_ordination(ps.B3, ordinate(ps.B3, "NMDS", "bray"), color = "Plant", shape = 'Depth') + theme_few() + guides(colour = guide_legend(override.aes = list(size = 8))) + guides(shape = guide_legend(override.aes = list(size = 8))) + geom_point(size = 4)
plot_ordination(ps.B3, ordinate(ps.B3, "NMDS", "bray"), color = "Plant", shape = 'Season') + theme_few() + guides(colour = guide_legend(override.aes = list(size = 8))) + guides(shape = guide_legend(override.aes = list(size = 8))) + geom_point(size = 4)

plot_ordination(ps.T3, ordinate(ps.T3, "NMDS", "bray"), color = "Plant", shape = 'Depth') + theme_few() + guides(colour = guide_legend(override.aes = list(size = 8))) + guides(shape = guide_legend(override.aes = list(size = 8))) + geom_point(size = 4)
plot_ordination(ps.T3, ordinate(ps.T3, "NMDS", "bray"), color = "Plant", shape = 'Season') + theme_few() + guides(colour = guide_legend(override.aes = list(size = 8))) + guides(shape = guide_legend(override.aes = list(size = 8))) + geom_point(size = 4)

# 2. PERMANOVA
#PERMANOVA
set.seed(1)

# Calculate bray curtis distance matrix
ps.A3_bray <- phyloseq::distance(ps.A3, method = "bray")
ps.B3_bray <- phyloseq::distance(ps.B3, method = "bray")
ps.T3_bray <- phyloseq::distance(ps.T3, method = "bray")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(physeq1))

# Adonis test
adonis(ps.A3_bray ~ Plant*Season*Depth, data = sampledf)
adonis(ps.B3_bray ~ Plant*Season*Depth, data = sampledf)
adonis(ps.T3_bray ~ Plant*Season*Depth, data = sampledf)

# if adonis is significant, move on to post-hoc test
pairwise.perm.manova(ps.A3_bray, sampledf$Plant, nperm = 999)
pairwise.perm.manova(ps.B3_bray, sampledf$Plant, nperm = 999)
pairwise.perm.manova(ps.T3_bray, sampledf$Plant, nperm = 999)


# 3. Beta-Diversity Estimates####
# there are two main things that you will likely want to analyze, a permanova and homogeneity of group dispersions 
# adonis runs a permutational MANOVA. One major assumption for this analysis is that within group similarity is similar across groups
# betadisper runs PERMDISP2 to calculate group dispersions
# we will start with permdisp and then run adonis, after creating a basic ordination plot

unifracA <- UniFrac(ps.A3, weighted=TRUE, normalized=TRUE)
unifracB <- UniFrac(ps.B3, weighted=TRUE, normalized=TRUE)
unifracT <- UniFrac(ps.T3, weighted=TRUE, normalized=TRUE)

# betadisper
disp.uniA <- betadisper(unifracA, sampledf$Plant)
anova(disp.uniA)

disp.uniB <- betadisper(unifracB, sampledf$Plant)
anova(disp.uniB)

disp.uniT <- betadisper(unifracT, sampledf$Plant)
anova(disp.uniT)

# if the anova had a p <0.05, perform post-hoc analysis like so:
TukeyHSD(disp.uniA) 
TukeyHSD(disp.uniB) 
TukeyHSD(disp.uniT) 

# creates confidence intervals between means of levels of a factor 
plot(disp.uniA, label=TRUE, main="Archaea")
plot(disp.uniB, label=TRUE, main="Bacteria")
plot(disp.uniT, label=TRUE, main="Bacteria + Archaea")

# what about phylogenetic distance?
adonis(unifracA ~ Plant*Season, data = sampledf)
pairwise.perm.manova(unifracA, sampledf$Plant, nperm= 999, p.method = "BH")

adonis(unifracB ~ Plant*Season, data = sampledf)
pairwise.perm.manova(unifracB, sampledf$Plant, nperm= 999, p.method = "BH")

adonis(unifracA ~ Plant*Season, data = sampledf)
pairwise.perm.manova(unifracT, sampledf$Plant, nperm= 999, p.method = "BH")

tax_table(ps.A3)
# 4. Plotting Relative Abundance Bar Charts####
# phylum-level
ps.compositional <- microbiome::transform(ps.perc, "compositional")
ps.phyla.perc <-taxa_level(ps.compositional, "Rank2")

# identify the 10 most abundant phylum
phylum.10 <- names(sort(taxa_sums(ps.phyla.perc), TRUE)[1:10])
# now we can use the list of the top phylum and subset those from the phyloseq object
ps.phylum.10 <- prune_taxa(phylum.10, ps.phyla.perc)

melt.phylum <- psmelt(ps.phylum.10)

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888",
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
#,"black", "black", "black", "black", "black","black", "black", "black", "black", "black","black", "black", "black", "black", "black")

ggplot(melt.phylum, aes(x = Sample, y = Abundance, fill = OTU)) + theme_few() +
  geom_bar(stat = "identity") + scale_fill_manual(values = safe_colorblind_palette) + labs(fill = "Phylum")

melt.phylum$Sample= factor(melt.phylum$Sample, levels=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,49,50,51,52,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72))
sample <- c('BW/Autumn/0-10',	'RO/Autumn/0-10',	'BS/Autumn/0-10','BW/Spring/0-10',	'RO/Spring/0-10',	'BS/Spring/0-10',	'BW/Autumn/10-20',	'RO/Autumn/10-20',	'BS/Autumn/10-20',	'BW/Spring/10-20',	'RO/Spring/10-20',	'BS/Spring/10-20')

###grafico 2
PhylumMean$Treatment = factor(PhylumMean$Treatment, levels=c('BW/Autumn/0-10',	'RO/Autumn/0-10',	'BS/Autumn/0-10','BW/Spring/0-10',	'RO/Spring/0-10',	'BS/Spring/0-10',	'BW/Autumn/10-20',	'RO/Autumn/10-20',	'BS/Autumn/10-20',	'BW/Spring/10-20',	'RO/Spring/10-20',	'BS/Spring/10-20'))
PhylumMean$Phylum = factor(PhylumMean$Phylum, levels=c("Other","Nitrospira","Chloroflexi","TM7","WS3","Armatimonadetes", "OD1", "Planctomycetes","Gemmatimonadetes", "Actinobacteria","Firmicutes","Verrucomicrobia", "Acidobacteria","Bacteroidetes", "Proteobacteria"))

ggplot(data=PhylumMean, aes(x=Treatment, y=Abundance, fill=Phylum)) + 
  geom_bar(stat="identity")+ theme_classic()+ ylab("Relative abundance (%)")+
  scale_fill_manual(values = safe_colorblind_palette) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5))

display.brewer.all()
display.brewer.pal(n = 16, name = 'Set1')

#5. Heatmap
lg10Genus <- read.table("clipboard", header = TRUE) 
Mean <- ddply(lg10Genus, "Treatment", colwise(mean)) 
Mean1 <- Mean
name <-Mean[,1]
Treatment <- as.factor(Mean1$Treatment)

Mean0_15<- data.frame(t(Mean1[-1]))
names (Mean0_15) = c('BW/Autumn/0-10','RO/Autumn/0-10','BS/Autumn/0-10',
                     'BW/Autumn/10-20','RO/Autumn/10-20',"BS/Autumn/10-20",
                     'BW/Spring/0-10','RO/Spring/0-10','BS/Spring/0-10',
                     'BW/Spring/10-20','RO/Spring/10-20','BS/Spring/10-20')
annotation_col = data.frame(
  Genera =substr(colnames(Mean0_15),1,3))
rownames(annotation_col)=colnames(Mean0_15)
colnames(Mean0_15)
row.names(Mean1$ID)



pheatmap(Mean0_15, show_rownames=TRUE,show_colnames=TRUE,
         annotation_col=annotation_col, annotation_legend = TRUE,
         scale = "none",clustering_method="ward.D2",
         fontsize_row = 6,
         clustering_distance_cols="euclidean")


##6.  Mantel Tests ####
asv.table <- data.frame(otu_table(ps.B4))

xdist.prod <- vegdist(asv.table, method = "bray")
ydist.prod <- vegdist(Metatabla$TN, method = "euclid")

# mantel test with bray-curtis dissimilarity
vegan::mantel(xdist.prod, ydist.prod, method = "spear")
# mantel test with weighted UniFrac distance
vegan::mantel(unifracB, ydist.prod, method = "spear")
# we can also include more than 1 environmental variable
ydist.prod <- vegdist(Metatabla$TN*Metatabla$pH, method = "euclid")
vegan::mantel(unifracT, ydist.prod, method = "spear")



##7.  Alfadiversity
ordered(sample_sums(ps.T3))

# let's calcuolate Shannon diversity and Richness (Observed)
alpha.div<-estimate_richness(ps.T3, measures=c("Shannon", "Observed",'Chao', 'Simpson','Richness'))
even <- evenness(ps.T3, 'pielou')
write.table(Data, "Alfadiversity en txt.txt")

# it is easiest to append the alpha-diversity estimates to the metadata file for downstream analyses
Data$Shannon <- paste(alpha.div$Shannon)
Data$Richness <- paste(alpha.div$Richness)
Data$Observed <- paste(alpha.div$Observed)
Data$Simpson <- paste(alpha.div$Simpson)
Data$Chao1 <- paste(alpha.div$Chao1)
Data$Evenness <- even$pielou
Data$Observed <- as.numeric(Data$Observed)
Data$Shannon <- as.numeric(Data$Shannon)
Data$Simpson <- as.numeric(Data$Simpson)
Data$Chao1 <- as.numeric(Data$Chao1)
Data$Plant <- as.factor(Data$Plant)
Data$Depth <- as.factor(Data$Depth)
Data$Season <- as.factor(Data$Season)
Data$Samples <- as.factor(Data$Samples)

##Shannon
# the below commands will perform anova's based on flavonoid production and plant genotype
# Note that I have included Block as a random factor here, thus I have to make a linear mixed effects model before running the anova
# make sure lmerTest is a loaded package or you will not see p-values
shannon.model<-lmer(Shannon ~ Plant*Depth*Season   + (1|Samples), data = Data)

# check assumptions for linear modeling
# 1. linearity
# plot fitted versus residuals
plot(shannon.model)
# plot the observed values versus the model values
# should see no pattern
# 2. homogeneity of variance
leveneTest(residuals(shannon.model) ~ Data$Plant*Data$Depth)
# if p <0.05, assumption is not met

# 3. normal distribution of residuals
qqmath(shannon.model)
# values should ~ follow the line

# assuming that the assumptions are satisfied, run the anova

anova(shannon.model) 
anova(shannon.model2) 
# note that this is a Type III Analysis of variance table

ggplot(Data, aes(Shannon, Season , fill = Plant)) + geom_boxplot() + theme_classic()

# post-hoc test using estimated marginal means
emmeans(shannon.model, pairwise ~ Plant, adjust = "none")
# indicates no differences between lines

#Chao
Chao.model<-lmer(Chao1 ~ Plant*Depth*Season   + (1|Samples), data = Data)

# check assumptions for linear modeling
# 1. linearity
# plot fitted versus residuals
plot(Chao.model)
# plot the observed values versus the model values
# should see no pattern
# 2. homogeneity of variance
leveneTest(residuals(Chao.model) ~ Data$Plant*Data$Depth)
# if p <0.05, assumption is not met

# 3. normal distribution of residuals
qqmath(Chao.model)
# values should ~ follow the line

# assuming that the assumptions are satisfied, run the anova

anova(Chao.model) 
# note that this is a Type III Analysis of variance table

ggplot(Data, aes(Chao1, Season , fill = Plant)) + geom_boxplot() + theme_classic()

# post-hoc test using estimated marginal means
emmeans(Chao.model, pairwise ~ Plant, adjust = "none")
emmeans(Chao.model, pairwise ~ Plant*Season*Depth, adjust = "none")
# indicates no differences between lines

#Observed
Observed.model<-lmer(Observed ~ Plant*Depth*Season   + (1|Samples), data = Data)

# check assumptions for linear modeling
# 1. linearity
# plot fitted versus residuals
plot(Observed.model)
# plot the observed values versus the model values
# should see no pattern
# 2. homogeneity of variance
leveneTest(residuals(Observed.model) ~ Data$Plant)
# if p <0.05, assumption is not met

# 3. normal distribution of residuals
qqmath(Observed.model)
# values should ~ follow the line

# assuming that the assumptions are satisfied, run the anova

anova(Observed.model) 
# note that this is a Type III Analysis of variance table

ggplot(Data, aes(Observed, Plant , fill = Plant)) + geom_boxplot() + theme_classic()

# post-hoc test using estimated marginal means
emmeans(Chao.model, pairwise ~ Plant, adjust = "none")
# indicates no differences between lines

### 7. Differential Abundance Analysis DESeq2 ####
ps.T3RB <- subset_samples(ps.T4, Plant == "RO" | Plant == "BW")
diagdds = phyloseq_to_deseq2(ps.T4, ~ Depth)
diagdds = DESeq(diagdds, test= "Wald" , fitType = "parametric")

#Investigate test reuslts table
res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.T3RB)[rownames(sigtab),], "matrix"))
head(sigtab)
dim(sigtab)


##8. DESeq2 Plot
scale_fill_discrete <- function(palname = "Set1",...){
  scale_fill_brewer(palette = palname, ...)
}
#Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Rank2, function(x) max(x))
x = sort(x, TRUE)
sigtab$Rank2 = factor(as.character(sigtab$Rank2), levels = names(x))
#Genus Order
x = tapply(sigtab$log2FoldChange, sigtab$Rank6, function(x) max(x))
x = sort(x, TRUE)
sigtab$Rank6 = factor(as.character(sigtab$Rank6), levels = names(x))

ggplot(sigtab, aes(x=Rank6, y=log2FoldChange, color=Rank2)) + geom_point(size=5) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) + 
  ggtitle("Plot of 10-20 vs. 0-10") 

## Indicator specie by multipatt (indicspecies)
# Write about your phyloseq OTY tabe and export it
write.csv(ps.T4@otu_table, 'ps.T4 otu_table.cvs')
write.csv(ps.perc.T2@otu_table, 'ps.perc.T2 otu_table.cvs')
write.csv(ps.perc.T2@tax_table, 'ps.perc.T2 tax_table.cvs')


# Import phyloseq OTU table as ab OTU table/dataframe
SpOTU <- read.csv('ps.T4 otu_table.cvs')
metadata <- read.csv('metadata.cvs')

#do some shuffling of the OTU table
SpOTUFlip <- as.data.frame(t(SpOTU)) #make it a dataframe and puts x into y and into x (flips it)
names(SpOTUFlip) <- as.matrix(SpOTUFlip[1, ]) #rename columns
SpOTUFlip <- SpOTUFlip[-1, ] #remove first row
SpOTUFlip_num <- as.data.frame(lapply(SpOTUFlip, as.numeric)) #convert from character to number
SpOTUFlip_num$SampleID <- row.names(SpOTUFlip) #puts row names as sample ID column

head(SpOTUFlip)
#read metadata
metadata <- meta[,1:5]
ncol(metadata)
nrow(metadata)

metadata$ID  <- as.factor(metadata$ID)
metadata$Plant  <- as.factor(metadata$Plant)

#Join based on SampleID
SpOTU_Final <- left_join(SpOTUFlip_num, metadata, by = c("SampleID" = "ID")) 
ncol(SpOTU_Final)
nrow(SpOTU_Final)
ps.T4

SPotus = SpOTU_Final[,1:5000] #Select just the OTU table part of the file
SPPlant = SpOTU_Final$Plant #The metadata column gruop you care about
SPDepth = SpOTU_Final$Depth #The metadata column gruop you care about
SPSeason = SpOTU_Final$Season #The metadata column gruop you care about

SPindS=multipatt(x=SPotus, cluster= metadata$Depth, func = "r.g", control = how(nperm=9999))
summary(SPindS)

###Spearman correlation
tablaCorBac<-cor(Data[,-(1:5)], method = 'spearman')
p.matBac <-cor_pmat(Data[,-(1:5)], method = 'spearman')
ggcorrplot(tablaCorBac, p.mat = p.matBac, hc.order = TRUE,
           type = 'lower', insig='blank')
write.table(tablaCorBac,"corr R Bac en txt.txt")
write.table(p.matBac, "corr P Bac en txt.txt")

### Soil properties analysis ANOVA 
aovTOC <- aov(TOC~Plant*Season*Depth, data=Data)
summary(aovTOC)
lsdaovTOC1 <- LSD.test(aovTOC, c("Plant"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaovTOC2 <- LSD.test(aovTOC, c("Season"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaovTOC3 <- LSD.test(aovTOC, c("Depth"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 

aovTC <- aov(TC~Plant*Season*Depth, data=Data)
summary(aovTC)
lsdaovTC1 <- LSD.test(aovTC, c("Plant"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaovTC2 <- LSD.test(aovTC, c("Season"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaovTC3 <- LSD.test(aovTC, c("Depth"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 

aovTN <- aov(TN~Plant*Season*Depth, data=Data)
summary(aovTN)
lsdaovTN1 <- LSD.test(aovTN, c("Plant"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaovTN2 <- LSD.test(aovTN, c("Season"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaovTN3 <- LSD.test(aovTN, c("Depth"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 

aovC.N <- aov(C.N~Plant*Season*Depth, data=Data)
summary(aovC.N)
lsdaovC.N1 <- LSD.test(aovC.N, c("Plant"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaovC.N2 <- LSD.test(aovC.N, c("Season"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaovC.N3 <- LSD.test(aovC.N, c("Depth"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 

aovHum <- aov(Hum~Plant*Season*Depth, data=Data)
summary(aovHum)
lsdaovHum1 <- LSD.test(aovHum, c("Plant"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaovHum2 <- LSD.test(aovHum, c("Season"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaovHum3 <- LSD.test(aovHum, c("Depth"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 

aovpH <- aov(pH~Plant*Season*Depth, data=Data)
summary(aovpH)
lsdaovpH1 <- LSD.test(aovpH, c("Plant"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaovpH2 <- LSD.test(aovpH, c("Season"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdaovpH3 <- LSD.test(aovpH, c("Depth"), alpha = 0.05, p.adj="bonferroni", group=TRUE) 

##Venn Diagram
ps.perc.T3.BW <- subset_samples(ps.perc, Plant == "BW")
ps.perc.T3.BS <- subset_samples(ps.perc, Plant == "BS")
ps.perc.T3.RO <- subset_samples(ps.perc, Plant == "RO")

ps.perc.T3.RO1 = prune_taxa(taxa_sums(ps.perc.T3.RO) > 0, ps.perc.T3.RO)
ps.perc.T3.BW1 = prune_taxa(taxa_sums(ps.perc.T3.BW) > 0, ps.perc.T3.BW)
ps.perc.T3.BS1 = prune_taxa(taxa_sums(ps.perc.T3.BS) > 0, ps.perc.T3.BS)

Bw <- rownames(otu_table(ps.perc.T3.BW1))
BS <- rownames(otu_table(ps.perc.T3.BS1))
RO <- rownames(otu_table(ps.perc.T3.RO1))

Venn <-venn(list(A=RO, B=Bw, c=BS))
venn(list(A=RO, B=Bw, c=BS))
