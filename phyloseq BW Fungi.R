packages <- c('ggcorrplot','ggthemes','dplyr', "ape", "ShortRead", "Biostrings", 
              "phyloseq", "corncob","DESeq2", "microbiome", "DECIPHER", "phangorn",
              "tibble", "lme4", "lmerTest", "ggplot2", "vegan", "car", "rcompanion", 'microbiomeSeq',
              "emmeans", "RVAideMemoire",'gplots','plotly','tidyr','VennDiagram','venneuler')
sapply(packages, require, character.only = TRUE)              


# Set working directory
setwd("C:/Users/juanp/Google Drive/labs/NCSU/Black Walnut project/ITS/")
  # Source code files

# Set plotting theme
theme_set(theme_bw())
# Assign variables for imported data
sharedfile = "BwFun.shared"
taxfile = "BwFun.taxonomy"

# Import mothur data
mothur_data <- import_mothur(mothur_shared_file = sharedfile,
                             mothur_constaxonomy_file = taxfile)



#Agregar la data
sampledata = sample_data(data.frame(
  Plant = c('BW', 'BS', 'BS', 'BS','BW', 'BW', 'BW', 'BW', 'BW', 'RO',
              'RO','BW', 'RO', 'RO','BS', 'BS', 'BS', 'BW','BW','BW','BW',
            'BW','BW','BW','BW','RO','RO','RO','RO','BS', 'BS','RO', 'BS',
            'BW','BW','BW','BW','BW','RO','RO','RO','RO','RO','BS','BS','BS','RO','RO'),
    Subsamples = c(1,1,2,3,1,2,3,4,5, 1,2,2,3,4,1,2,3,3,4,1,5,2,3,4,5,1,2,3,4,1,2,1,3,1,2,3,4,5,1,2,3,4,2,1,2,3,3,4),
  Season = c('Autumn','Autumn','Autumn','Autumn','Autumn','Autumn','Autumn','Autumn',
             'Autumn','Autumn','Autumn','Autumn','Autumn','Autumn','Autumn','Autumn',
             'Autumn','Autumn','Autumn','Spring','Autumn','Spring','Spring','Spring',
             'Spring','Spring','Spring','Spring','Spring','Spring','Spring','Autumn',
             'Spring','Spring','Spring','Spring','Spring','Spring','Spring','Spring',
             'Spring','Spring','Autumn','Spring','Spring','Spring','Autumn','Autumn'),
  row.names=sample_names(mothur_data),
  stringsAsFactors=FALSE,
Depth = c('0-10','0-10','0-10','0-10','10-20','10-20','10-20','10-20','10-20',
          '10-20','10-20','0-10','10-20','10-20','10-20','10-20','10-20','0-10',
          '0-10','0-10','0-10','0-10','0-10','0-10','0-10','0-10','0-10',
          '0-10','0-10','0-10','0-10','0-10','0-10','10-20','10-20','10-20',
          '10-20','10-20','10-20','10-20','10-20','10-20','0-10','10-20',
          '10-20','10-20','0-10','0-10')))

meta <- read.table("Metadata en txt.txt", header = TRUE, row.names = 1)

#Random phylogenetic tree
random_tree = rtree(ntaxa(mothur_data), rooted=TRUE, tip.label=taxa_names(mothur_data))
plot(random_tree)

physeq1 = merge_phyloseq(mothur_data, sampledata)

# Now, let's make sure that the correct information is included in our phyloseq object
# we want to duplicate the phyloseq object, so that we can have 1 original and edit the other

physeq1.summary <- summary(physeq1@otu_table) # Should include asv info
physeq1@tax_table # Should include taxonomic info
physeq1@sam_data # Should reflect the mapping file that we imported
write.table(ps.3@tax_table, "ps.3 tax_table en txt.txt")

ps.2 = subset_taxa(physeq1, Rank1 ==  "Fungi") #
ps.perc2 <- transform_sample_counts(ps.3, function(x) x / sum(x) * 100) 
sample_sums(ps.2)
ps.pruned <- prune_samples(sample_sums(ps.2)>=00, ps.2) 
ps.3 = prune_taxa(taxa_sums(ps.2) >= 5, ps.2)
ps.4 = prune_taxa(taxa_sums(ps.2) >= 100, ps.2)
ps.perc2 = prune_taxa(taxa_sums(ps.perc2) >= 0.01, ps.perc2)

# 1. NMDS plot####
plot_ordination(ps.3, ordinate(ps.3, "NMDS", "bray"), color = "Plant", shape = 'Depth') + theme_few() + guides(colour = guide_legend(override.aes = list(size = 8))) + guides(shape = guide_legend(override.aes = list(size = 8))) + geom_point(size = 4)
plot_ordination(ps.3, ordinate(ps.3, "NMDS", "bray"), color = "Plant", shape = 'Season') + theme_few() + guides(colour = guide_legend(override.aes = list(size = 8))) + guides(shape = guide_legend(override.aes = list(size = 8))) + geom_point(size = 4)
 plot_ordination(ps.3, ordinate(ps.3, "NMDS", "bray"), color = 'Season') + theme_few() + guides(colour = guide_legend(override.aes = list(size = 8))) + guides(shape = guide_legend(override.aes = list(size = 8))) + geom_point(size = 4)

pca_res <- ordinate(ps.3, "NMDS", "bray")
fit2 <- envfit(pca_res ~ TOC + TN, metadata, perm = 199)

plot_ordination(pca_res, geom = 'label_repel')
plot(fit2)


# 18. Beta-Diversity Estimates####
# there are two main things that you will likely want to analyze, a permanova and homogeneity of group dispersions 
# adonis runs a permutational MANOVA. One major assumption for this analysis is that within group similarity is similar across groups
# betadisper runs PERMDISP2 to calculate group dispersions
# we will start with permdisp and then run adonis, after creating a basic ordination plot


unifrac <- UniFrac(ps.3, weighted=TRUE, normalized=TRUE)

# betadisper
disp.uni <- betadisper(unifrac, DataITS$Plant)
anova(disp.uni)

# if the anova had a p <0.05, perform post-hoc analysis like so:
TukeyHSD(disp.uni) 

# creates confidence intervals between means of levels of a factor 
plot(disp.uni, label=TRUE, main="Fungi")


# 3. Plotting Relative Abundance Bar Charts####
# phylum-level
ps.compositional <- microbiome::transform(ps.perc2, "compositional")
ps.phyla.perc <-taxa_level(ps.perc2, "Rank2")

# identify the 10 most abundant phylum
phylum.10 <- names(sort(taxa_sums(ps.phyla.perc), TRUE)[1:10])
# now we can use the list of the top phylum and subset those from the phyloseq object
ps.phylum.10 <- prune_taxa(phylum.10, ps.phyla.perc)

melt.phylum <- psmelt(ps.phylum.10)

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
#,"black", "black", "black", "black", "black","black", "black", "black", "black", "black","black", "black", "black", "black", "black")

melt.phylum$Sample= factor(melt.phylum$Sample, levels=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72))

ggplot(data=melt.phylum, aes(x=Sample, y=Abundance, fill=OTU)) + 
  geom_bar(stat="identity")+ theme_classic()+ ylab("Relative abundance (%)")+
  scale_fill_manual(values=brewer.pal(n = 10, name = "Paired"))

# do you notice anything ...unpleasant about this plot? - relative abundances are used and should sum to 1...

phyl.means <- aggregate(Abundance~Plant, melt.phylum, FUN=mean)

ggplot(phyl.means, aes(x = Plant, y = Abundance, fill = Season)) + theme_classic() + #, fill = OTU
  geom_bar(stat = "identity") + scale_fill_manual(values = safe_colorblind_palette)

# bar plot at the genus-level
ps.gen <- taxa_level(ps.3, "Rank6")

gen.10 <- names(sort(taxa_sums(ps.gen), TRUE)[1:15])
ps.gen.10 <- prune_taxa(gen.10, ps.gen)
melt.gen <- psmelt(ps.gen.10)

gen.means <- aggregate(Abundance~Plant+OTU, melt.gen, FUN=mean)

ggplot(gen.means, aes(x = Plant, y = Abundance, fill = OTU))  +
  geom_bar(stat = "identity") + theme_classic()

# is there a difference in abundance of Pseudomonas between the Lines?
# extract Pseudomonas abundances from melt.gen file
pseud <- melt.gen[ which(melt.gen$OTU=='Pseudomonas'), ]

# create model
pseud.lm <- lm(Abundance ~ Season, data = pseud)
anova(pseud.lm)

emmeans(pseud.lm, pairwise ~ Line)

# 20. Mantel Tests ####
meta <- read.table("Metadata en txt.txt", header = TRUE, row.names = 1)
asv.table <- data.frame(otu_table(ps.3))

xdist.prod <- vegdist(asv.table, method = "bray")
ydist.prod <- vegdist(meta$AS, method = "euclid")

# mantel test with bray-curtis dissimilarity
vegan::mantel(xdist.prod, ydist.prod, method = "spear")
# mantel test with weighted UniFrac distance
vegan::mantel(unifrac, ydist.prod, method = "spear")
# we can also include more than 1 environmental variable
ydist.prod <- vegdist(meta$TOC*meta$TN, method = "euclid")
vegan::mantel(unifrac, ydist.prod, method = "spear")


# 21. Differential Abundance Analysis Corn Cob ####
# differential abundance analyses compare 2 groups 
# Since we have three groups, we need to split the data 
set.seed(1)

ps.CH <- subset_samples(ps.2, Season == "Autumn")
ps.CH2 <- subset_samples(ps.CH, Plant == "RO" | Plant == "BW")
ps.RN <- subset_samples(ps.pruned, Treatment == "Agriculture" | Treatment == "Grassland")
ps.LV <- subset_samples(ps.pruned, Treatment == "Agriculture" | Treatment == "Grassland")

ex1 <- differentialTest(formula = ~ Plant,
                        phi.formula ~ Plant,
                        formula_null = ~ 1,
                        phi.formula_null = ~ Plant,
                        data = ps.CH,
                        test = "wald", boot = FALSE,
                        fdr_cutoff = 0.05)
plot(ex1)

otu_table(ps.3)
otu_to_taxonomy(OTU = "Otu000001", data = ibd)
check_GN04 <- bbdml(formula = Otu000001 ~ ibd,
                    phi.formula = ~ ibd,
                    data = ps.CH)
lrtest(mod_null = check_GN04, mod = check_GN04)
plot(check_GN04, B = 50, color = "Plant")
summary(check_GN04)


ex2 <- differentialTest(formula = ~ Biochar,
                        phi.formula = ~ Biochar,
                        formula_null = ~ 1,
                        phi.formula_null = ~ Biochar,
                        data = ps.N321.N336,
                        test = "Wald", boot = FALSE,
                        fdr_cutoff = 0.05)

plot(ex2)

ex3 <- differentialTest(formula = ~ Biochar,
                        phi.formula = ~ Biochar,
                        formula_null = ~ 1,
                        phi.formula_null = ~ Biochar,
                        data = ps.N334.N336,
                        test = "Wald", boot = FALSE,
                        fdr_cutoff = 0.05)
plot(ex3)

# how do you know which Line has more of ASV107?

plot(asv.table$ASV107 ~ meta$Line)
aggregate(asv.table$ASV107~Line, meta, FUN=mean)

#Alfadiversity
ordered(sample_sums(ps.3))

# let's calcuolate Shannon diversity and Richness (Observed)
alpha.div<-estimate_richness(ps.3, measures=c("Shannon", "Observed",'Chao', 'Simpson','Richness'))
even <- evenness(ps.3, 'pielou')
write.table(alpha.div, "Alfadiversity en txt.txt")


# it is easiest to append the alpha-diversity estimates to the metadata file for downstream analyses
DataITS$Shannon <- paste(alpha.div$Shannon)
DataITS$Richness <- paste(alpha.div$Richness)
DataITS$Observed <- paste(alpha.div$Observed)
DataITS$Simpson <- paste(alpha.div$Simpson)
DataITS$Chao1 <- paste(alpha.div$Chao1)
DataITS$Evenness <- even$pielou
DataITS$Observed <- as.numeric(DataITS$Observed)
DataITS$Shannon <- as.numeric(DataITS$Shannon)
DataITS$Simpson <- as.numeric(DataITS$Simpson)
DataITS$Chao <- as.numeric(DataITS$Chao)
DataITS$Plant <- as.factor(DataITS$Plant)
DataITS$Sample <- as.factor(DataITS$Sample)
DataITS$Depth <- as.factor(DataITS$Depth)
DataITS$Season <- as.factor(DataITS$Season)

# the below commands will perform anova's based on flavonoid production and plant genotype
# Note that I have included Block as a random factor here, thus I have to make a linear mixed effects model before running the anova
# make sure lmerTest is a loaded package or you will not see p-values
shannon.model<-lmer(Shannon ~ Plant*Season*Depth  + (1|Samples), data = DataITS)
shannon.model2<-lmer(Shannon ~ Treatment  + (1|Sample), data = Data)

# check assumptions for linear modeling
# 1. linearity
# plot fitted versus residuals
plot(shannon.model)
# plot the observed values versus the model values
# should see no pattern
# 2. homogeneity of variance
leveneTest(residuals(shannon.model) ~ DataITS$Plant)
# if p <0.05, assumption is not met

# 3. normal distribution of residuals
qqmath(shannon.model)
# values should ~ follow the line

# assuming that the assumptions are satisfied, run the anova

anova(shannon.model) 
# note that this is a Type III Analysis of variance table

ggplot(DataITS, aes(Shannon, Plant, fill = Plant)) + geom_boxplot() + theme_classic()


###Observed
Observed.model<-lmer(Observed ~ Plant*Season*Depth  + (1|Samples), data = DataITS)

ps.3a <- subset_samples(ps.3, Plant == "RO" | Plant == "BW")
ps.3b <- subset_samples(ps.3a, Subsamples != "5")
ps.3c <- subset_samples(ps.3b, Season == "Spring")

plot_anova_diversity(ps.3b, method = c("richness", "shannon", "simpson"), grouping_column = "Plant")
ps.3b@sam_data

# check assumptions for linear modeling
# 1. linearity
# plot fitted versus residuals
plot(Observed.model)
# plot the observed values versus the model values
# should see no pattern
# 2. homogeneity of variance
leveneTest(residuals(Observed.model) ~ DataITS$Plant)
# if p <0.05, assumption is not met

# 3. normal distribution of residuals
qqmath(Observed.model)
# values should ~ follow the line

# assuming that the assumptions are satisfied, run the anova

anova(Observed.model) 
# note that this is a Type III Analysis of variance table

ggplot(DataITS, aes(Observed, Plant, fill = Plant)) + geom_boxplot() + theme_classic()

# post-hoc test using estimated marginal means
emmeans(Observed.model, pairwise ~ Plant, adjust = "none")
# indicates no differences between lines

# if covariate, then we can use estimated marginal means of linear trends
emtrends(even.model, pairwise ~ Season, var = "Plant")
# indicates that Lines do not have a different trend of Total Phenolics ~ Richness 
# we can plot this like so:
emmip(even.model, Line ~ Total_Phenolics, cov.reduce = range)

####PERMANOVA
set.seed(1)

# Calculate bray curtis distance matrix
ps.2_bray <- phyloseq::distance(ps.3, method = "bray")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(ps.3))

# Adonis test
adonis(ps.2_bray ~ Plant*Season*Depth, data = sampledf)
adonis(ps.2_bray ~ Plant, data = sampledf)
adonis(ps.2_bray ~ Season, data = sampledf)
adonis(ps.2_bray ~ Depth, data = sampledf)

# adonis
adonis(bc.phyl ~ Plant*Depth*Season, data = DataITS)

adonis(bc.asv ~ Plant, data = DataITS)
# if adonis is significant, move on to post-hoc test
pairwise.perm.manova(bc.asv, DataITS$Plant, nperm = 999)
# what about phylogenetic distance?
adonis(unifrac ~ Plant, data = DataITS)
pairwise.perm.manova(unifrac, DataITS$Plant, nperm= 999, p.method = "BH")

ps.3Au <- subset_samples(ps.3, Season == "Autumn")
ps.3Sp <- subset_samples(ps.3, Season == "Spring")
ps.3Au_bray <- phyloseq::distance(ps.3Au, method = "bray")
sampledf <- data.frame(sample_data(ps.3Au))
adonis(ps.3Au_bray ~ Plant, data = sampledf)

ps.3Sp_bray <- phyloseq::distance(ps.3Sp, method = "bray")
sampledf <- data.frame(sample_data(ps.3Sp))
adonis(ps.3Au_bray~Plant, data = sampledf)

### 22. Differential Abundance Analysis DESeq2 ####
ps.3a <- subset_samples(ps.4, Plant == "RO" | Plant == "BW")
diagdds = phyloseq_to_deseq2(ps.4, ~ Plant)
diagdds = DESeq(diagdds, test= "Wald" , fitType = "parametric")

#Investigate test reuslts table
res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.3a)[rownames(sigtab),], "matrix"))
  head(sigtab)
dim(sigtab)

##DESeq2 Plot
scale_fill_discrete <- function(palname = "Set1",...){
  scale_fill_brewer(palette = palname, ...)
}
#Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Rank2, function(x) max(x))
x = sort(x, TRUE)
sigtab$Rank2 = factor(as.character(sigtab$Rank2), levels = names(x))
#Genus Order
x = tapply(sigtab$log2FoldChange, sigtab$Rank7, function(x) max(x))
x = sort(x, TRUE)
sigtab$Rank7 = factor(as.character(sigtab$Rank7), levels = names(x))

ggplot(sigtab, aes(x=Rank7, y=log2FoldChange, color=Rank2)) + geom_point(size=5) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) + 
  ggtitle("Plot of RO vs. BW") 


## Indicator specie by multipatt (Interspecie)
# Write about your phyloseq OTY tabe and export it
write.csv(ps.3@otu_table, 'ps.3 otu_table.cvs')
write.csv(ps.3@tax_table, 'ps.3 tax_table.cvs')
write.csv(metadata, 'metadata.cvs')


# Import phyloseq OTU table as ab OTU table/dataframe
SpOTU <- read.csv('ps.3 otu_table.cvs')
metadata <- read.csv('metadata.cvs')

#do some shuffling of the OTU table
SpOTUFlip <- as.data.frame(t(SpOTU)) #make it a dataframe and puts x into y and into x (flips it)
names(SpOTUFlip) <- as.matrix(SpOTUFlip[1, ]) #rename columns
SpOTUFlip <- SpOTUFlip[-1, ] #remove first row
SpOTUFlip_num <- as.data.frame(lapply(SpOTUFlip, as.numeric)) #convert from character to number
SpOTUFlip_num$SampleID <- row.names(SpOTUFlip) #puts row names as sample ID column

head(SpOTUFlip)
#read metadata
metadata <- read.table("clipboard", header = TRUE)
metadata$ID  <- as.factor(metadata$ID)
metadata$Plant  <- as.factor(metadata$Plant)

#Join based on SampleID
SpOTU_Final <- left_join(SpOTUFlip_num, metadata, by = c("SampleID" = "ID")) 
SpOTU_Final <- cbind(SpOTUFlip_num, metadata[,1])

SPotus = SpOTU_Final[,1:7897] #Select just the OTU table part of the file
SPPlant = SpOTU_Final$Plant #The metadata column gruop you care about
SPDepth = SpOTU_Final$Depth #The metadata column gruop you care about
SPSeason = SpOTU_Final$Season #The metadata column gruop you care about

SPindFD=multipatt(x=SPotus, cluster= metadata$Depth, func = "r.g", control = how(nperm=999))
summary(SPindFD)

###Spearman correlation
tablaCor<-cor(metadata[,-(1:5)], method = 'spearman')
p.mat <-cor_pmat(metadata[,-(1:5)], method = 'spearman')
ggcorrplot(tablaCor, p.mat = p.mat, hc.order = TRUE,
           type = 'lower', insig='blank')
write.table(tablaCor,"corr R BW en txt.txt")
write.table(p.mat, "corr P BW en txt.txt")

###Venn diagram
##ps.perc2
ps.perc2.BW <- subset_samples(ps.perc2, Plant == "BW")
ps.perc2.BS <- subset_samples(ps.perc2, Plant == "BS")
ps.perc2.RO <- subset_samples(ps.perc2, Plant == "RO")

ps.perc2.RO1 = prune_taxa(taxa_sums(ps.perc2.RO) > 0, ps.perc2.RO)
ps.perc2.BW1 = prune_taxa(taxa_sums(ps.perc2.BW) > 0, ps.perc2.BW)
ps.perc2.BS1 = prune_taxa(taxa_sums(ps.perc2.BS) > 0, ps.perc2.BS)

Bw1 <- rownames(otu_table(ps.perc2.BW1))
BS <- rownames(otu_table(ps.perc2.BS1))
RO <- rownames(otu_table(ps.perc2.RO1))

venn(list(RO, Bw1, BS))
Venn
#Create a venneuler object
age.venn=venneuler(c('A' = 566+485+1613+740, 'B' = 1805+626+1613+70, 'C' = 545+626+485+1613, 'A&B' = 740+1613, 'B&C' = 626+1613, 'A&C' = 485+1613, 'A&B&C' = 1613))

#Add group names
age.venn$labels = c("BW", "BS", "RO")

#Plot
plot(age.venn)

## CC networks based on MicrobiomeSeq
physeq <- taxa_level(ps.3, which_level = "Rank6")

co_occr <- co_occurence_network(physeq, grouping_column = "Plant", rhos =  "0.75", 
                                scale.vertex.size=3, scale.edge.width=15, 
                                plotBetweennessEeigenvalue=TRUE)

require(visNetwork)

g <- co_occr$net$graph
data <- toVisNetworkData(g)
visNetwork(nodes = data$nodes, edges = data$edges, width = 900)%>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)    

taxa.roles <- module.roles(co_occr$net$graph)

