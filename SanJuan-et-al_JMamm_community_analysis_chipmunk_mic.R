# Microbiome analysis -----------------------------------------------------
#   Author:       Priscilla San Juan
#   Date edited:  23 Feb 2024
#   Topic:        Microbiome pipeline for chipmunk microbiome 

# Table of contents -------------------------------------------------------
# 1  -  Load packages and color vectors 
# 2  -  Subset
# 3  -  Alpha diversity
# 4  -  Microbial Abundance 
# 5  -  Ordination
# 6  -  PermANOVA
# 7  -  Beta diversity
# 8  -  


# LOAD PACKAGES -----------------------------------------------------------
req_pkg <- c("readr","microbiome","dplyr","phyloseq","vegan","ggplot2","utils", 
             "DESeq2","ape","forcats","DescTools","scales","tidyr","ResourceSelection", 
             "gridExtra","ggbeeswarm","mvabund","RColorBrewer","randomcoloR","MASS", 
             "viridis","ggpubr","decontam","ggbeeswarm","fantaxtic","plotly") 

# Load all required packages and show version
for(i in req_pkg){
  print(i)
  print(packageVersion(i))
  library(i, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE, character.only=TRUE)
}

if(!"devtools" %in% installed.packages()){
  install.packages("devtools")
}
devtools::install_github("gmteunisse/fantaxtic")
devtools::install_github("gmteunisse/ggnested")

library(beepr)
require(microbiome)
library(ggrepel)
library(SRS)
library(ade4)        # aux
library(agricolae)   # aux
library(lmerTest)    # aux


# Set color palettes ------------------------------------------------------
n<-22
z<-24
library(RColorBrewer)
library(randomcoloR)
qual_col_pals=brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector=unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))
palette_qual <- distinctColorPalette(n)
pie(rep(1, n), col=palette_qual)
palette_qual_all <- distinctColorPalette(z)
pie(rep(1, z), col=palette_qual_all)
colors_scheme_hosts <- c("Tamias amoenus"="#BEAED4")
#   "Callospermophilus lateralis"="#fdb462",
#   "Tamias amoenus"="#BEAED4",
#   "Peromyscus maniculatus"="#7FC97F"
colors_chip_pin <- c("feces"="#BEAED4",
                     "pinworm"="#fc8d59")
colors_site <- c("Mill Creek"="#91bfdb",
                 "Whiskey Flats"="#d8b365")

# phyloseq
ps.noncontam

# Remove OTUs with zero counts
ps.noncontam.noZeros=prune_taxa(taxa_sums(ps.noncontam)>0, ps.noncontam) 

# Save an object to a file
# saveRDS(ps.noncontam.noZeros, file = "microbiome_phyloseq.rds")
# readRDS(file="microbiome_phyloseq.rds")

# Subset ------------------------------------------------------------------

chipmunk.only <- prune_samples(sample_data(ps.noncontam.noZeros)$Species == "Tamias amoenus", ps.noncontam.noZeros)
feces.only <- prune_samples(sample_data(ps.noncontam.noZeros)$Parts == "feces", ps.noncontam.noZeros)
pw.only <- prune_samples(sample_data(ps.noncontam.noZeros)$Parts == "pinworm", ps.noncontam.noZeros)
pw.only2 <- prune_taxa(taxa_sums(pw.only) > 0, pw.only)
chipfeces <- prune_samples(sample_data(feces.only)$Species == "Tamias amoenus", feces.only)


# Quick look at stats -----------------------------------------------------
# Most abundant OTUs within all samples
common.taxa.bac <- sort(taxa_sums(ps.noncontam.noZeros),T)
common.taxa.bac <- tax_table(ps.noncontam.noZeros)[names(common.taxa.bac[1:50])]
common.taxa.bac

# Statistics of read coverage across samples
max(sample_sums(chipmunk.only))
range(sample_sums(chipmunk.only))
sample_sums(chipmunk.only)
hist(sample_sums(chipmunk.only), breaks = 1000)


# Alpha diversity ---------------------------------------------------------
library(agricolae)
library(forcats)

# Host vs pinworms --------------------------------------------------------
# Reorder following the value of another column:
alpha.samples <- plot_richness(chipmunk.only, 
                            x="Parts", 
                            measures="Shannon",
                            #shape = "Parts",
                            color="Parts") + 
  geom_boxplot() + 
  #geom_point(size=3) + 
  scale_colour_manual(values = colors_chip_pin) + 
  geom_quasirandom(size = 3.0) + 
  #ggtitle("")+
  theme_bw() + xlab("Parts"); alpha.samples 

alpha.samples.edit <-   alpha.samples$data %>%
    ggplot(aes(x=Parts, y=value, col=Parts)) +
    geom_boxplot() + 
    geom_quasirandom(size = 3.0) + 
    scale_colour_manual(values = colors_chip_pin) +
    theme(legend.position = "none") +
    xlab("Parts") + ylab("Shannon Diversity Index") +
    theme_light() 

# Krustal-Wallis test (non-parametric)
kw_parts <- kruskal.test(value~Parts, data=alpha.samples$data)
print(kw_parts)

# load the dunn.test package for post-hoc test
library(dunn.test)
dunn.test(alpha.samples$data$value, alpha.samples$data$Parts, method = "bonferroni")

# Parts at phyla level ----------------------------------------------------
chipmunk.only.phyla <- tax_glom(chipmunk.only, taxrank = "Phylum")
chipmunk.only.phyla=prune_taxa(taxa_sums(chipmunk.only.phyla)>0, chipmunk.only.phyla) 

test <- plot_richness(chipmunk.only.phyla, 
              x="Parts", 
              measures="Shannon",
              #shape = "Parts",
              #measures=c("Observed","Shannon","Chao1"),
              color="Parts") + 
  geom_boxplot() + 
  #geom_point(size=3) + 
  scale_colour_manual(values = colors_chip_pin) + 
  geom_quasirandom(size = 3.0) + 
  #ggtitle("")+
  theme_bw() + xlab("Parts")
kruskal.test(value~Parts, data=test$data)
print(kw_parts)


# Site --------------------------------------------------------------------
alpha.site <- plot_richness(chipfeces, 
                           x="Collection_site", 
                           measures="Shannon", 
                           #measures=c("Observed","Shannon","Chao1"), 
                           color="Collection_site") + 
  # 
  geom_quasirandom(size = 3.0) + 
  scale_colour_manual(values=colors_site) +
  theme_light(); alpha.site
alpha.site$data %>%
  mutate(LAF_species = fct_relevel(LAF_number, "14610","14611","14628","14638")) %>%
  ggplot( aes(x=Collection_site, y=value, col=Collection_site, shape=Parts)) +
  # geom_boxplot() + 
  geom_quasirandom(size = 3.0) +
  scale_colour_manual(values = colors_site) +
  theme_light() +
  theme(legend.position = "none") +
  xlab("Sampling Locality") + ylab("Shannon Diversity Index") 

# Krustal-Wallis test (non-parametric)
kw_site <- kruskal.test(value~Collection_site, data=alpha.site$data)
print(kw_site)

# load the dunn.test package for post-hoc test
library(dunn.test)
dunn.test(alpha.site$data$value, alpha.site$data$Collection_site, method = "bonferroni")


# Site all samples ---------------------------------------------------------
all.site <- plot_richness(chipmunk.only, 
                            x="Collection_site", 
                            measures="Shannon", 
                            #measures=c("Observed","Shannon","Chao1"), 
                            color="Collection_site") + 
  # 
  geom_quasirandom(size = 3.0) + 
  scale_colour_manual(values=colors_site) +
  theme_light(); all.site
all.site$data %>%
  mutate(LAF_species = fct_relevel(LAF_number, "14610","14611","14628","14638")) %>%
  ggplot( aes(x=Collection_site, y=value, col=Collection_site, shape=Parts)) +
  # geom_boxplot() + 
  geom_quasirandom(size = 3.0) +
  scale_colour_manual(values = colors_site) +
  theme_light() +
  theme(legend.position = "none") +
  xlab("Sampling Locality") + ylab("Shannon Diversity Index") 

# Krustal-Wallis test (non-parametric)
kw_all.site <- kruskal.test(value~Collection_site, data=all.site$data)
print(kw_all.site)

# load the dunn.test package for post-hoc test
library(dunn.test)
dunn.test(alpha.site$data$value, alpha.site$data$Collection_site, method = "bonferroni")

# Bacterial Abundance -----------------------------------------------------
# Merge rare taxa to speed up examples
PS.rel <- microbiome::transform(chipmunk.only, "compositional")
PS.rel1 <- aggregate_rare(PS.rel, level = "Phylum", detection = 0/100, prevalence = 0/100)

PS.rel.melt <- PS.rel %>%
  psmelt
all_phylum.bac <- PS.rel1 %>%
  psmelt      
all_order.bac <- PS.rel2 %>%
  psmelt      

all_order.bac$Correct <- factor(all_order.bac$LACM.number, 
                                levels = c("LACM.99839.pw",
                                           "LACM.99839.feces",
                                           "LACM.99848.feces",
                                           "LACM.99821.pw",
                                           "LACM.99821.feces",
                                           "LACM.99822.feces"))

# Phyla
more_phyla_pal <- c("Actinobacteriota"="#6EE6D7","Bacteroidota"="#E0E9A3", 
                    "Campylobacterota"="#E1C5DF", "Deinococcota"="#DC8053",
                    "Firmicutes"="#C5EAD2", "Proteobacteria"="#5E6C92", 
                    "Deferribacterota"="#DB4DAE", "Patescibacteria"="#9BD371",
                    "Desulfobacterota"="#7B7FD9","Gemmatimonadota"="#7BACDA",
                    "Myxococcota"="#7A3EDF","Other"="#75E7D0",
                    "Planctomycetota"="#B9828C","Verrucomicrobiota"="#E3DA4A")

ggplot(all_phylum.bac, aes(x = Correct, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(~Collection_site, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = more_phyla_pal) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_light()+
  theme(axis.title.x = element_blank(), 
        strip.text = element_text(size = 10, face = "bold", color = "black",hjust = 0)) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance") +
  xlab("Sample name") + 
  ggtitle("Bacterial Phyla Composition") 

# Order
ggplot(all_order.bac, aes(x = Correct, y = Abundance, fill = Order)) + 
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(~Collection_site, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = palette_qual_all) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_light()+
  theme(axis.title.x = element_blank(), 
        strip.text = element_text(size = 10, face = "bold", color = "black",hjust = 0)) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance") +
  xlab("Sample name") +
  ggtitle("Bacterial Order Composition") 

# -------------------------------------------------------------------------
# BETA DIVERSITY: ORDINATION & PERMANOVA 
#Set seed for reproducibility
set.seed(2)

# Ordination --------------------------------------------------------------
chip.comp <- microbiome::transform(chipmunk.only, "compositional")
data.ordinaNMDS.bray <- ordinate(chip.comp, method="NMDS", distance="bray"); beep()

NMDS.bray <- plot_ordination((chip.comp),
                             data.ordinaNMDS.bray,
                             color="Collection_site",
                             shape="Parts", 
                             title="") + 
  theme_bw() + ggtitle("NDMS with Bray-Curtis Distance") +
  scale_colour_manual(values = colors_site) + 
  geom_label_repel(aes(label = LACM.number), 
                   label.padding = 0.1,
                   nudge_y = 0.05,
                   box.padding=0.15,
                   force=1) +
  geom_point(size=5, alpha=1) + 
  theme(axis.text = element_text(size=12), 
        axis.title = element_text(size=16)); NMDS.bray 


# PermANOVA ---------------------------------------------------------------
# chipmunk only samples - distance metric
df.chip <- as(sample_data(chipmunk.only),"data.frame")
dist.matrix.wha <- phyloseq::distance(chipmunk.only,"bray")
dist.matrix.uni <- phyloseq::distance(chipmunk.only,"unifrac")

# chipmunk poop only - distance metric
chippoo <- as(sample_data(chipfeces),"data.frame")
dist.matrix.wha2 <- phyloseq::distance(chipfeces,"bray")

# One-factor PERMANOVA
permANOVA.par <- adonis2(dist.matrix.wha~df.chip$Parts,
                         df.chip,
                         permutations = 100000); permANOVA.par
permANOVA.par2 <- adonis2(dist.matrix.uni~df.chip$Parts,
                         df.chip,
                         permutations = 100000); permANOVA.par
permANOVA.site <- adonis2(dist.matrix.wha~df.chip$Collection_site,
                          df.chip,
                          permutations = 100000); permANOVA.site
permANOVA.site2 <- adonis2(dist.matrix.wha2~chippoo$Collection_site,
                         chippoo,
                         permutations = 100000); permANOVA.site2


# DESeq2 ------------------------------------------------------------------
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
chip_DESeq <- phyloseq_to_deseq2(chipmunk.only, ~ Parts)
geoMeans = apply(counts(chip_DESeq), 1, gm_mean)
chip_DESeq = estimateSizeFactors(chip_DESeq, geoMeans = geoMeans)
chip_DESeq = DESeq(chip_DESeq, fitType="local")

res = results(chip_DESeq, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(chipmunk.only)[rownames(sigtab), ], "matrix"))
allsigtab = cbind(as(res, "data.frame"), as(tax_table(chipmunk.only)[rownames(res), ], "matrix"))
head(sigtab)
head(allsigtab)

# -------------------------------------------------------------------------
theme_set(theme_bw())
sigtabgen = subset(sigtab, !is.na(Genus))
sigtabfam = subset(sigtab, !is.na(Family))

# Phylum colors
phyla_pal <- c("Actinobacteriota"="#6EE6D7","Bacteroidota"="#E0E9A3", "Campylobacterota"="#7A3EDF", 
               "Deinococcota"="#DC8053","Firmicutes"="#C5EAD2", "Proteobacteria"="#5E6C92", 
               "Deferribacterota"="#DB4DAE", "Patescibacteria"="#9BD371")

taxa_overrep_cap <- ggplot(sigtabfam, aes(y=Family, x=log2FoldChange, color=Phylum, alpha=0.5, stroke=1)) + 
  geom_vline(xintercept = 0, linetype = "dotted", size = 1) +
  geom_point(size=6) + ylab("Bacterial Family") + xlab("log2FoldChange") +
  scale_color_manual(values = phyla_pal) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0, vjust=0.5), 
        axis.text = element_text(size=14), axis.title = element_text(size=15),
        legend.title.align = 0.5, legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 15)) + guides(color = guide_legend("Bacterial Phyla")) +
  ggtitle("Differential abundance of bacterial families")
taxa_overrep_cap


taxa_overrep_asvs <- ggplot(sigtabfam, aes(y=row.names(sigtabfam), x=log2FoldChange, color=Phylum, alpha=0.5, stroke=1)) + 
  geom_vline(xintercept = 0, linetype = "dotted", size = 1) +
  geom_point(size=6) + ylab("Bacterial Family") + xlab("log2FoldChange") +
  scale_color_manual(values = phyla_pal) +
        axis.text = element_text(size=14), axis.title = element_text(size=15),
        legend.title.align = 0.5, legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 15)) + guides(color = guide_legend("Bacterial Phyla")) +
  ggtitle("Differential abundance of bacterial families")
taxa_overrep_asvs

taxa_overrep_asvs$data$log2FoldChange

taxa_overrep_chip <- taxa_overrep_asvs$data[taxa_overrep_asvs$data$log2FoldChange < 0,]
taxa_overrep_pinw <- taxa_overrep_asvs$data[taxa_overrep_asvs$data$log2FoldChange > 0,]

# -------------------------------------------------------------------------
# Volcano plot with ASVs
allsigtab
ggplot(allsigtab, aes(y=-log10(pvalue), x=log2FoldChange, color=Phylum, stroke=1)) + 
  geom_label_repel(aes(label = rownames(allsigtab))) +
  ylab(label=("-log(10)pvalue")) +
  geom_vline(xintercept = 0, linetype = "dotted", size = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", size = 0.5) +
  xlab("log2FoldChange") +
  geom_point(size=4, alpha=0.4) + 
  scale_color_manual(values = phyla_pal)


# End of working ------------------------------------------------------------



