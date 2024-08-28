
# Load Libraries ----------------------------------------------------------
library(dada2)
library(phyloseq)
library(phangorn)
library(DECIPHER)
library(tidyr)
library(dplyr)
library(ggplot2)
library(vegan)
library(DESeq2)
library(microbiome)
#library(microbiomeutilities)
#library(ggpubr)

# Setting Path for input files --------------------------------------------
setwd('/Users/psanjuan/Documents/Zymo_submission/data/zr15679.rawdata.240217/raw_seqs/')
path='/Users/psanjuan/Documents/Zymo_submission/data/zr15679.rawdata.240217/raw_seqs/'
list.files(path)

# Specify forward & reverse read fastqs -----------------------------------
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))

# split filenames with '_', pick first character string
sample.names <- sapply(strsplit(basename(fnFs), "_R"), `[`, 1) 
sample.names

# Data quality metrics ----------------------------------------------------
# looking at forward fastq files, Phred scores
plotQualityProfile(fnFs[1:9]) # recommended to trim last 10 nucleotides (@ ~311)
plotQualityProfile(fnRs[1:9]) # reverse reads, trim @ ~271

# Filter and trim ---------------------------------------------------------
# Filtering reads based on sequence quality scores
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz")) # creating filtered folder
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

# Pulling sample names from filtered fasta files
names(filtFs) <- sample.names 
names(filtRs) <- sample.names

# set your primer sequence
# got primers from Zymo doc
# https://files.zymoresearch.com/protocols/_d6421_quick-16s_plus_ngs_library_prep_kit_(v3-v4_udi).pdf
FWD <-  "CCTACGGGDGGCWGCAG"  
#FWD2 <- "CCTAYGGGGYGCWGCAG"
REV <- "GACTACNVGGGTMTCTAATCC"
trimLeft = c(FWD,REV)

# Use known primer sequences to trim from your amplicon sequences
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(311,271),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE,trimLeft = c(17,24))
View(out)


# Learning Error rates ----------------------------------------------------
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# Plotting out the errors
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)


# Sample inference --------------------------------------------------------
# from forward reads
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)

# from reverse reads
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

# Inspect the dada objects produced above
dadaFs[[1]]
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)

# View mergers 
head(mergers[[1]])

# Making your ASV table ---------------------------------------------------
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Length of ASV's and their frequencies
table(nchar(getSequences(seqtab)))


# Identifying and removing chimeras ---------------------------------------
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus",multithread=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab) # tells us that chimeras account for ~10% of merged sequences


# Tracks reads through the pipeline ---------------------------------------
# Good checkpoint to ensure you did not lose too many reads
getN <- function(x) sum(getUniques(x))
track<- cbind(out, sapply(dadaFs,getN),sapply(dadaRs, getN), sapply(mergers,getN), rowSums(seqtab.nochim))
colnames(track) <- c("input","filtered","denoisedF","denoisedR","merged","nonchim")
rownames(track) <- sample.names
head(track)
write.csv(track, "Sequencing_Statistics_16S.csv")

# Assign taxonomy for 16S using native bayesian classifier built ----------
# Define your path correctly [path to be removed for practice]
taxa <- assignTaxonomy(seqtab.nochim, "/Users/psanjuan/Documents/Ref_databases/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)

# Attach taxonomic/species labels
taxa <- addSpecies(taxa, "/Users/psanjuan/Documents/Ref_databases/silva_species_assignment_v138.1.fa.gz")

# Check your taxonomy assignment
taxa.print <- taxa # Removing sequence rownames for display only
dim(taxa.print)
head(taxa.print)

# Giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

##  Giving taxonomy table corresponding names as above (ASV_1, ASV_2...)
row.names(taxa.print) <- sub(">", "", asv_headers)
write.table(taxa.print, "ASVs_named_correctly.tsv", sep="\t", quote=F, col.names=NA)
# Check the creation of output ASV files in your folder!!

# OPTIONAL: Check your taxonomic assignment. If you have unknowns or low level atx assignments that need to be updated by other means, you can do that here. and attach updated taxonomic labels as follows.

# END OF BIOINFORMATICS STEPS ---------------------------------------------
# NOW WE MOVE INTO PHYLOSEQ OBJECT CREATION & COMMUNITY STATISTICS

# Now we will create a sequence alignment to help improve clustering steps & phylogeny informed distances
sequences<- getSequences('ASVs.fa') # create sequence object

# make alignment, details here: https://rdrr.io/bioc/DECIPHER/man/AlignSeqs.html
alignment <- AlignSeqs(DNAStringSet(sequences), anchor=NA) 
phang.align <- phyDat(as(alignment, "matrix"), type="DNA") # Alignment matrix

# make maximum likelihood distance matrix
dm <- dist.ml(phang.align) 

# Create a neighbour joining tree as your starting tree
treeNJ <- NJ(dm) 

# Create model fit with pml, many options for model fitting are available, we will use GTR model as below.
fit <- pml(treeNJ, data = phang.align) 

# adjust model fit with the GTR model
fitGTR <- update(fit, k=4, inv=0.2) 


# Phyloseq object creation ------------------------------------------------
# Import Sample Metadata
map <- import_qiime_sample_data("/Users/psanjuan/Documents/Zymo_submission/sample_data_phyloseq.tsv") 

# Combine OTU table, Tax table, Sample metadata file and Phylogenetic tree to create a four dimensional phyloseq object
ps <- phyloseq(otu_table(asv_tab, taxa_are_rows = TRUE),
               tax_table(taxa.print),phy_tree(fitGTR$tree))
ps <- merge_phyloseq(ps, map)

# Check your phyloseq object, sample names, taxa names, sample sums,  taxa sums, rank names, number of samples, number of taxa, etc
# Self learning
ps
sample_names(ps)

# Check your phyla for fun & filter for uncharacterised reads
table(tax_table(ps)[, "Phylum"], exclude = NULL)
ps0 <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
ps0 #this new ps object should have 10 fewer taxa as 10 NA has been removed!

# Compute prevalence of each feature, store as data.frame 
prevdf = apply(X = otu_table(ps0), 
               MARGIN = ifelse(taxa_are_rows(ps0), yes = 1, no = 2), 
               FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame 
prevdf = data.frame(Prevalence = prevdf, TotalAbundance = taxa_sums(ps0), tax_table(ps0))

# Visualise prevalence of each phyla
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

# Filter for nuisance taxa, such as chloroplast or other known contaminants
# Taxa cleaning to remove sequences that aligned to chloroplast and mitochondria
ps1 <- ps0 %>%
  subset_taxa(
    Kingdom == "Bacteria" &
      Family != "Mitochondria" &
      Class != "Chloroplast" &
      Phylum != "Cyanobacteria" &
      Phylum != "Chloroplast" &
      Phylum != "Chloroflexi")

#Check ps object & if everything looks OK, rename to ps
ps1 

# Re-name original ps as new refined one
ps=ps1

# Create a sample sum table to look at coverage metrics
sample_sums(ps)
write.csv(sample_sums(ps),"16S_Reads_Per_Sample.csv") #outputs a coverage csv file

# Look for skew in coverage across sample types
set.seed(711)
level_order <- c('feces', 'pinworm', 'negative_control_DNA') #set your own variable order

DATA.2 <- ps  

df = as.data.frame(sample_data(DATA.2))
df$LibrarySize = sample_sums(DATA.2)
df = df[order(df$LibrarySize),]
df$Index = seq(nrow(df))


# Plot ordered library size across 'batch' & coloured by 'replicate'
ggplot(data=df, aes(x=Index, y=LibrarySize, colour= Species))+
  geom_point()+
  facet_wrap(~ factor(Parts, level = level_order)) +
  scale_y_continuous(trans='sqrt')


