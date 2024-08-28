# DECONTAM: FILTER CONTAMS USING NEG CONTROL ------------------------------
library(decontam)
library(plotly)

# Inspect library sizes ---------------------------------------------------
df <- as.data.frame(sample_data(ps))
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
tem <- ggplot(data=df, aes(x=Index, y=LibrarySize, color=Parts)) + geom_point()
ggplotly(tem)

# Use prevalence method to filter out suspected contaminants --------------
sample_data(ps)$is.neg <- sample_data(ps)$Parts == "negative_control_DNA"
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg", threshold = 0.5)
hist(contamdf.prev$p) # looked at histogram to determine cut-off, range from 0.1-0.7 is the same
table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant))

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Parts == "negative_control_DNA", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Parts == c("feces", "pinworm"), ps.pa)

# Make data.frame of prevalence in positive and negative samples ----------
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), 
                    pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)

ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

# Remove contams from phyloseq obj ----------------------------------------
ps.noncontam <- prune_taxa(!contamdf.prev$contaminant, ps)
ps.noncontam

# Remove the negative control from analysis -------------------------------
ps.noncontam <- prune_samples(sample_data(ps.noncontam)$Parts != "negative_control_DNA", ps.noncontam)

# ready to use in phyloseq ------------------------------------------------
ps.noncontam
