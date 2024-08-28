# source: https://microsud.github.io/microbiomeutilities/articles/microbiomeutilities.html
# https://microbiome.github.io/tutorials/core_venn.html
# install.packages("devtools")
# devtools::install_github("microsud/microbiomeutilities")
library(microbiomeutilities)
# --------------------------------------------------------------------------
chip_LACM_99821 <- prune_samples(sample_data(chipmunk.only)$LACM == "99821", chipmunk.only)
chip_LACM_99839 <- prune_samples(sample_data(chipmunk.only)$LACM == "99839", chipmunk.only)
pins_only <- prune_samples(sample_data(chipmunk.only)$Parts == "pinworm", chipmunk.only)

table(meta(chip_LACM_99821)$Parts, useNA = "always")
table(meta(chip_LACM_99839)$Parts, useNA = "always")
table(meta(pins_only)$LACM, useNA = "always")

# convert to relative abundances ------------------------------------------
pseq.rel <- microbiome::transform(chip_LACM_99821, "compositional")
sampletype <- unique(as.character(meta(pseq.rel)$LACM.number))
print(sampletype)
library(eulerr)
list_core <- c() # an empty object to store information
for (n in sampletype){ # for each variable n in DiseaseState
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(pseq.rel, sampletype == n) # Choose sample from Parts by n
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.001, # 0.001 in at least 90% samples 
                         prevalence = 99/100)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each variable
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}

print(list_core)


# plot venn diagram --------------------------------------------------------
chip_LACM_99821_plot <- plot(venn(list_core), fills = c("#BEAED4","#fc8d59"))
chip_LACM_99839_plot <- plot(venn(list_core), fills = c("#BEAED4","#fc8d59"))
venpinplot <- plot(venn(list_core), fills = c("#fdc1a4","#fc8d59"))
vennALLplot <- plot(venn(list_core), fills = c("#fdc1a4","#fc8d59"))

# getting rows ------------------------------------------------------------
rows <- c(list_core$LACM.99839.feces, list_core$LACM.99839.pw) 
data_frame <- as.data.frame(tax_table(chipmunk.only))
setwd("/Users/psanjuan/Documents/Zymo_submission/output/figs")

# LACM_99821 --------------------------------------------------------------
# extracting data frame rows to get taxa names
data_mod <- data_frame[rownames(data_frame) %in% rows, ]
data_mod$ASVname <- row.names(data_mod)
print(data_mod)
write.csv(data_mod, "LACM_99821_core_0.001det_0.50prev.csv")

# extracting data frame rows to get taxa names for feces
LACM_99821_row <- list_core$LACM.99821.feces
data_mod_LACM_99821 <- data_frame[rownames(data_frame) %in% LACM_99821_row, ]  
data_mod_LACM_99821$ASVname <- row.names(data_mod_LACM_99821)
print(data_mod_LACM_99821)
write.csv(data_mod_LACM_99821, "LACM_99821_feces_core_0.001det_0.50prev.csv")

# extracting data frame rows to get taxa names for pinworms
LACM_99821_row_pin <- list_core$LACM.99821.pw
data_mod_LACM_99821_pin <- data_frame[rownames(data_frame) %in% LACM_99821_row_pin, ]  
data_mod_LACM_99821_pin$ASVname <- row.names(data_mod_LACM_99821_pin)
print(data_mod_LACM_99821_pin)
write.csv(data_mod_LACM_99821_pin, "LACM_99821_pin_core_0.001det_0.50prev.csv")

# LACM_99839 --------------------------------------------------------------
# extracting data frame rows to get taxa names
data_mod <- data_frame[rownames(data_frame) %in% rows, ]
data_mod$ASVname <- row.names(data_mod)
print(data_mod)
write.csv(data_mod, "LACM_99839_core_0.001det_0.50prev.csv")

# extracting data frame rows to get taxa names for feces
LACM_99839_row <- list_core$LACM.99839.feces
data_mod_LACM_99839 <- data_frame[rownames(data_frame) %in% LACM_99839_row, ]  
data_mod_LACM_99839$ASVname <- row.names(data_mod_LACM_99839)
print(data_mod_LACM_99839)
write.csv(data_mod_LACM_99839, "LACM_99839_feces_core_0.001det_0.50prev.csv")

# extracting data frame rows to get taxa names for pinworms
LACM_99839_row_pin <- list_core$LACM.99839.pw
data_mod_LACM_99839_pin <- data_frame[rownames(data_frame) %in% LACM_99839_row_pin, ]  
data_mod_LACM_99839_pin$ASVname <- row.names(data_mod_LACM_99839_pin)
print(data_mod_LACM_99839_pin)
write.csv(data_mod_LACM_99839_pin, "LACM_99839_pin_core_0.001det_0.50prev.csv")

# pinworms only -----------------------------------------------------------
# extracting data frame rows to get taxa names
data_mod <- data_frame[rownames(data_frame) %in% rows, ]
data_mod$ASVname <- row.names(data_mod)
print(data_mod)
write.csv(data_mod, "pinworm_core_0.001det_0.50prev.csv")

# extracting data frame rows to get taxa names for LAF.14610.pw
LACM_99821pw_row <- list_core$LACM.99821.pw
data_mod_LACM_99821pw <- data_frame[rownames(data_frame) %in% LACM_99821pw_row, ]  
data_mod_LACM_99821pw$ASVname <- row.names(data_mod_LACM_99821pw)
print(data_mod_LACM_99821pw)
write.csv(data_mod_LACM_99821pw, "pinworm_LACM_99821_core_0.001det_0.50prev.csv")

# extracting data frame rows to get taxa names for LAF.14628.pw
LACM_99839_row_pin <- list_core$LACM.99839.pw
data_mod_LACM_99839_pin <- data_frame[rownames(data_frame) %in% LACM_99839_row_pin, ]  
data_mod_LACM_99839_pin$ASVname <- row.names(data_mod_LACM_99839_pin)
print(data_mod_LACM_99839_pin)
write.csv(data_mod_LACM_99839_pin, "pinworm_LACM_99839_core_0.001det_0.50prev.csv")

