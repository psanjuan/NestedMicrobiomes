# Pinworm microbiomes are distinct from their chipmunk host gut microbiota
## Priscilla A. San Juan<sup>1(0000-0001-5341-6622)</sup>, Lizbeth Palma<sup>2</sup>, and Kayce Bell<sup>1(0000-0001-6835-0021)</sup>
1: Department of Mammalogy, Natural History Museum of Los Angeles, Los Angeles, California, USA
2: California State University, Dominguez Hills, Los Angeles, California, USA


## R code used to analyze chipmunk and pinworm gut microbiomes
### Files included are:
- *SanJuan-et-al_JMamm_dada2_bioinformatics.R*
  - Bioinformatics code, get from raw fastq files to a phyloseq object
- *SanJuan-et-al_JMamm_decontam_chipmunk_mic.R*
  - Processing step to remove any potential contaminants from the DNA extraction negative control
- *SanJuan-et-al_JMamm_community_analysis_chipmunk_mic.R*
  - Community ecology/microbiome analysis of chipmunk and pinworm gut microbiomes
- *SanJuan-et-al_JMamm_venn_chip_v_pw.R*
  - Comparing the core microbiome of chipmunks and pinworms to analyze shared and distinct taxa
