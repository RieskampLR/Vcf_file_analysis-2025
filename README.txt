Readme

Lea Rachel Rieskamp
Project: VCF file based analysis of population differentiation and selection
Stepwise description and code (BASH)

Tools and Versions used:

VCFtools (0.1.16)
bcftools (1.3.1)
Shapeit (v2.r904)
SweepFinder2
R (4.2.1 (2022-06-23 ucrt))

# Data Processing:

# Data exploration
# Individuals
cat ProjTaxa.vcf | sed -n '2343p' | cut -f 10- |  awk '{print NF}'
# --> 16 + 1 ref
# Categories
# 8N, K, Lesina, Naxos2
# --> 4 (based on shared names and naming conventions)
# Chromosomes
cat ProjTaxa.vcf | sed -n '2344,$p' | cut -f 1 | sort | uniq
# --> chr5 & chrZ --> 2
# Data is not phased (/)
# Variants before filtering:
bcftools view -H ProjTaxa.vcf | wc -l
# --> 3816977

# Data Preparation:
# Download the VCF file (ProjTaxa.vcf.gz) from the server

# Setting filters:
# Create subsample using vcfrandomsample from vcflib and bcftools:
bcftools view ProjTaxa.vcf | vcfrandomsample -r 0.012 > ProjTaxa_subset.vcf
# Compress and index the subset VCF:
bgzip ProjTaxa_subset.vcf
bcftools index ProjTaxa_subset.vcf.gz
# Statistics of the vcf subset using vcftools
# Calculate allele frequency for each variant
vcftools --gzvcf ProjTaxa_subset.vcf.gz --freq2 --out ./vcftools_analysis/ProjTaxa_subset  --max-alleles 2
# --freq2 outputs the frequencies without info about the alleles (what we want), --freq would return their identity.
# Calculate mean depth of coverage per individual
vcftools --gzvcf ProjTaxa_subset.vcf.gz --depth --out ./vcftools_analysis/ProjTaxa_subset
# Calculate mean depth of coverage per site
vcftools --gzvcf ProjTaxa_subset.vcf.gz --site-mean-depth --out ./vcftools_analysis/ProjTaxa_subset
# Calculate site quality score for each site
vcftools --gzvcf ProjTaxa_subset.vcf.gz --site-quality --out ./vcftools_analysis/ProjTaxa_subset
# Calculate proportion of missing data per individual
vcftools --gzvcf ProjTaxa_subset.vcf.gz --missing-indv --out ./vcftools_analysis/ProjTaxa_subset
# Calculate proportion of missing data per site
vcftools --gzvcf ProjTaxa_subset.vcf.gz --missing-site --out ./vcftools_analysis/ProjTaxa_subset
# Calculate heterozygosity and inbreeding coefficient per individual (F)
vcftools --gzvcf ProjTaxa_subset.vcf.gz --het --out ./vcftools_analysis/ProjTaxa_subset

# R:
# Codes TaxaProject BINP28 - Filtering the vcf file - R
# Libraries
library(tidyverse)
### Variant based stats ###
# Load site quality vcftools output
var_qual <- read_delim("C:/Users/lea/OneDrive/Dokumente/UniLund/BINP28/Project/vcftools_analysis/ProjTaxa_subset.lqual", delim = "\t",
                       col_names = c("chr", "pos", "qual"), skip = 1)
# Plot distribution of site quality scores
a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light() + xlim(0, 5000)
# Checking around the standard minimum threshold of 30:
a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light() + xlim(0, 50)

# Load mean depth per site vcftools output
var_depth <- read_delim("C:/Users/lea/OneDrive/Dokumente/UniLund/BINP28/Project/vcftools_analysis/ProjTaxa_subset.ldepth.mean", delim = "\t",
                        col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
# Checking around main bulk of variants:
a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light() + xlim(0, 50)
# Check basic stats:
summary(var_depth$mean_depth)
# Load missing data per site vcftools output
var_miss <- read_delim("C:/Users/lea/OneDrive/Dokumente/UniLund/BINP28/Project/vcftools_analysis/ProjTaxa_subset.lmiss", delim = "\t",
                       col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(var_miss, aes(fmiss)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
# Check basic stats:
summary(var_miss$fmiss)
# Load allele frequency vcftools output
var_freq <- read_delim("C:/Users/lea/OneDrive/Dokumente/UniLund/BINP28/Project/vcftools_analysis/ProjTaxa_subset.frq", delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
# Minor allele frequencies
var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))
a <- ggplot(var_freq, aes(maf)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
# Check basic stats:
summary(var_freq$maf)

### Individual based stats ###
# Load mean depth per individual vcftools output
ind_depth <- read_delim("C:/Users/lea/OneDrive/Dokumente/UniLund/BINP28/Project/vcftools_analysis/ProjTaxa_subset.idepth", delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1)
a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()

# Load proportion of missing data per individual vcftools output
ind_miss <- read_delim("C:/Users/lea/OneDrive/Dokumente/UniLund/BINP28/Project/vcftools_analysis/ProjTaxa_subset.imiss", delim = "\t",
                       col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()

# Load heterozygosity and inbreeding coefficient per individual vcftools output
ind_het <- read_delim("C:/Users/lea/OneDrive/Dokumente/UniLund/BINP28/Project/vcftools_analysis/ProjTaxa_subset.het", delim = "\t",
                      col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)
a <- ggplot(ind_het, aes(f)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()



# Applying filters to whole vcf file:
# Define variables for filters
QUAL=30
MIN_DEPTH=7
MAX_DEPTH=18
MISS=0.9
MAF=0.1
# Filtering using vcftools (in vcftool_analysis directory)
vcftools --vcf ../ProjTaxa.vcf --remove-indels --min-alleles 2 --max-alleles 2 --maf $MAF --max-missing $MISS --minQ $QUAL --min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH --minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | gzip -c > ProjTaxa_prefiltered.vcf.gz
# Rename out.log file
mv out.log out_prefiltering.log

# Variants left after filtering:
bcftools view -H ProjTaxa_prefiltered.vcf.gz | wc -l
# --> 166556

# Find individuals with >10% missing data after the filtering:
gunzip ProjTaxa_prefiltered.vcf.gz
vcftools --vcf ProjTaxa_prefiltered.vcf --missing-indv --out indivs_with_high_missing_data
cat indivs_with_high_missing_data.imiss

# Remove outgroup (Naxos2) and Individuals with high missing data (K006):
vcftools --vcf ProjTaxa_prefiltered.vcf --recode --out ProjTaxa_filtered --remove-indv Naxos2 --remove-indv K006
# --> outputfile: ProjTaxa_filtered.recode.vcf (same number of variants)

# FST analysis:
# Prep population file (pop.txt):
# Populations are defined based on naming conventions: 8N, K, Lesina, & ref:
# 8N05240 pop1
# 8N05890 pop1
# 8N06612 pop1
# 8N73248 pop1
# 8N73604 pop1
# K006 pop2
# K010 pop2
# K011 pop2
# K015 pop2
# K019 pop2
# Lesina_280 pop3
# Lesina_281 pop3
# Lesina_282 pop3
# Lesina_285 pop3
# Lesina_286 pop3
ref_genome ref
# --> pop1: 8N, pop2: K, pop3: Lesina, ref: ref_genome

# Prep file for each population:
awk '$2 == "pop1" {print $1}' pop.txt > pop1.txt
awk '$2 == "pop2" {print $1}' pop.txt > pop2.txt
awk '$2 == "pop3" {print $1}' pop.txt > pop3.txt


# FST calculations:
# Make directory for outputs
mkdir Fst_analysis
# Run analyses for each population pair (pop1/2, pop1/3, pop2/3)
vcftools --vcf ProjTaxa_filtered.recode.vcf --weir-fst-pop pop1.txt --weir-fst-pop pop2.txt --fst-window-size 20000 --out ./Fst_analysis/Fst_pop1_pop2
vcftools --vcf ProjTaxa_filtered.recode.vcf --weir-fst-pop pop1.txt --weir-fst-pop pop3.txt --fst-window-size 20000 --out ./Fst_analysis/Fst_pop1_pop3
vcftools --vcf ProjTaxa_filtered.recode.vcf --weir-fst-pop pop2.txt --weir-fst-pop pop3.txt --fst-window-size 20000 --out ./Fst_analysis/Fst_pop2_pop3
# Analysis is done with sliding windows of 20kb size

R:
# Visualizing FSTs
# pop1: 8N
# pop2: K
# pop3: Lesina

# Libraries
library(tidyverse)

# Load FST data
fst_data_pop1_pop2 <- read.csv("C:/Users/lea/OneDrive/Dokumente/UniLund/BINP28/Project/Fst_analysis/Fst_pop1_pop2.windowed.weir.fst", header = TRUE, sep = "\t")
fst_data_pop1_pop3 <- read.csv("C:/Users/lea/OneDrive/Dokumente/UniLund/BINP28/Project/Fst_analysis/Fst_pop1_pop3.windowed.weir.fst", header = TRUE, sep = "\t")
fst_data_pop2_pop3 <- read.csv("C:/Users/lea/OneDrive/Dokumente/UniLund/BINP28/Project/Fst_analysis/Fst_pop2_pop3.windowed.weir.fst", header = TRUE, sep = "\t")


# Plots

# Pop1/2 - chr5:
ggplot(fst_data_pop1_pop2 %>% filter(CHROM == "chr5"), aes(x = BIN_START, y = WEIGHTED_FST)) +
  geom_line(size=0.3) +
  labs(x = "Position on chr5", y = "FST (8N/K)") +
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  theme(axis.line = element_line(color = "black")) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 6e+07, by = 1e+07)) +
  scale_y_continuous(limits = c(-0.2, 1), breaks = seq(0, 1, by = 0.2))

# Pop1/2 - chrZ:
ggplot(fst_data_pop1_pop2 %>% filter(CHROM == "chrZ"), aes(x = BIN_START, y = WEIGHTED_FST)) +
  geom_line(size=0.3) +
  labs(x = "Position on chrZ", y = "FST (8N/K)") +
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  theme(axis.line = element_line(color = "black")) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 6e+07, by = 1e+07)) +
  scale_y_continuous(limits = c(-0.2, 1), breaks = seq(0, 1, by = 0.2))

# Pop1/3 - chr5:
ggplot(fst_data_pop1_pop3 %>% filter(CHROM == "chr5"), aes(x = BIN_START, y = WEIGHTED_FST)) +
  geom_line(size=0.3) +
  labs(x = "Position on chr5", y = "FST (8N/Lesina)") +
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  theme(axis.line = element_line(color = "black")) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 6e+07, by = 1e+07)) +
  scale_y_continuous(limits = c(-0.2, 1), breaks = seq(0, 1, by = 0.2))

# Pop1/3 - chrZ:
ggplot(fst_data_pop1_pop3 %>% filter(CHROM == "chrZ"), aes(x = BIN_START, y = WEIGHTED_FST)) +
  geom_line(size=0.3) +
  labs(x = "Position on chrZ", y = "FST (8N/Lesina)") +
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  theme(axis.line = element_line(color = "black")) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 6e+07, by = 1e+07)) +
  scale_y_continuous(limits = c(-0.2, 1), breaks = seq(0, 1, by = 0.2))

# Pop2/3 - chr5:
ggplot(fst_data_pop2_pop3 %>% filter(CHROM == "chr5"), aes(x = BIN_START, y = WEIGHTED_FST)) +
  geom_line(size=0.3) +
  labs(x = "Position on chr5", y = "FST (K/Lesina)") +
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  theme(axis.line = element_line(color = "black")) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 6e+07, by = 1e+07)) +
  scale_y_continuous(limits = c(-0.2, 1), breaks = seq(0, 1, by = 0.2))

# Pop2/3 - chrZ:
ggplot(fst_data_pop2_pop3 %>% filter(CHROM == "chrZ"), aes(x = BIN_START, y = WEIGHTED_FST)) +
  geom_line(size=0.3) +
  labs(x = "Position on chrZ", y = "FST (K/Lesina)") +
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  theme(axis.line = element_line(color = "black")) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 6e+07, by = 1e+07)) +
  scale_y_continuous(limits = c(-0.2, 1), breaks = seq(0, 1, by = 0.2))


# SweepFinder2:
# Make directory for outputs:
mkdir Sweep_analysis
# Phase data with shapeit2:
# Separate vcf file for the 2 chromosomes
vcftools --vcf ProjTaxa_filtered.recode.vcf --chr chr5 --recode --out Sweep_analysis/ProjTaxa_filtered_chr5
vcftools --vcf ProjTaxa_filtered.recode.vcf --chr chrZ --recode --out Sweep_analysis/ProjTaxa_filtered_chrZ
# Phasing
# In Sweep_analysis directory
shapeit --input-vcf ProjTaxa_filtered_chr5.recode.vcf -O chr5_phased.vcf --window 0.5 -T 4
shapeit --input-vcf ProjTaxa_filtered_chrZ.recode.vcf -O chrZ_phased.vcf --window 0.5 -T 4
# chr5 number of haplotypes identified
cat chr5_phased.vcf.haps | wc -l
# --> 102775
# chrZ number of haplotypes identified
cat chrZ_phased.vcf.haps | wc -l
# --> 63781

# Convert outputs to vcf:
shapeit -convert --input-haps chr5_phased.vcf --output-vcf chr5_phased.vcf.vcf
shapeit -convert --input-haps chrZ_phased.vcf --output-vcf chrZ_phased.vcf.vcf

# Compress and index the vcfs:
bgzip chr5_phased.vcf.vcf
bgzip chrZ_phased.vcf.vcf
bcftools index chr5_phased.vcf.vcf.gz
bcftools index chrZ_phased.vcf.vcf.gz
# View
bcftools view -H chr5_phased.vcf.vcf.gz | head | cut -f 1-12
bcftools view -H chrZ_phased.vcf.vcf.gz | head | cut -f 1-12

# Reformat
bcftools +fill-tags chr5_phased.vcf.vcf.gz -Oz -o chr5_phased_with_tags.vcf.gz
bcftools +fill-tags chrZ_phased.vcf.vcf.gz -Oz -o chrZ_phased_with_tags.vcf.gz
bcftools query -f '%POS\t%ALT\t%AC\t%AN\n' chr5_phased_with_tags.vcf.gz > chr5_allele_counts.txt
bcftools query -f '%POS\t%ALT\t%AC\t%AN\n' chrZ_phased_with_tags.vcf.gz > chrZ_allele_counts.txt
awk '{print $1 "\t" $3 "\t" $4 "\t1"}' chr5_allele_counts.txt > reformatted_chr5_allele_counts.txt
awk '{print $1 "\t" $3 "\t" $4 "\t1"}' chrZ_allele_counts.txt > reformatted_chrZ_allele_counts.txt
# Add required headers
sed -i '1i position\tx\tn\tfolded' reformatted_chr5_allele_counts.txt
sed -i '1i position\tx\tn\tfolded' reformatted_chrZ_allele_counts.txt

# SweepFinder2
SweepFinder2 -s 1000 reformatted_chr5_allele_counts.txt chr5_sweepfinder_out.txt
SweepFinder2 -s 1000 reformatted_chrZ_allele_counts.txt chrZ_sweepfinder_out.txt

R:
# Visualizing Sweeps

# Libraries
library(ggplot2)

# Load SweepFinder2 data
sweep_data_chr5 <- read.table("C:/Users/lea/OneDrive/Dokumente/UniLund/BINP28/Project/Sweep_analysis/chr5_sweepfinder_out.txt", header = TRUE)
sweep_data_chrZ <- read.table("C:/Users/lea/OneDrive/Dokumente/UniLund/BINP28/Project/Sweep_analysis/chrZ_sweepfinder_out.txt", header = TRUE)


# Plots

# Chr5
ggplot(sweep_data_chr5, aes(x = location, y = LR)) +
  geom_line() +
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  theme(axis.line = element_line(color = "black")) +
  scale_x_continuous(expand = c(0, 0)) +
  ylim(-5, 20) +
  labs(x = "Position", y = "Likelihood Ratio (LR)", title = "Selective Sweep Analysis")

# ChrZ
ggplot(sweep_data_chrZ, aes(x = location, y = LR)) +
  geom_line() +
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  theme(axis.line = element_line(color = "black")) +
  scale_x_continuous(expand = c(0, 0)) +
  ylim(-5, 20) +
  labs(x = "Position", y = "Likelihood Ratio (LR)", title = "Selective Sweep Analysis")


# Tajimas D:
# Unzip phased vcf files
gunzip -c chr5_phased.vcf.vcf.gz > chr5_phased.vcf.vcf
gunzip -c chrZ_phased.vcf.vcf.gz > chrZ_phased.vcf.vcf
# Tajimas D with vcftools
vcftools --vcf chr5_phased.vcf.vcf --TajimaD 20000 --out chr5_Tajima
vcftools --vcf chrZ_phased.vcf.vcf --TajimaD 20000 --out chrZ_Tajima


R:
# Visualizing Tajimas D

# Libraries
library(ggplot2)

# Load data
tajimas_d_chr5 <- read.table("C:/Users/lea/OneDrive/Dokumente/UniLund/BINP28/Project/TajimasD_analysis/chr5_Tajima.Tajima.D", header = TRUE)
tajimas_d_chrZ <- read.table("C:/Users/lea/OneDrive/Dokumente/UniLund/BINP28/Project/TajimasD_analysis/chrZ_Tajima.Tajima.D", header = TRUE)

# Plots

# Chr5
ggplot(tajimas_d_chr5, aes(x = BIN_START, y = TajimaD)) +
  geom_line() +
  theme_minimal() +
  labs(x = "Genomic Position", y = "Tajima's D", title = "Tajima's D Across Genome")

# ChrZ
ggplot(tajimas_d_chrZ, aes(x = BIN_START, y = TajimaD)) +
  geom_line() +
  theme_minimal() +
  labs(x = "Genomic Position", y = "Tajima's D", title = "Tajima's D Across Genome")


# Combined plots (LD and Tajima´s D)

# Load SweepFinder2 data
sweep_data_chr5 <- read.table("C:/Users/lea/OneDrive/Dokumente/UniLund/BINP28/Project/Sweep_analysis/chr5_sweepfinder_out.txt", header = TRUE)
sweep_data_chrZ <- read.table("C:/Users/lea/OneDrive/Dokumente/UniLund/BINP28/Project/Sweep_analysis/chrZ_sweepfinder_out.txt", header = TRUE)


# Plot chr5
ggplot(tajimas_d_chr5, aes(x = BIN_START, y = TajimaD)) +
  geom_point(size=0.5, col="lightblue") +
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  theme(axis.line = element_line(color = "black")) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 6e+07, by = 1e+07)) +
  labs(x = "Genomic Position", y = "Tajima's D and LR/5", title = "Tajima's D and LR Across Chr5") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_line(data=sweep_data_chr5, aes(x = location, y = LR/5), col="darkblue", size=0.45) +
  ylim(-0.7, 4)


# Plot chrZ
ggplot(tajimas_d_chrZ, aes(x = BIN_START, y = TajimaD)) +
  geom_point(size=0.5, col="lightblue") +
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  theme(axis.line = element_line(color = "black")) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 6e+07, by = 1e+07)) +
  labs(x = "Genomic Position", y = "Tajima's D and LR/5", title = "Tajima's D and LR Across ChrZ") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_line(data=sweep_data_chrZ, aes(x = location, y = LR/5), col="darkblue", size=0.45) +
  ylim(-0.7, 4)















This script provides a complete workflow for analyzing the population structure and selection in House Sparrows. It processes the VCF data, applies filtering, performs FST and selection analyses, and generates plots for interpretation.
