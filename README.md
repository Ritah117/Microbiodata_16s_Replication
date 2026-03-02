# Microbiodata_16s_Replication
This repository contains a robust R-based workflow for processing microbial sequences. It extends the standard DADA2 pipeline by integrating Cutadapt for precise primer removal and a Multi-step Taxonomic Rescue strategy using NCBI BLAST and taxonomizr.
🧬 Methodology Overview
The pipeline is designed to handle the specific challenges of bee gut microbiomes, such as high host DNA contamination (Mitochondria) and dietary DNA (Chloroplasts).

Pipeline Overview
Preprocessing: Quality control and primer trimming (Cutadapt).

DADA2 Core: Filtering, error learning, denoising, and merging.

Chimera Removal: Identification and removal of bimeras.

Taxonomy Assignment: Initial classification using the SILVA v138 database.

Taxonomic Rescue: Using NCBI BLAST and taxonomizr to identify ASVs unclassified at the Genus level.

Data Integration: Creating phyloseq objects for downstream ecological analysis.

Step 1: Preprocessing & Primer Trimming
Before DADA2 can accurately model sequencing errors, non-biological sequences (primers) must be removed. We use Cutadapt to perform "anchored" trimming, ensuring we only keep reads where primers are found in the correct orientation.

Logic: Forward and reverse primers are converted into all possible orientations (complement, reverse, and reverse-complement) to catch any flipped reads.

Execution:

Forward Primer: CCTACGGGNGGCWGCAG

Reverse Primer: GACTACHVGGGTATCTAATCC

Parameters: -n 2 (search for 2 primers), --discard-untrimmed (removes reads without primers).

Step 2: The DADA2 Core Workflow
This is the "heart" of the pipeline where raw reads are converted into ASVs (Amplicon Sequence Variants). Unlike traditional OTU picking, DADA2 resolves differences of a single nucleotide.

Sub-Step,Process,Metric for Success
Filter & Trim,Truncates reads at quality drops (truncLen).,High % of reads retained after filtering.
Error Learning,"Machine learning of error rates (A→C, etc.).",Convergence of the error model (visualized via plotErrors).
Denoising,Application of the core DADA2 algorithm.,"Identification of ""true"" biological variants."
Merging,Overlapping forward and reverse reads.,Overlap length matches expected V3-V4 insert size.

Chimera Removal
PCR can create "Frankenstein" sequences where two different biological templates fuse. We use the consensus method to identify and remove these artifacts.

Impact: In this dataset, we typically retain ~90-95% of reads as non-chimeric.

Track Table: We generate a track dataframe to monitor read loss at every single stage of the pipeline to ensure no step is too aggressive.

Initial Taxonomy (SILVA v138)
We first assign taxonomy using the SILVA v138 database. This provides a high-confidence "baseline" classification from Kingdom down to Species.

taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138_train_set.fa.gz", multithread=TRUE)

 Step 5: Taxonomic Rescue (NCBI BLAST)One of the most advanced features of this pipeline is the Taxonomic Rescue. Many bee-specific microbes are poorly represented in SILVA but present in the broader NCBI database.Isolation: Filter for all ASVs where Genus == NA.Export: Convert these specific sequences to a .fasta file.BLAST search: Query the NCBI GenBank database.Re-integration: Use the taxonomizr package to convert Accession numbers into a full taxonomic path (Phylum $\to$ Genus) and fill the "NAs" in the original table.


Step 6: Data Integration & Phyloseq
The final step assembles all the separate pieces of information into a single, unified phyloseq object.

Components:

OTU Table: Cleaned ASV counts.

Taxonomy Table: The "Rescued" SILVA + BLAST classifications.

Sample Data: Metadata (Species, Site, Treatment).

Final Cleaning: We perform a final "Taxonomic Scrub" to remove Chloroplasts, Mitochondria, and Archaea, ensuring the analysis focuses strictly on the bacterial microbiome.

Ensure you have R installed along with the following packages:
```bash
install.packages(c("tidyverse", "vegan", "ggpubr", "UpSetR", "devtools"))
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("dada2", "phyloseq", "DECIPHER", "phangorn", "ggtree", "Biostrings"))
```
# Initial Data Loading & Parsing
First, load the necessary library and verify the version to ensure reproducibility.
```bash
library(dada2)
# Check version: Ensure you are using DADA2 1.20+ for optimal performance
packageVersion("dada2")
```
1.Defining File Paths

Specify the local directory containing your paired-end FASTQ files.
```bash
# Define the path to your raw data
file_path <- "C:/Users/mukam/Desktop/test_data"

# List files to verify the path is correct
list.files(file_path)
```
2.Sorting Forward and Reverse Reads

We use pattern matching to separate forward (_R1) and reverse (_R2) files. Sorting is crucial to ensure that the $i^{th}$ file in dataF corresponds to the $i^{th}$ file in dataR.
```bash
# Extract forward read file paths
dataF <- sort(list.files(file_path, pattern="_R1_001.fastq.gz", full.names = TRUE))
dataF

# Extract reverse read file paths
dataR <- sort(list.files(file_path, pattern="_R2_001.fastq.gz", full.names = TRUE))
dataR
```
3. Extracting Sample Names
   
We parse the filenames to extract unique sample identifiers. The sample name is the second part of the filename string after splitting by an underscore (_).

```bash
# Extract sample names from filenames
# Adjust the number in `[` to match your file naming convention
list.sample.names <- sapply(strsplit(basename(dataF), "_"), `[`, 2)                
list.sample.names
```
4.Generating Quality Plots

We use DADA2 to generate quality profiles. This allows us to visualize the quality scores of the raw reads, which informs the decision for parameter settings in the filterAndTrim step (e.g., where to truncate reads).

```bash
# Plot quality profiles for the first 5 samples
# Green line: Mean quality score
# Orange line: Median quality score
# Dashed blue line: Percentage of reads reaching that position
plotQualityProfile(dataR[1:5])
```
# Primer Identification
This step prepares the pipeline to identify and remove primer sequences. Removing primers is critical because they contain degenerate bases ("N", "W", "R", etc.) that do not match the biological sequence, which can disrupt the DADA2 denoising algorithm.

We load essential bioinformatics libraries and define the target primer sequences for the V3-V4 region of the 16S rRNA gene.
```bash
library(ShortRead)
library(Biostrings)
library(stringr)
```
#Defining Primer Sequences

These sequences are defined based on the specific PCR primers used in the library preparation.

```bash
# Define Forward Primer
fwd_primer <- "CCTACGGGNGGCWGCAG"

# Define Reverse Primer
rev_primer <- "GACTACHVGGGTATCTAATCC"
```
#Defining Primer Orientations

We define a function that takes a primer sequence and generates all four possible orientations. This ensures that no matter how the sequence was read by the sequencer, the primer will be identified and removed.

```bash
# Function to generate all possible orientations of a primer
allOrients <- function(primer) {
  require(Biostrings)
  # Convert character string to Biostrings DNAString object
  dna <- DNAString(primer)
  
  # Generate orientations
  orients <- c(Forward = dna, 
               Complement = Biostrings::complement(dna), 
               Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  
  # Convert back to character vector for easier processing later
  return(sapply(orients, toString))
}
```
#Generating Primer Orientations

We utilize the allOrients function to create character vectors for both the forward and reverse primers in all four possible orientations: Forward, Complement, Reverse, and Reverse Complement.

```bash
# Generate orientations for the forward primer
fwd_primer_orients <- allOrients(fwd_primer)

# Generate orientations for the reverse primer
rev_primer_orients <- allOrients(rev_primer)

# Display the orientations for the reverse primer
rev_primer_orients
```
#Reverse Complement Specifics

While the function generates all orientations, we specifically generate and store the Reverse Complement (revComp) of both primers to be used later as specific arguments for the Cutadapt trimming tool to ensure precise anchored trimming.
```bash
# Reverse complement of the forward primer
fwd_primer_rev <- as.character(reverseComplement(DNAStringSet(fwd_primer)))
fwd_primer_rev

# Reverse complement of the reverse primer
rev_primer_rev <- as.character(reverseComplement(DNAStringSet(rev_primer)))
rev_primer_rev
````
#sanity check before removing primers

This step acts as a sanity check. Before attempting to remove primers, we must verify that they exist in the raw data and determine in which orientations they appear. This ensures our trimming strategy is correct.

We define a function to count how many reads contain specific primer orientations and apply this function to the first sample in our dataset (both forward and reverse reads).

This function reads a FASTQ file, searches for a primer sequence pattern, and counts the number of reads that contain at least one instance of that pattern.

```bash
# Function to count reads with specific primer pattern
count_primers <- function(primer, filename) {
  # Load reads and count patterns
  num_hits <- vcountPattern(primer, sread(readFastq(filename)), fixed = FALSE)
  # Return sum of reads with >0 hits
  return(sum(num_hits > 0))
}
```
We use sapply to run this function across all defined orientations for a representative sample file.

```bash
# Count primers in the FIRST FORWARD read file (dataF[[1]])
rbind(R1_fwd_primer = sapply(fwd_primer_orients, count_primers, filename = dataF[[1]]), 
      R2_fwd_primer = sapply(fwd_primer_orients, count_primers, filename = dataF[[1]]), 
      R1_rev_primer = sapply(rev_primer_orients, count_primers, filename = dataF[[1]]), 
      R2_rev_primer = sapply(rev_primer_orients, count_primers, filename = dataF[[1]]))

# Count primers in the FIRST REVERSE read file (dataR[[1]])
rbind(R1_fwd_primer = sapply(fwd_primer_orients, count_primers, filename = dataR[[1]]), 
      R2_fwd_primer = sapply(fwd_primer_orients, count_primers, filename = dataR[[1]]), 
      R1_rev_primer = sapply(rev_primer_orients, count_primers, filename = dataR[[1]]), 
      R2_rev_primer = sapply(rev_primer_orients, count_primers, filename = dataR[[1]]))
```
# Primer Removal & Quality Filtering
#Configuring Cutadapt

We specify the path to the cutadapt executable and create a directory to store the output files.
```bash
# Specify the path to cutadapt on your machine
cutadapt <- "C:/Users/mukam/miniconda3/envs/cutadapt_env/Scripts/cutadapt.exe"

# Test if cutadapt is found and working
system2(cutadapt, args = "--version")

# Define path for trimmed files
file_path <- "C:/Users/mukam/Desktop/test_data"
cut_dir <- file.path(file_path, "cutadapt")
if (!dir.exists(cut_dir)) dir.create(cut_dir)

# Create paths for output files
fwd_cut <- file.path(cut_dir, basename(dataF)) 
rev_cut <- file.path(cut_dir, basename(dataR))
names(fwd_cut) <- list.sample.names
names(rev_cut) <- list.sample.names

# Path for log files
cut_logs <- path.expand(file.path(cut_dir, paste0(list.sample.names, ".log")))
```
#Running Cutadapt

We define the arguments for trimming and loop through all samples to apply them.
```bash
# Define cutadapt arguments
# -g: forward primer, -a: reverse primer complement
# -G: reverse primer, -A: forward primer complement
cutadapt_args <- c("-g", fwd_primer, "-a", rev_primer_rev, 
                   "-G", rev_primer, "-A", fwd_primer_rev,
                   "-n", 2,"-m",1, "-j",32, "--discard-untrimmed") 

# Execute cutadapt
for (i in seq_along(dataF)) {
  system2(cutadapt, 
          args = c(cutadapt_args,
                   "-o", fwd_cut[i], "-p", rev_cut[i], 
                   dataF[i], dataR[i]),
          stdout = cut_logs[i])  
}
```
#Sanity Check & Pre-filtering Quality Plot

Verify that primers were removed and check the quality of the trimmed reads.
```bash
# Sanity check: count primers in trimmed files (should be 0)
rbind(R1_fwd_primer = sapply(fwd_primer_orients, count_primers, filename = fwd_cut[[1]]), 
      R2_fwd_primer = sapply(fwd_primer_orients, count_primers, filename = fwd_cut[[1]]), 
      R1_rev_primer = sapply(rev_primer_orients, count_primers, filename = fwd_cut[[1]]), 
      R2_rev_primer = sapply(rev_primer_orients, count_primers, filename = fwd_cut[[1]]))

# Visualize quality after trimming
plotQualityProfile(fwd_cut[1:5])
plotQualityProfile(rev_cut[1:5])
```
#DADA2 Quality Filtering (filterAndTrim)

We filter for quality and remove low-quality reads.
```bash
# Define paths for filtered files
filt.dataF <- file.path(file_path, "filtered", paste0(list.sample.names, "_F_filt.fastq.gz"))
filt.dataR <- file.path(file_path, "filtered", paste0(list.sample.names, "_R_filt.fastq.gz"))
names(filt.dataF) <- list.sample.names
names(filt.dataR) <- list.sample.names

# Run DADA2 Filtering
out <- filterAndTrim(fwd_cut, filt.dataF, rev_cut, filt.dataR, 
                     truncLen=c(240,200),
                     maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
head(out)
```
#Learning Error Rates

DADA2 requires an error model to distinguish between sequencing errors and true biological variants.
```bash
# Learn error rates
errF <- learnErrors(filt.dataF, multithread=TRUE)
errR <- learnErrors(filt.dataR, multithread=TRUE)

# Visualize error model
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
```
# Denoising

We apply the error model (errF, errR) to the filtered reads (filt.dataF, filt.dataR) to infer the true sequences.

```bash
# Apply DADA2 algorithm to infer ASVs
dadaFs <- dada(filt.dataF, err=errF, multithread=FALSE)
dadaRs <- dada(filt.dataR, err=errR, multithread=FALSE)
```
#Merging Paired-End Reads

We merge the denoised forward and reverse reads to reconstruct the full V3-V4 amplicon.

```bash
# Merge denoised forward and reverse reads
merge.reads <- mergePairs(dadaFs, filt.dataF, dadaRs, filt.dataR, verbose=TRUE)
```
#Sequence Table Construction

We construct a sequence table—a high-resolution analog of the traditional "OTU table"—which records the abundance of each ASV in each sample.

```bash
# Construct sequence table
seqtab <- makeSequenceTable(merge.reads)

# Inspect dimensions of the table (samples x ASVs)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
```
# Chimera Removal & Pipeline Tracking
Chimeras are sequences formed by the fusion of two different templates during PCR. Although DADA2 significantly reduces these, a final identification and removal step is necessary.

#Identifying and Removing Chimeras

We use the consensus method to identify chimeric sequences and remove them from our sequence table.

```bash
# Remove chimeric sequences
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# Calculate the percentage of reads that are non-chimeric
sum(seqtab.nochim)/sum(seqtab)
```
#Tracking Read Retention

It is crucial to monitor how many reads are lost at each step (filtering, denoising, merging, chimera removal).

```bash
# 1. Recover sample names
sample.names <- rownames(seqtab.nochim)

# 2. Define a helper function to count unique sequences
getN <- function(x) sum(getUniques(x))

# 3. Assemble the tracking table
track <- cbind(out, # From filterAndTrim
               sapply(dadaFs, getN), # Denoised Forward
               sapply(dadaRs, getN), # Denoised Reverse
               sapply(merge.reads, getN), # Merged reads
               rowSums(seqtab.nochim)) # Non-chimeric reads

# 4. Label columns and rows
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
```

#Downloading Taxonomic Databases

Before assigning taxonomy, we must download the formatted SILVA v138 databases. These files allow DADA2 to map our ASV sequences to known bacterial lineages.

```bash
# Download SILVA training set for genus-level assignment
download.file("https://zenodo.org/record/3986799/files/silva_nr99_v138_train_set.fa.gz", "silva_nr99_v138_train_set.fa.gz")

# Download SILVA species assignment file
download.file("https://zenodo.org/record/3986799/files/silva_species_assignment_v138.fa.gz", "silva_species_assignment_v138.fa.gz")

# Download SILVA training set with species information
download.file("https://zenodo.org/record/3986799/files/silva_nr99_v138_wSpecies_train_set.fa.gz","silva_nr99_v138_wSpecies_train_set.fa.gz")
```
#Taxonomic Assignment & Custom Curation

We assign taxonomic ranks (Kingdom to Genus) to our ASVs using the SILVA database. We also prepare a customized FASTA header list for future phylogenetic tree construction and load a finalized taxonomy table that likely includes specific classifications for bee-associated microbes.

#Assigning Taxonomy (SILVA)
```bash
# Assign taxonomy using the SILVA training set
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138_train_set.fa.gz", multithread=TRUE)

# Optional: Add species-level assignment if desired
# taxa <- addSpecies(taxa, "silva_species_assignment_v138.fa.gz")

# Create a copy for printing/viewing
taxa.print <- taxa
rownames(taxa.print) <- NULL # Remove rownames for cleaner display
head(taxa.print)
```
#Creating ASV Headers

We create generic, numbered headers (>ASV_1, >ASV_2, etc.) for our sequences. This is a crucial step for maintaining traceability when exporting sequences for BLAST or phylogenetic tree building.

```bash
# Initialize character vector for ASV headers
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

# Generate headers based on the number of ASVs
for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# Preview the first few headers
head(asv_headers)
```
#Importing Custom Taxonomy

```bash
# Load the finalized, custom-curated taxonomy table
tax_final <- read.csv("Final_StinglessBee_Taxonomy_Complete.csv", row.names = 1)
```
#Data Alignment & Phyloseq Object Assembly

Before we can analyze the data, we must align the sample identifiers across our three main data tables: Counts (OTU), Taxonomy (TAX), and Metadata (SAM).

#Aligning Tables

We check the number of samples in each table and filter them to include only the intersection—samples present in both tables.
```bash
# Check how many samples are in each (to see the difference)
ncol(count_asv_tab)          # Number of samples in your sequence table
nrow(sdata1)                 # Number of samples in your metadata

# Filter the Count Table to ONLY include samples present in your metadata
common_samples <- intersect(colnames(count_asv_tab), sdata1$sample_id)
count_asv_tab_filtered <- count_asv_tab[, common_samples]

# Filter the Taxonomy to match the ASVs present in the filtered counts
# (This prevents 'taxa mismatch' errors later)
taxa_filtered <- taxa.print[rownames(count_asv_tab_filtered), ]
```
# Check how many samples are in each (to see the difference)
ncol(count_asv_tab)          # Number of samples in your sequence table
nrow(sdata1)                 # Number of samples in your metadata

# Filter the Count Table to ONLY include samples present in your metadata
common_samples <- intersect(colnames(count_asv_tab), sdata1$sample_id)
count_asv_tab_filtered <- count_asv_tab[, common_samples]

# Filter the Taxonomy to match the ASVs present in the filtered counts
# (This prevents 'taxa mismatch' errors later)
taxa_filtered <- taxa.print[rownames(count_asv_tab_filtered), ]
```
#Building the Phyloseq Object

We convert our R objects into phyloseq compatible formats and assemble the final object.

```bash
# Set the sample names on the metadata safely
SAM <- sample_data(sdata1)
rownames(SAM) <- sdata1$sample_id

# Build the components
OTU <- otu_table(count_asv_tab_filtered, taxa_are_rows = TRUE)
TAX <- tax_table(taxa_filtered)

# ASSEMBLE the final object
physeq_final <- phyloseq(OTU, TAX, SAM)

# Check the result
physeq_final
```
#Finalizing & Exporting the Taxonomy Table

We take the curated BLAST results, format them as a matrix, align the row names with our ASV headers, and export them for documentation and downstream analysis.

```bash
# Load the curated taxonomic data (likely from BLAST/Taxonomizr step)
tax_final <- read.csv("Final_StinglessBee_Taxonomy_Complete.csv", row.names = 1)

# Convert to a matrix for compatibility with phyloseq
taxa.print <- as.matrix(tax_final)

# Align row names with ASV IDs (removing the ">" symbol from headers)
rownames(taxa.print) <- sub(">", "", asv_headers)

# Export the finalized taxonomy table as a CSV
write.csv(taxa.print, file="ASVs_taxonomy.csv")
```
#Phyloseq Object Assembly

We load the necessary libraries, format our tables to align perfectly, and assemble the phyloseq object.

#Preparing Components
```bash
library(phyloseq)
library(dplyr)
library(tibble)

# 1. Create the Taxonomy and OTU tables from our cleaned data
TAX = tax_table(taxa.print) # Taxonomy table created in Part 13
OTU = otu_table(count_asv_tab, taxa_are_rows = TRUE) # Feature table from Part One

# 2. Reading the sample metadata into R
sdata <- read.csv("stingless_bee_sample_metadata.csv", sep = ',', header = TRUE) 

# 3. Clean the metadata column names and set Row Names
# Note: Ensure 'sample_id' matches the column names in your OTU table
colnames(sdata) <- c("sample_id", "species") 

sdata1 <- sdata %>% 
  remove_rownames() %>% 
  column_to_rownames(var = "sample_id")

samdata = sample_data(sdata1)
```

#Creating the Phyloseq Object

```bash
# 4. Creating the final phyloseq object
physeq = phyloseq(OTU, TAX, samdata)

# 5. Verify the object
physeq
```
# Data Decontamination & Filtering

We utilize phyloseq to subset our taxa and samples based on specific taxonomic ranks and abundance thresholds.

#Taxonomic Filtering

We remove unwanted lineages by subsetting the tax_table. We use is.na() to ensure that sequences that are unclassified at a high taxonomic rank are not accidentally removed.

```bash
# 1. Filter out Chloroplasts (Plant DNA from pollen)
physeq1 <- subset_taxa(physeq, (Order != "Chloroplast") | is.na(Order))
cat("Taxa remaining after removing Chloroplasts:", ntaxa(physeq1), "\n")

# 2. Filter out Chloroflexi (Often considered contaminants in these samples)
physeq2 <- subset_taxa(physeq1, (Phylum != "Chloroflexi") | is.na(Phylum))
cat("Taxa remaining after removing Chloroflexi:", ntaxa(physeq2), "\n")

# 3. Filter out Mitochondria (Host/Bee DNA)
physeq3 <- subset_taxa(physeq2, (Family != "Mitochondria") | is.na(Family))
cat("Taxa remaining after removing Mitochondria:", ntaxa(physeq3), "\n")

# 4. Filter out Archaea (Focusing strictly on Bacteria)
physeq4 <- subset_taxa(physeq3, (Kingdom != "Archaea") | is.na(Kingdom))
cat("Taxa remaining after removing Archaea:", ntaxa(physeq4), "\n")
```

#Sample & Abundance Filtering

We remove the negative control sample and prune ASVs with extremely low read counts (potential noise).

```bash
# 5. Remove the Negative Control (NC68) 
clean_physeq = subset_samples(physeq4, sample_names(physeq4) != "NC68")

# Final summary of your clean dataset
clean_physeq

# Filtering ASVs that have fewer than 5 total reads across all samples
filtered_physeq <- prune_taxa(taxa_sums(clean_physeq) > 5, clean_physeq)

# Check how many taxa were removed
message("Original taxa count: ", ntaxa(clean_physeq))
message("Filtered taxa count: ", ntaxa(filtered_physeq))

# View the final object
filtered_physeq
```
#Data Preparation for Visualization

To create meaningful taxonomic barplots, we must normalize the data to relative abundance (proportions) so that samples with different sequencing depths can be directly compared.

#Data Transformation and Merging

We extract the tables from our filtered phyloseq object, perform the transformation, and merge them.

```bash
# 1. Get the Taxonomy table
tax_df <- as.data.frame(tax_table(filtered_physeq))

# 2. Get the Count/Abundance table and transform to Relative Abundance
# This converts raw counts to proportions (0 to 1)
ps_rel <- transform_sample_counts(filtered_physeq, function(x) x / sum(x))
otu_df <- as.data.frame(otu_table(ps_rel))

# 3. Merge them together by their Row Names (ASV IDs)
# This creates a comprehensive dataframe for plotting
final_table <- cbind(tax_df, otu_df)

# 4. View the result to verify structure
head(final_table)
```
#Taxonomic Aggregation and Averages

We utilize the tidyverse suite of packages to transform the long-format data into a summarized table showing the mean relative abundance of each genus across different bee species.

Data Manipulation Pipeline
```bash
library(tidyverse)

# 1. Aggregate ASV data to Genus level
# Group by Genus and sum abundances for each sample
genus_data <- final_table %>%
  group_by(Genus) %>%
  summarise(across(where(is.numeric), sum))

# 2. Pivot to Long format for joining with metadata
long_data <- genus_data %>%
  pivot_longer(cols = -Genus, names_to = "sample_id", values_to = "abundance")

# 3. Join with metadata to connect abundances to bee species
metadata <- as.data.frame(as.matrix(sample_data(filtered_physeq)))
metadata$sample_id <- rownames(metadata)

merged_species <- long_data %>%
  left_join(metadata, by = "sample_id")

# 4. Calculate Mean Abundance per Species and reshape
species_averages <- merged_species %>%
  group_by(Genus, species) %>%
  summarise(mean_abundance = mean(abundance, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = species, values_from = mean_abundance)

# Replace NAs (e.g., if a genus was absent in a species) with 0
species_averages[is.na(species_averages)] <- 0
```
#Visualizing Taxonomic Composition

We use ggplot2 to create a stacked barplot representing the relative abundance of the top 20 genera in each sample.

#Data Preparation for Plotting

We ensure the abundance data and taxonomic classifications are merged and then filter for the most abundant taxa.

```bash
library(tidyverse)

# 1. Convert row names to a column called ASV_ID and reshape to long format
bar <- otu_df %>%
  rownames_to_column(var = "ASV_ID") %>% 
  pivot_longer(cols = -ASV_ID, names_to = "SampleID", values_to = "abundance")

# 2. Ensure taxonomic data has ASV_ID column for merging
if(!"ASV_ID" %in% colnames(tax_df)){
  tax_df <- tax_df %>% rownames_to_column(var = "ASV_ID")
}

# 3. Join abundance data with taxonomy
bar <- bar %>% left_join(tax_df, by = "ASV_ID")

# 4. Find the Top 20 Genera based on total abundance
top_20_genera_names <- bar %>%
  group_by(Genus) %>%
  summarise(total = sum(abundance)) %>%
  arrange(desc(total)) %>%
  slice(1:20) %>%
  pull(Genus)

# 5. Filter data for top 20 and aggregate abundance by Genus
species_averages <- bar %>%
  filter(Genus %in% top_20_genera_names) %>%
  group_by(Genus, SampleID) %>%
  summarise(abundance = sum(abundance), .groups = 'drop')

# Quick check of prepared data
head(species_averages)
```
#Plotting with ggplot2

```bash
library(ggplot2)

# Create the stacked barplot
ggplot(species_averages, aes(x = SampleID, y = abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "fill") + # "fill" creates relative abundance (0-1)
  theme_minimal() +
  labs(title = "Top 20 Genera Relative Abundance",
       x = "Sample Name",
       y = "Relative Abundance (%)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  # Apply a color palette for distinct genera
  scale_fill_manual(values = grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(20))
```
#Refining the Taxonomic Composition Plot

We now update our data structure to specifically isolate the top 20 genera, consolidate the rest, and order the factors for logical plotting.

```bash
# 1. Identify the Top 20 Genera by total abundance
top_genera <- bar %>%
  group_by(Genus) %>%
  summarise(total_abundance = sum(abundance)) %>%
  arrange(desc(total_abundance)) %>%
  slice(1:20) %>%
  pull(Genus)

# 2. Re-calculate abundances, grouping rare taxa into "Other"
bar_filtered <- bar %>%
  mutate(Genus = as.character(Genus)) %>%
  mutate(Genus = ifelse(Genus %in% top_genera, Genus, "Other")) %>%
  group_by(Genus, SampleID) %>% 
  summarise(abundance = sum(abundance), .groups = 'drop')

# 3. Order the Genus factor to put "Other" at the bottom of the plot
ordered_genera <- c(setdiff(top_genera, "Other"), "Other")
bar_filtered$Genus <- factor(bar_filtered$Genus, levels = rev(ordered_genera))

# 4. Define a palette for 21 categories (20 genera + 1 "Other")
myPalette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(21)
```
#Generating the Refined Plot

```bash
# 5. Generate the Plot
ggplot(bar_filtered, aes(x = SampleID, y = abundance, fill = Genus)) +
  geom_col(position = "fill") + # "fill" ensures 0 to 1 relative abundance
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = myPalette) + 
  labs(title = "Taxonomic Composition by Bee Species (Top 20 Genera)",
       subtitle = "Rare taxa grouped into 'Other'",
       x = "Sample ID",
       y = "Relative Abundance")
```

# Preparing Data for Beta Diversity (Venn Diagram)

To analyze beta diversity, we first need to restructure our aggregated abundance data into a list of "present" genera for each sample. This allows us to compare the core microbiome versus group-specific microbes.

#Data Restructuring Pipeline

```bash
library(tidyverse)

# 1. Convert species_averages back to a list of 'Present' genera
# This filters out genera with 0 abundance for each sample
species_lists <- species_averages %>%
  filter(abundance > 0) %>%
  group_split(SampleID) # This creates a list dataframe for each sample

# 2. Name the list elements by the SampleID for reference
names(species_lists) <- sapply(species_lists, function(x) as.character(unique(x$SampleID)))

# 3. Extract just the Genus names for each group
# This creates the clean list of names needed for Venn diagram tools
venn_input <- lapply(species_lists, function(x) x$Genus)

# 4. View the result to ensure it is structured correctly
str(venn_input)
```

# UpSet Plot Visualization & Core Microbiome Analysis
An UpSet plot clearly shows the intersection of genera across all 10 bee samples.

Data Preparation for UpSet Plot

We convert the species_averages dataframe into a list format where each list item is a character vector of genera present in that sample.

```bash
library(tidyverse)

# 1. Identify sample columns (excluding the 'Genus' column)
sample_names <- setdiff(colnames(species_averages), "Genus")

# 2. Automatically create lists of present Genus for each sample
venn_list <- lapply(sample_names, function(sp) {
  # Get Genus names where abundance > 0
  as.character(species_averages$Genus[species_averages[[sp]] > 0 & !is.na(species_averages$Genus)])
})

# 3. Name the list elements by Sample Name
names(venn_list) <- sample_names

# 4. Check the counts of unique Genera per sample
print("Number of Genera per sample:")
print(sapply(venn_list, length))
```
#Generating the UpSet Plot
```bash
# 1. Load UpSetR library
if (!require("UpSetR")) install.packages("UpSetR")
library(UpSetR)

# 2. Convert the list to the binary matrix format required by UpSetR
upset_input <- fromList(venn_list)

# 3. Generate the plot
upset(upset_input, 
      nsets = length(venn_list), 
      order.by = "freq", 
      decreasing = TRUE,
      main.bar.color = "#2c3e50", # Color of the intersection bars
      sets.bar.color = "#e74c3c", # Color of the total set size bars
      point.size = 3.5, 
      line.size = 1,
      mb.ratio = c(0.6, 0.4), # Ratio of intersection bar chart to matrix
      text.scale = c(1.3, 1.3, 1.2, 1.2, 1.5, 1.2)) # Scale text size
```

#Identifying the Core Microbiome

We find the intersection of all lists to identify genera shared by every single sample

```bash
# 4. Find the Core Genera (Shared by every single sample)
core_genera <- Reduce(intersect, venn_list)
core_genera <- core_genera[!is.na(core_genera)]

# 5. Print the results clearly
cat("\n--- Core Microbiome Analysis ---\n")
cat("Total number of shared genera:", length(core_genera), "\n")
cat("The Core Genera shared by all species are:\n")
print(core_genera)
```
# Comprehensive Beta Diversity Analysis (Full Taxa Venn)

We now use the complete dataset to build a Venn diagram that shows the intersection of all microbial genera between the main host categories.

#Data Aggregation and Venn Preparation

```bash
# 1. Create a wide table for all 581 taxa, aggregating at the Genus level
full_species_wide <- bar %>%
  group_by(Genus, SampleID) %>%
  summarise(abundance = sum(abundance), .groups = 'drop') %>%
  pivot_wider(names_from = SampleID, values_from = abundance, values_fill = 0)

# 2. Create lists of present genera for each individual sample
full_venn_list <- lapply(setdiff(colnames(full_species_wide), "Genus"), function(samp) {
  as.character(full_species_wide$Genus[full_species_wide[[samp]] > 0])
})
names(full_venn_list) <- setdiff(colnames(full_species_wide), "Genus")

# 3. Group lists by Species category (DS, HA, HG)
# This merges all samples within a category together
vd_full <- list(
  DS = unique(unlist(full_venn_list[grep("DS", names(full_venn_list))])),
  HA = unique(unlist(full_venn_list[grep("HA", names(full_venn_list))])),
  HG = unique(unlist(full_venn_list[grep("HG", names(full_venn_list))]))
)
```
#Plotting and Interpretation

```bash
# 4. Check the total counts of genera per category
cat("Full Microbiome Genera counts:\n")
print(sapply(vd_full, length))

# 5. Plot the full Venn diagram
# Note: Ensure the 'gplots' or 'venn' package is installed and loaded
venn(vd_full)
```
# phylogenetic analysis
#Exporting Sequences to FASTA Format

To build a phylogenetic tree (e.g., for FastTree or IQ-TREE), we need the raw DNA sequences of the ASVs associated with their corresponding ASV IDs and taxonomic identification.

#Generating FASTA Headers and Sequences

We create a structured dataframe linking ASV IDs, raw sequences, and Genus-level taxonomy, then format them into FASTA lines.
```bash
library(phyloseq)

# 1. Get the 'Taxa Names' from your phyloseq object
# (These are typically the actual DNA sequences in DADA2)
raw_ids <- taxa_names(filtered_physeq)

# 2. Get the Taxonomy
asv_tax <- as.data.frame(as.matrix(tax_table(filtered_physeq)))

# 3. Create a data frame linking IDs, Sequences, and Taxonomy
seq_data <- data.frame(
  OTU = paste0("ASV_", 1:length(raw_ids)), # Create clean names like ASV_1, ASV_2
  Sequence = raw_ids,                      # Use the raw ID as the sequence
  Genus = asv_tax$Genus,
  stringsAsFactors = FALSE
)

# 4. Safety Check: Ensure IDs are actually DNA sequences (A,T,C,G)
if(!grepl("^[ATCGatcg]+$", seq_data$Sequence[1])) {
  cat("IDs are not sequences. Pulling from original 'seqs' object instead...\n")
  # seq_data$Sequence <- as.character(seqs[1:nrow(seq_data)]) # Uncomment if needed
}
```
#Writing the FASTA File
```bash
# 5. Build the FASTA lines
fasta_lines <- c()
for (i in 1:nrow(seq_data)) {
  # Clean Genus name (handle NAs)
  genus_name <- ifelse(is.na(seq_data$Genus[i]), "Unknown", seq_data$Genus[i])
  
  # Create descriptive header: >ASV_ID_Genus
  header <- paste0(">", seq_data$OTU[i], "_", genus_name)
  sequence <- seq_data$Sequence[i]
  fasta_lines <- c(fasta_lines, header, sequence)
}

# 6. Save and Verify
writeLines(fasta_lines, "bee_microbiome_phylo_final.fasta")

library(Biostrings)
phylo_sequences <- readDNAStringSet("bee_microbiome_phylo_final.fasta")

cat("\n--- Final Check: Total Sequences ---\n")
print(phylo_sequences)
```
# Phylogenetic Tree Construction
#Formatting and Saving Sequences

```bash
# 1. Initialize the FASTA lines
fasta_lines <- c()

# 2. Build the FASTA format
for (i in 1:nrow(seq_data)) {
  # Header: >ASV_ID_Genus (no spaces)
  header <- paste0(">", seq_data$OTU[i], "_", seq_data$Genus[i])
  sequence <- seq_data$Sequence[i]
  
  fasta_lines <- c(fasta_lines, header, sequence)
}

# 3. Save the file
writeLines(fasta_lines, "bee_microbiome_phylo.fasta")

# 4. Preview the first two sequences (4 lines)
print("Preview of the FASTA file:")
print(readLines("bee_microbiome_phylo.fasta", n = 4))
```
#Sequence Alignment (DECIPHER)

Raw sequences differ in length and contain insertions. Alignment matches homologous bases across all sequences.
```bash
library(Biostrings)
library(DECIPHER)

# Read the file
phylo_sequences <- readDNAStringSet("bee_microbiome_phylo.fasta")

# 1. Multiple sequence alignment
alignment <- AlignSeqs(phylo_sequences, anchor = NA)
```
# Tree Construction (phangorn)
We construct a Neighbor-Joining (NJ) tree as a starting point, then optimize it using a Maximum Likelihood (GTR+G+I) model.

```bash
library(phangorn)

# 2. Prepare data for Phangorn
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")

# 3. Calculate distances and build a starting tree (Neighbor-Joining)
dm <- dist.ml(phang.align)
treeNJ <- nj(dm) 

# 4. Fit the Maximum Likelihood model (GTR+G+I)
fit = pml(treeNJ, data=phang.align)

# 5. Optimize the tree (Rearrangement and substitution parameters)
fitGTR <- update(fit, k=4, inv=0.2) 
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0)) 

# 6. Save and extract the finalized tree
saveRDS(fitGTR, "stingless_bee_phangorn_tree.RDS")
phylo_tree <- fitGTR$tree

# ***Save workspace for future use
save(phang.align, treeNJ, fit, file = "tree_data.RData")
```
# Visualizing the Phylogenetic Tree
We will use a circular ("fan") layout to visualize the 581 ASVs, as this format is efficient for displaying large trees with many tips.

#Loading and Setting Up Visualization

```bash
# Set your working directory to where you saved the file
setwd("C:/Users/mukam/Desktop/")

# Load the optimized tree result
fitGTR <- readRDS("fitGTR_final.rds")

# Extract the tree object
phylo_tree <- fitGTR$tree

# Install and load required visualization packages
if (!require("ggtree")) BiocManager::install("ggtree")
library(ggtree)
library(ggplot2)
```
#Generating and Saving the Plot

```bash
# 1. Create the base circular tree
p <- ggtree(phylo_tree, layout="fan", size=0.3) +
  # Add labels for each ASV (keep size small to avoid overcrowding)
  geom_tiplab(size=0.5, aes(label=label), offset=0.01) + 
  theme_tree2() +
  labs(title = "Phylogenetic Tree of Stingless Bee Microbiome",
       subtitle = "Maximum Likelihood (GTR+G+I) - 581 ASVs")

# 2. Save as a high-resolution PDF for detailed viewing
ggsave("bee_microbiome_tree_highres.pdf", p, width=12, height=12)

# 3. Alternative base R plotting (circular layout)
pdf("Bee_Tree_Plot.pdf", width = 15, height = 15)
plot(phylo_tree, type = "fan", cex = 0.4)
dev.off()
```
# PCoA Analysis (Beta Diversity)
This section covers the beta diversity analysis using Principal Coordinate Analysis (PCoA). PCoA allows us to visualize the overall similarity or dissimilarity between the microbial communities of different samples based on the Bray-Curtis distance metric.

#Distance Calculation and Ordination

```bash
library(vegan)

# 1. Transpose the ASV table (Samples as Rows, ASVs as Columns)
features_pcoa <- t(count_asv_tab)

# 2. Verify sample count
cat("Number of samples for PCoA:", nrow(features_pcoa), "\n")

# 3. Calculate Bray-Curtis Distance
# Standard for microbiome abundance data
dist_matrix <- vegdist(features_pcoa, method = "bray")

# 4. Run the PCoA
pcoa_res <- cmdscale(dist_matrix, k = 2, eig = TRUE)

# 5. Calculate percentage of variation explained for plot axes
pcoa_variance <- round(pcoa_res$eig / sum(pcoa_res$eig) * 100, 1)
```
#Preparing Data for Plotting

We create a dataframe containing the coordinates for the first two principal components (PC1 and PC2) and map the sample IDs to bee species names.

```bash
# 6. Create the plotting data frame
pcoa_df <- data.frame(
  PC1 = pcoa_res$points[,1],
  PC2 = pcoa_res$points[,2],
  SampleID = rownames(features_pcoa)
)

# 7. Add Bee Species labels based on SampleID prefixes
pcoa_df$Species <- "Unknown"
pcoa_df$Species[grepl("DS", pcoa_df$SampleID)] <- "D. schmidti"
pcoa_df$Species[grepl("HA", pcoa_df$SampleID)] <- "H. araujo"
pcoa_df$Species[grepl("HG", pcoa_df$SampleID)] <- "H. gribodoi"

# Check the prepared data
head(pcoa_df)
```
# Phylogenetic Integration & Alpha Diversity Analysis

We will assemble a final phyloseq object that includes the phylogenetic tree, normalize the data, perform various ordinations (UniFrac and Bray-Curtis), and calculate diversity metrics on rarefied data to ensure comparability between samples.

#Finalizing the Phyloseq Object

We ensure that the ASV identifiers match across the OTU table, taxonomy table, and phylogenetic tree before assembly.

```bash
library(stringr)
library(phyloseq)
library(vegan)

# 1. Clean the Tree Labels to match ASV table row names
phylo_tree$tip.label <- str_extract(phylo_tree$tip.label, "ASV_\\d+")

# 2. Identify the matching ASVs across all components
common_asvs <- intersect(rownames(my_feature_table), rownames(my_taxonomy))
common_asvs <- intersect(common_asvs, phylo_tree$tip.label)
common_asvs <- common_asvs[!is.na(common_asvs)]

# 3. Filter and Build the finalized object
if(length(common_asvs) > 0) {
  OTU_final <- otu_table(my_feature_table[common_asvs, ], taxa_are_rows = TRUE)
  TAX_final <- tax_table(my_taxonomy[common_asvs, ])
  TREE_final <- keep.tip(phylo_tree, common_asvs)
  
  physeq3 <- phyloseq(OTU_final, TAX_final, samdata, TREE_final)
  
  cat("Success! Your phyloseq object is ready.\n")
  print(physeq3)
}
```
#Normalization & Ordination (Beta Diversity)

```bash
# Normalize to Relative Abundance
standardized_physeq = transform_sample_counts(physeq3, function(x) x / sum(x))

# Define colors for plotting
myPalette <- c("#E41A1C", "#377EB8", "#4DAF4A")

# --- Weighted UniFrac PCoA ---
ordu = ordinate(standardized_physeq, "PCoA", "unifrac", weighted=TRUE)
plot_ordination(standardized_physeq, ordu, color="Species", shape="Species") + 
  geom_point(size=5, alpha=0.7) +
  theme_bw() +
  scale_color_manual(values = myPalette) +
  labs(title = "Weighted UniFrac PCoA", subtitle = "Phylogenetic distances")

# --- Bray-Curtis PCoA ---
ordu_bray = ordinate(standardized_physeq, "PCoA", "bray")
plot_ordination(standardized_physeq, ordu_bray, color="Species", shape="Species") + 
  geom_point(size=4, alpha=0.8) +
  scale_color_manual(values = myPalette) +
  theme_bw() +
  labs(title = "Bray-Curtis PCoA", subtitle = "Community dissimilarity")
```
#Alpha Diversity Estimation (Rarefied)

To compare diversity fairly, we must rarefy samples to an even depth to remove bias from differing sequencing efforts.

```bash
# 1. Plot Rarefaction Curve to decide sampling depth
otu_tab <- t(as(otu_table(physeq4), "matrix"))
library(vegan)
rarecurve(otu_tab, step = 2000, label = TRUE, main = "Rarefaction Curves")

# 2. Rarefy to an even depth (example depth: 106927)
set.seed(9242)
rarefied_physeq <- rarefy_even_depth(physeq4, sample.size = 106927, replace = FALSE)

# 3. Calculate richness metrics
alpha_meas <- estimate_richness(rarefied_physeq, measures = c("Observed", "Shannon"))
```
#Statistical Testing & Visualization

We test for normality of the Shannon index and compare diversity between species using ANOVA.

```bash
library(tidyverse)
library(ggpubr)

# Prepare data for plotting
alpha_df <- alpha_meas %>% rownames_to_column("sample_id")
sdata1 <- as.data.frame(as(sample_data(physeq4), "matrix")) %>% rownames_to_column("sample_id")
shannon_editted <- merge(alpha_df %>% select(sample_id, Shannon), sdata1, by = "sample_id")

# Test for normality
shapiro.test(shannon_editted$Shannon)

# Generate Publication-Quality Plot
ggplot(shannon_editted, aes(x = species, y = Shannon)) + 
  geom_boxplot(aes(fill = species), alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5) + 
  theme_bw() +
  labs(title = "Alpha Diversity Comparison", y = "Shannon Diversity Index", x = "Bee Species") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"), legend.position = "none") +
  stat_compare_means(method = "anova", label.y = 2.5) # Compares means across groups
```

































































































































































































































