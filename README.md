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














