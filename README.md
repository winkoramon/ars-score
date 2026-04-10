# ARS Pipeline: A Framework for Functional lncRNA Variant Prioritization

Developed by **Korawich Uthayopas**, Ancestry and Health Genomics Laboratory, University of Sydney.

The **Accessibility, Regulation, and Structure (ARS)** framework is designed to identify and prioritize functional somatic variants within long non-coding RNAs (lncRNAs), specifically optimized for geo-ancestrally diverse prostate cancer cohorts.

## Framework Components
The pipeline calculates the **Driver Variant Impact Score (DVIS)** by integrating:
- **Structure (S):** Structural disruption probabilities via RNAsnp (global/local p-values).
- **Accessibility (A):** Changes in RNA base-pairing probabilities via RNAplfold.
- **Regulation (R):** Overlap with regulatory elements (promoters, enhancers) and strong TF motif disruptions (e.g., AR, FOXA1, HOXB13).

## Key Features
- Automated representative transcript selection (MANE-Select priority).
- Genomic-to-cDNA coordinate mapping for non-coding transcripts.
- Quantitative scoring (ARS) for variant prioritization.

## Citation
If you use this pipeline, please cite:
*Uthayopas, K., et al. (2026). 

# How to run
#1. Download all required data files (https://drive.google.com/drive/folders/1RJ4trF80l1EB1QPKqEd6tzfh85SNfbxg?usp=share_link) and put the folder 'Data' into root folder.

#2. Run the script
python 01_generate_scoring.py
python 02_combine_ars_score.py