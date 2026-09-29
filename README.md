# ARS Pipeline: A Framework for Functional lncRNA Variant Prioritization

Developed by **Korawich Uthayopas**, Ancestry and Health Genomics Laboratory, University of Sydney.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)


The **Accessibility, Regulation, and Structure (ARS)** framework is designed to identify and prioritize functional somatic variants within long non-coding RNAs (lncRNAs), specifically optimized for geo-ancestrally diverse prostate cancer cohorts.

## Contents

- [Overview](#overview)
- [Repository structure](#repository-structure)
- [1. System requirements](#1-system-requirements)
- [2. Installation guide](#2-installation-guide)
- [3. Demo](#3-demo)
- [4. Instructions for use](#4-instructions-for-use)
- [5. Reproducing the manuscript results](#5-reproducing-the-manuscript-results)
- [License](#license)
- [Citation](#citation)
- [Contact](#contact)

---

## Overview

For each somatic single-nucleotide variant (SNV) in a driver lncRNA, the pipeline integrates three kinds of evidence into a single integer **ARS score**:

| Component | Evidence | Rule | Points |
|---|---|---|---|
| **Structure (S)** | RNAsnp global p-value | p < 0.05 | +2 |
| | RNAsnp local p-value | p < 0.05 | +2 |
| **Accessibility (A)** | RNAplfold ΔPu (single nucleotide, `Pu1_delta`) | \|Δ\| > 0.1 | +2 |
| | RNAplfold ΔPu (8-nt window coverage, `Pu8cover_delta`) | \|Δ\| > 0.1 | +2 |
| **Regulation (R)** | Overlap with promoter | yes | +1 |
| | Overlap with proximal enhancer-like signature (pELS) | yes | +1 |
| | Overlap with distal enhancer-like signature (dELS) | yes | +1 |
| | Overlap with CTCF-bound element | yes | +1 |
| | Strong TF-motif disruption, per TF family (AR, FOXA, HOXB13, GATA, NKX3-1, ETS, TWIST1, ZEB1, SNAI1, SOX4), counted only for variants in a promoter or enhancer | yes | +1 each |

The maximum possible score is 22. Higher scores indicate stronger combined evidence of a functional effect.

Key features:

- Automated representative transcript selection (MANE Select priority).
- Genomic-to-cDNA coordinate mapping for non-coding transcripts.
- Quantitative ARS scoring for variant prioritization, with a publication-ready score-distribution plot.

## Repository structure

```
ars-score/
├── Code/
│   ├── 01_generate_scoring.py      # SNV → transcript mapping, sequence windows, regulatory & motif inputs
│   └── 02_combine_ars_score.py     # Combines S, A and R evidence into the ARS score and plots it
├── demo/
│   ├── make_demo_data.py           # Regenerates the simulated demo dataset
│   ├── Data/                       # Small SIMULATED dataset (no patient data)
│   └── expected_output_ARS_scored_driver_variants.csv
├── requirements.txt
├── LICENSE
└── README.md
```

---

## 1. System requirements

### Operating systems

The software is pure Python and runs on Linux and macOS.

| Tested on | Version |
|---|---|
| Linux | Ubuntu 24.04 LTS |
| macOS | *(add the macOS version you used)* |

Windows is not officially supported because `pysam` does not provide native Windows builds; Windows users can run the pipeline via WSL2.

### Software dependencies

| Dependency | Tested version |
|---|---|
| Python | 3.11.15 |
| numpy | 2.4.6 |
| pandas | 3.0.6 |
| pysam | 0.24.1 |
| scipy | 1.17.1 |
| seaborn | 0.13.2 |
| matplotlib | 3.11.2 |
| tqdm | 4.70.1 |

The exact versions are pinned in [`requirements.txt`](requirements.txt).

The structural and accessibility inputs used in the manuscript were produced with the following external tools. They are **not** needed to run the demo or to re-score the provided data.

| Tool | Used for | Version |
|---|---|---|
| RNAsnp | SNV effects on RNA secondary structure | *(add version)* |
| ViennaRNA (RNAplfold) | Base-pairing / accessibility probabilities | *(add version)* |

### Hardware

No non-standard hardware is required. The demo runs on any standard desktop or laptop. For the full dataset, which includes the GRCh38 reference genome (~3 GB), at least 8 GB of RAM and 10 GB of free disk space are recommended.

---

## 2. Installation guide

### Instructions

```bash
# 1. Clone the repository
git clone https://github.com/winkoramon/ars-score.git
cd ars-score

# 2. Create and activate a virtual environment (recommended)
python3 -m venv .venv
source .venv/bin/activate          # Windows (WSL): same command

# 3. Install the dependencies
pip install -r requirements.txt
```

Conda users can do the same with:

```bash
conda create -n ars python=3.11 -y
conda activate ars
pip install -r requirements.txt
```

### Typical install time

Under **2 minutes** on a normal desktop computer with a broadband connection. Installing the dependencies into a clean virtual environment took about 20 seconds in our tests.

---

## 3. Demo

A small, fully **simulated** dataset is provided in [`demo/Data`](demo/Data): 60 SNVs in three synthetic lncRNAs on a 20 kb synthetic contig, with simulated RNAsnp, RNAplfold, regulatory and TF-motif results. It contains no patient data. It can be regenerated with `python demo/make_demo_data.py` (fixed random seed).

### Instructions to run on the demo data

From the repository root:

```bash
python Code/01_generate_scoring.py --data-dir demo/Data
python Code/02_combine_ars_score.py --data-dir demo/Data --out-dir demo/output
```

### Expected output

The console should end with:

```
ARS Scoring Complete. Output saved to: demo/output/ARS_scored_driver_variants.csv
Generating publication-ready ARS distribution plot...
Plot saved to: demo/output/ars_distribution_1000dpi.png
```

Two files are written to `demo/output/`:

| File | Description |
|---|---|
| `ARS_scored_driver_variants.csv` | One row per SNV with the input columns plus `lncRNA_ancestry_type`, `snv_ancestry_type` and `ARS_Score` |
| `ars_distribution_1000dpi.png` | Histogram of ARS scores |

For the demo, the 60 SNVs have ARS scores from 0 to 8 (median 3) with the following distribution:

| ARS score | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|---|
| Number of SNVs | 6 | 1 | 16 | 12 | 14 | 7 | 2 | 1 | 1 |

You can check your result against the reference file:

```bash
python -c "import pandas as pd; a=pd.read_csv('demo/output/ARS_scored_driver_variants.csv'); b=pd.read_csv('demo/expected_output_ARS_scored_driver_variants.csv'); print('Demo output matches:', a.equals(b))"
```

### Expected run time

Under **10 seconds** on a normal desktop computer (about 2 seconds in our tests).

---

## 4. Instructions for use

### Running on your own data

Put your inputs into a data folder that follows the layout below, then run:

```bash
python Code/01_generate_scoring.py  --data-dir /path/to/Data
python Code/02_combine_ars_score.py --data-dir /path/to/Data [--out-dir /path/to/output]
```

If `--data-dir` is omitted, `./Data` is used. The data folder can also be set with the `ARS_DATA_DIR` environment variable. If `--out-dir` is omitted, results are written to `<data-dir>/VariantImpactScore/`.

### Input data layout

```
Data/
├── pca+_447_samples_snv_v6.csv                         # somatic SNVs
├── driver_lncRNA_dict.pkl                              # driver lncRNA lists
├── Driver_lncRNA/multivariate_snv_esults.pkl           # ancestry-associated SNVs
├── Genomic_ref/
│   ├── GRCh38.primary_assembly.genome.fa               # reference genome
│   ├── gencode.v47.long_noncoding_RNAs.gtf             # GENCODE lncRNA annotation
│   └── MANE_lncRNAs__MANE-Select__with_ENST_ENSG_IDs.csv
├── RNAstructure/input/
│   ├── RNAsnp_combined_results.csv
│   └── RNAplfold_accessibility_changes.csv
├── regulatory_summary/Table_S1_SNV_annotations.csv     # optional location, see below
└── TF_motif/tf_motif_disruption/per_snv_tf_calls_STRONG.csv
```

Required columns:

| File | Required columns |
|---|---|
| SNV table | `Gene ID`, `Chromosome`, `POS`, `REF`, `ALT` |
| `driver_lncRNA_dict.pkl` | dict with `compd_dlnc_list` (list of Ensembl gene IDs), `compd_dlnc_af_sp` and `compd_dlnc_eu_sp` (DataFrames with a `Gene` column) |
| `multivariate_snv_esults.pkl` | dict with `mtv_pos_af_snv` (list of `snv_<row index>` names) |
| MANE table | `ENSG_core`, `ENST_core` |
| RNAsnp results | `gene_id`, `chrom`, `gpos`, `ref_tx`, `alt_tx`, `pvalue_global`, `pvalue_local`, `user_ref_match`, `status_warning` |
| RNAplfold results | `gene_id`, `chrom`, `gpos`, `ref_tx`, `alt_tx`, `Pu1_delta`, `Pu8cover_delta` |
| Regulatory table | `gene_id`, `chrom`, `pos`, `ref`, `alt`, `in_promoter`, `hit_pELS`, `hit_dELS`, `hit_CTCF`, `hit_ELS` |
| TF-motif calls | `tf_family`, `gene_id`, `chrom`, `pos`, `ref`, `alt`, `strong_change` |

Gene IDs are matched without version suffixes (e.g. `ENSG00000123456`). The regulatory table is read from `<data-dir>/regulatory_summary/` if present; otherwise from `Code/regulatory_summary/`. A different path can be given with `--regmap`.

---

## 5. Reproducing the manuscript results

1. Download the full input data from [Google Drive](https://drive.google.com/drive/folders/1RJ4trF80l1EB1QPKqEd6tzfh85SNfbxg?usp=share_link).
2. Place the downloaded `Data` folder in the repository root, so that the layout matches [Input data layout](#input-data-layout).
3. Run:

   ```bash
   python Code/01_generate_scoring.py
   python Code/02_combine_ars_score.py
   ```

4. The ARS scores for all driver-lncRNA SNVs are written to `Data/VariantImpactScore/ARS_scored_driver_variants.csv`, and the ARS distribution figure is written to `Data/VariantImpactScore/ars_distribution_1000dpi.png`.

---

## License

This project is released under the [MIT License](LICENSE), an [Open Source Initiative](https://opensource.org/license/mit)-approved license.

## Citation

If you use this pipeline, please cite:

> Uthayopas, K., *et al.* (2026). *[Manuscript title].* *[Journal]*. DOI: *[to be added]*

## Contact

For questions or issues, please open an [issue](https://github.com/winkoramon/ars-score/issues) or contact Korawich Uthayopas (Ancestry and Health Genomics Laboratory, University of Sydney).

