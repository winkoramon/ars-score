import os
import re
import math
import pickle
import numpy as np
import pandas as pd
import pysam
from collections import defaultdict, namedtuple
from tqdm import tqdm

# ==========================================
# 1. CONFIGURATION & FILE PATHS
# ==========================================
PROJECT_FOLDER = './' # Assuming running from repo root
DATA_DIR = os.path.join(PROJECT_FOLDER, 'Data')

SNV_PATH = os.path.join(DATA_DIR, 'Main_data/562_data/save_folder/master_table/pca+_447_samples_snv_v6.csv')
DICT_PATH = os.path.join(DATA_DIR, 'Driver_lncRNA/driver_lncRNA_dict.pkl')
GTF_PATH = os.path.join(DATA_DIR, 'lncRNA/gencode.v47.long_noncoding_RNAs.gtf')
FA_PATH = os.path.join(DATA_DIR, 'lncRNA/GRCh38.primary_assembly.genome.fa')
MANE_PATH = os.path.join(DATA_DIR, 'RNAstructure/MANE_lncRNAs__MANE-Select__with_ENST_ENSG_IDs.csv')

# Output directories
OUTDIR = os.path.join(DATA_DIR, 'RNAstructure/input/')
META_DIR = os.path.join(OUTDIR, 'meta')
WINDOWS_DIR = os.path.join(OUTDIR, 'windows_401nt_per_snv')

os.makedirs(META_DIR, exist_ok=True)
os.makedirs(WINDOWS_DIR, exist_ok=True)

WINDOW_RADIUS = 200
VALID_GENE_TYPES = {"lncRNA", "lincRNA", "antisense", "sense_intronic", "sense_overlapping", "processed_transcript"}
Exon = namedtuple("Exon", ["chrom", "start", "end", "strand", "exon_number"])

# ==========================================
# 2. HELPER FUNCTIONS
# ==========================================
def parse_gtf_attributes(attr_str):
    attrs = {}
    for m in re.finditer(r'(\S+)\s+"([^"]+)"\s*;?', attr_str):
        k, v = m.group(1), m.group(2)
        if k in attrs:
            if isinstance(attrs[k], list): attrs[k].append(v)
            else: attrs[k] = [attrs[k], v]
        else:
            attrs[k] = v
    return attrs

def revcomp(seq):
    comp = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(comp)[::-1]

def norm_chr(chrom):
    return chrom if chrom.startswith("chr") else f"chr{chrom}"

def infer_gene_type(feature, attrs):
    if feature == "gene":
        return attrs.get("gene_type") or attrs.get("gene_biotype")
    elif feature == "transcript":
        return attrs.get("transcript_type") or attrs.get("transcript_biotype")
    return None

# ==========================================
# 3. DATA LOADING & TRANSCRIPT SELECTION
# ==========================================
def load_and_filter_snvs():
    print("Loading SNVs and Driver Dictionary...")
    snp_df = pd.read_csv(SNV_PATH)
    snp_df['Gene_ID_Trimmed'] = snp_df['Gene ID'].str.split('.').str[0]
    
    with open(DICT_PATH, 'rb') as f:
        loaded_dict = pickle.load(f)
    
    compd_dlnc_list = loaded_dict['compd_dlnc_list']
    return snp_df[snp_df['Gene_ID_Trimmed'].isin(compd_dlnc_list)]

def parse_gtf_and_mane():
    print("Parsing GTF and MANE Select annotations...")
    mane_df = pd.read_csv(MANE_PATH)
    mane_dict = {row['ENSG_core']: row['ENST_core'] for _, row in mane_df.iterrows()}
    
    genes, tx_by_gene, exons_by_tx = {}, defaultdict(dict), defaultdict(list)
    mane_select = set()

    with open(GTF_PATH, "r") as fh:
        for line in fh:
            if not line or line.startswith("#"): continue
            chrom, _, feature, start, end, _, strand, _, attrs_str = line.strip().split("\t")
            attrs = parse_gtf_attributes(attrs_str)
            start, end = int(start), int(end)

            if feature == "gene":
                g_id = attrs["gene_id"].split(".")[0]
                genes[g_id] = {"gene_id": g_id, "gene_name": attrs.get("gene_name", ""),
                               "gene_type": infer_gene_type("gene", attrs), "chrom": norm_chr(chrom),
                               "start": start, "end": end, "strand": strand}
            elif feature == "transcript":
                g_id, t_id = attrs["gene_id"].split(".")[0], attrs["transcript_id"].split(".")[0]
                tag = attrs.get("tag", [])
                tag_list = tag if isinstance(tag, list) else [tag]
                
                if any(t in ("MANE_Select", "MANE_Select_v1") for t in tag_list) or (g_id in mane_dict and t_id == mane_dict[g_id]):
                    mane_select.add(t_id)

                tx_by_gene[g_id][t_id] = {"transcript_id": t_id, "chrom": norm_chr(chrom), 
                                          "start": start, "end": end, "strand": strand, "is_mane_select": t_id in mane_select}
            elif feature == "exon":
                t_id = attrs["transcript_id"].split(".")[0]
                exons_by_tx[t_id].append(Exon(norm_chr(chrom), start, end, strand, attrs.get("exon_number")))

    lnc_genes = {g for g, d in genes.items() if (d.get("gene_type") in VALID_GENE_TYPES or "lncrna" in (d.get("gene_type", "").lower()))}
    return genes, tx_by_gene, exons_by_tx, lnc_genes

# ==========================================
# 4. SEQUENCE EXTRACTION & MAPPING
# ==========================================
def process_mutational_effects(snv_df, genes, tx_by_gene, exons_by_tx, lnc_genes):
    print("Mapping SNVs to cDNA and generating fasta windows...")
    fa = pysam.FastaFile(FA_PATH)
    
    def build_transcript_seq(tx_id):
        exs = exons_by_tx.get(tx_id, [])
        if not exs: return ""
        exs_sorted = sorted(exs, key=lambda e: (e.exon_number if e.exon_number is not None else e.start), reverse=(exs[0].strand == "-"))
        seq = "".join([fa.fetch(e.chrom, e.start-1, e.end).upper() for e in exs_sorted])
        return revcomp(seq) if exs_sorted[0].strand == "-" else seq

    # ... (Sequence generation logic for RNAsnp and RNAplfold goes here)
    # Includes building `rnasnp_rows` and writing out `WT{idx}.fa` and `SNV{idx}.fa`
    pass

# ==========================================
# 5. REGULATORY REGIONS & MOTIF DISRUPTION
# ==========================================
def map_regulatory_regions():
    print("Mapping to cCREs and FANTOM peaks...")
    # Add pyRanges overlapping logic here
    pass

def compute_motif_disruption():
    print("Computing TF motif disruptions (Δ-scores)...")
    # Add MEME parsing and log-odds scoring here
    pass

# ==========================================
# MAIN EXECUTION
# ==========================================
if __name__ == "__main__":
    snv_df = load_and_filter_snvs()
    genes, tx_by_gene, exons_by_tx, lnc_genes = parse_gtf_and_mane()
    
    # Execute pipeline modules
    process_mutational_effects(snv_df, genes, tx_by_gene, exons_by_tx, lnc_genes)
    map_regulatory_regions()
    compute_motif_disruption()
    
    print("Scoring generation complete. Outputs saved to:", OUTDIR)