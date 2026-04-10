import os
import pickle
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

# ==========================================
# 1. CONFIGURATION & FILE PATHS
# ==========================================
PROJECT_FOLDER = './'
DATA_DIR = os.path.join(PROJECT_FOLDER, 'Data')
OUTDIR = os.path.join(DATA_DIR, 'VariantImpactScore')
os.makedirs(OUTDIR, exist_ok=True)

# Input files
SNV_PATH = os.path.join(DATA_DIR, 'Main_data/562_data/save_folder/master_table/pca+_447_samples_snv_v6.csv')
DICT_PATH = os.path.join(DATA_DIR, 'Driver_lncRNA/driver_lncRNA_dict.pkl')
MTV_DLNC_PATH = os.path.join(DATA_DIR, 'Driver_lncRNA/multivariate_dlncrna_esults.pkl')
MTV_SNV_PATH = os.path.join(DATA_DIR, 'Driver_lncRNA/multivariate_snv_esults.pkl')

RNASNP_PATH = os.path.join(DATA_DIR, 'RNAstructure/input/RNAsnp_combined_results.csv')
RNAPLFOLD_PATH = os.path.join(DATA_DIR, 'RNAstructure/input/RNAplfold_accessibility_changes.csv')
REG_MAP_PATH = os.path.join(PROJECT_FOLDER, 'Code/regulatory_summary/Table_S1_SNV_annotations.csv')
TF_MOTIF_PATH = os.path.join(DATA_DIR, 'TF_motif/tf_motif_disruption/per_snv_tf_calls_STRONG.csv')

# TF Families defined in the pipeline
TF_FAMILIES = ['AR', 'FOXA_FAM', 'HOXB13', 'GATA_FAM', 'NKX3-1', 
               'ETS_FAM', 'TWIST1', 'ZEB1', 'SNAI1', 'SOX4']

# ==========================================
# 2. DATA LOADING & FILTERING
# ==========================================
def load_base_data():
    print("Loading SNVs and filtering for driver lncRNAs...")
    snp_df = pd.read_csv(SNV_PATH)
    snp_df['Gene_ID_Trimmed'] = snp_df['Gene ID'].str.split('.').str[0]
    
    with open(DICT_PATH, 'rb') as f:
        loaded_dict = pickle.load(f)
        
    compd_dlnc_list = loaded_dict['compd_dlnc_list']
    sel_df = snp_df[snp_df['Gene_ID_Trimmed'].isin(compd_dlnc_list)].copy()
    
    # Assign ancestry types based on lncRNA lists
    compd_dlnc_af_sp = loaded_dict['compd_dlnc_af_sp']['Gene'].values
    compd_dlnc_eu_sp = loaded_dict['compd_dlnc_eu_sp']['Gene'].values
    
    conditions = [
        sel_df['Gene_ID_Trimmed'].isin(compd_dlnc_af_sp),
        sel_df['Gene_ID_Trimmed'].isin(compd_dlnc_eu_sp)
    ]
    sel_df['lncRNA_ancestry_type'] = np.select(conditions, ['African-specific', 'European-specific'], default='Common')
    
    # Load Clinical
    with open(MTV_SNV_PATH, "rb") as f:
        mtv_snv_dict = pickle.load(f)
        
    mtv_pos_af_snv = set(mtv_snv_dict['mtv_pos_af_snv'])
    
    # Assign SNV ancestry type
    sel_df['snv_name'] = ['snv_' + str(i) for i in sel_df.index]
    sel_df['snv_ancestry_type'] = sel_df['snv_name'].apply(
        lambda x: 'African' if x in mtv_pos_af_snv else 'Non-African'
    )
    
    return sel_df, mtv_snv_dict

# ==========================================
# 3. FEATURE MAPPING (ARS Components)
# ==========================================
def build_feature_dictionaries():
    print("Building functional annotation mapping dictionaries...")
    
    # RNAsnp
    rnasnp_df = pd.read_csv(RNASNP_PATH)
    rnasnp_df = rnasnp_df[(rnasnp_df['user_ref_match'] == True) & (rnasnp_df['status_warning'].isna())]
    rnasnp_pglobal = rnasnp_df.set_index(['gene_id', 'chrom', 'gpos', 'ref_tx', 'alt_tx'])['pvalue_global'].to_dict()
    rnasnp_plocal = rnasnp_df.set_index(['gene_id', 'chrom', 'gpos', 'ref_tx', 'alt_tx'])['pvalue_local'].to_dict()

    # RNAplfold
    rnaplfold_df = pd.read_csv(RNAPLFOLD_PATH)
    # Filter strictly to variants successfully scored in RNAsnp
    rnasnp_keys = set(rnasnp_pglobal.keys())
    rnaplfold_df['key'] = list(zip(rnaplfold_df.gene_id, rnaplfold_df.chrom, rnaplfold_df.gpos, rnaplfold_df.ref_tx, rnaplfold_df.alt_tx))
    rnaplfold_df = rnaplfold_df[rnaplfold_df['key'].isin(rnasnp_keys)]
    
    r_pu1 = rnaplfold_df.set_index('key')['Pu1_delta'].to_dict()
    r_pu8 = rnaplfold_df.set_index('key')['Pu8cover_delta'].to_dict()

    # Regulatory
    regmap = pd.read_csv(REG_MAP_PATH)
    reg_keys = ['gene_id', 'chrom', 'pos', 'ref', 'alt']
    reg_promoter = regmap.set_index(reg_keys)['in_promoter'].to_dict()
    reg_pels = regmap.set_index(reg_keys)['hit_pELS'].to_dict()
    reg_dels = regmap.set_index(reg_keys)['hit_dELS'].to_dict()
    reg_ctcf = regmap.set_index(reg_keys)['hit_CTCF'].to_dict()

    # TFs
    tf_df = pd.read_csv(TF_MOTIF_PATH)
    reg_hit_keys = set(regmap[(regmap['in_promoter']==True) | (regmap['hit_ELS']==True)].set_index(reg_keys).index)
    tf_df['key'] = list(zip(tf_df.gene_id, tf_df.chrom, tf_df.pos, tf_df.ref, tf_df.alt))
    tf_df = tf_df[tf_df['key'].isin(reg_hit_keys)]
    
    tf_strong = tf_df.set_index(['tf_family', 'gene_id', 'chrom', 'pos', 'ref', 'alt'])['strong_change'].to_dict()

    return {
        'rnasnp_pglobal': rnasnp_pglobal, 'rnasnp_plocal': rnasnp_plocal,
        'r_pu1': r_pu1, 'r_pu8': r_pu8,
        'reg_promoter': reg_promoter, 'reg_pels': reg_pels, 'reg_dels': reg_dels, 'reg_ctcf': reg_ctcf,
        'tf_strong': tf_strong
    }

def map_features_to_snvs(sel_df, dicts):
    print("Mapping features and computing ARS scores...")
    
    def get_ars_score(row):
        score = 0
        coor = (row['Gene_ID_Trimmed'], row['Chromosome'], row['POS'], row['REF'], row['ALT'])
        
        # 1-2. Structure (RNAsnp)
        pglobal = dicts['rnasnp_pglobal'].get(coor, np.nan)
        plocal = dicts['rnasnp_plocal'].get(coor, np.nan)
        if pd.notna(pglobal) and pglobal < 0.05: score += 2
        if pd.notna(plocal) and plocal < 0.05: score += 2
            
        # 3-4. Accessibility (RNAplfold)
        pu1 = dicts['r_pu1'].get(coor, np.nan)
        pu8 = dicts['r_pu8'].get(coor, np.nan)
        if pd.notna(pu1) and abs(pu1) > 0.1: score += 2
        if pd.notna(pu8) and abs(pu8) > 0.1: score += 2
            
        # 5-8. Regulation (Epigenetic Features)
        if dicts['reg_promoter'].get(coor) is True: score += 1
        if dicts['reg_pels'].get(coor) is True: score += 1
        if dicts['reg_dels'].get(coor) is True: score += 1
        if dicts['reg_ctcf'].get(coor) is True: score += 1
            
        # 9. Regulation (TF Motifs)
        for fam in TF_FAMILIES:
            tf_coor = (fam, coor[0], coor[1], coor[2], coor[3], coor[4])
            if dicts['tf_strong'].get(tf_coor) is True:
                score += 1
                
        return score

    sel_df['ARS_Score'] = sel_df.apply(get_ars_score, axis=1)
    return sel_df

# ==========================================
# 4. PLOTTING & REPORTING
# ==========================================
def generate_publication_plot(scored_df):
    print("Generating publication-ready ARS distribution plot...")
    mpl.rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": 1000,          
        "font.size": 8.5,             
        "axes.linewidth": 0.6,
        "axes.titlesize": 9.5,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.minor.visible": False,
        "ytick.minor.visible": True,
        "grid.linewidth": 0.4,
    })

    x = scored_df["ARS_Score"].astype(int)
    bins = np.arange(x.min() - 0.5, x.max() + 1.5, 1)
    counts, edges = np.histogram(x, bins=bins)
    centers = 0.5 * (edges[:-1] + edges[1:])

    fig = plt.figure(figsize=(5, 3))
    ax = plt.gca()

    ax.bar(centers, counts, width=0.9, edgecolor="black", linewidth=0.5, color="darkorchid", alpha=0.9)

    for xc, c in zip(centers, counts):
        if c >= 3:
            ax.text(xc, c + max(0.02 * counts.max(), 0.5), f"{c}", ha="center", va="bottom", fontsize=7)

    ax.set_xlabel("ARS (Accessibility, Regulation & Structure) Impact Score")
    ax.set_ylabel("Number of Somatic SNVs")
    ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins='auto'))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=6))
    ax.grid(axis="y", linestyle="-", alpha=0.25)
    ax.set_ylim(0, max(1, counts.max()) * 1.18)

    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    fig.tight_layout(pad=0.3)
    out_png = os.path.join(OUTDIR, "ars_distribution_1000dpi.png")
    fig.savefig(out_png, dpi=1000, bbox_inches="tight")
    plt.close(fig)
    print(f"Plot saved to: {out_png}")

if __name__ == "__main__":
    sel_df, mtv_dict = load_base_data()
    dicts = build_feature_dictionaries()
    scored_df = map_features_to_snvs(sel_df, dicts)
    
    # Save the combined scores
    out_csv = os.path.join(OUTDIR, 'ARS_scored_driver_variants.csv')
    scored_df.to_csv(out_csv, index=False)
    print(f"ARS Scoring Complete. Output saved to: {out_csv}")
    
    # Generate Plots
    generate_publication_plot(scored_df)