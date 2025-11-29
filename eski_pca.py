import marimo

__generated_with = "0.17.7"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import statsmodels.api as sm
    import os 
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    import numpy as np
    import hail as hl

    hl.init()
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    from scipy.stats import fisher_exact
    import requests
    import re
    return PCA, hl, np, pd, plt, sns


@app.cell
def _(hl):
    from pprint import pprint

    from hail.plot import show

    hl.plot.output_notebook()
    return


@app.cell
def _(hl):
    mt = hl.import_vcf('merge.vep.vcf.gz',reference_genome= "GRCh38", force_bgz=True, array_elements_required=False)
    return (mt,)


@app.cell
def _(mt):
    mt.s.show(5)
    return


@app.cell
def _(mt):
    mt.GT.show()
    return


@app.cell
def _(mt):
    mt.describe()
    return


@app.cell
def _(mt):
    # Get the existing IDs as strings
    # Adjust IDs to match metadata (remove the '-DNA' suffix)
    combined_mt = mt.key_cols_by(s_corrected=mt.s.replace('-DNA', ''))

    # Check the result
    print("\n--- Combined MatrixTable ---")
    print("Number of variants and samples:", combined_mt.count())
    combined_mt.describe()
    return (combined_mt,)


@app.cell
def _(combined_mt):
    combined_mt.describe()
    return


@app.cell
def _(pd):
    #Read the TSV file
    meta_data_df = pd.read_csv("meta.tsv", sep="\t", header=0)

    print("--- Original Column Names ---")
    print(meta_data_df.columns.tolist())

    #Clean column names and make them more usable
    cleaned_columns_step1 = meta_data_df.columns.str.replace(r'[\(\)\n]', '', regex=True).str.strip().str.replace(' ', '_').str.replace('__', '_')
    meta_data_df.columns = cleaned_columns_step1

    meta_data_df = meta_data_df.rename(columns={
    'No': 'Record_No',
    'Case_Number': 'Case_Number',
    'Class': 'Class',
    'CAA': 'CAA_status',
    'Diagnostic_Age': 'Diagnostic_Age_Status',
    'Sequelae': 'Sequelae_Status',
    'Family_History': 'Family_History_Status',
    'Sex': 'Sex',
    'Consanguineous_marriage': 'Consanguineous_marriage_status',
    'Degree_of_C.M.': 'Degree_of_CM',
    'Age_at_diagnosis_year': 'Age_at_diagnosis_years',
    'Family_history_of_CHD': 'Family_history_of_CHD_status',
    'KD_in_siblings': 'KD_in_siblings_status'
    })

    print("\n--- Cleaned Column Names (Final) ---")
    print(meta_data_df.columns.tolist())

    #Format IDs in the 'Case_Number' column by adding 'M_' prefix
    id_column_name = 'Case_Number'
    meta_data_df['formatted_sample_id'] = 'M_' + meta_data_df[id_column_name].astype(str)

    print("\n--- Meta Data IDs Formatted with 'M_' Prefix (First 5 Rows) ---")
    print(meta_data_df[['Case_Number', 'formatted_sample_id']].head())
    return (meta_data_df,)


@app.cell
def _(combined_mt, hl, meta_data_df, pd):
    print("--- Detailed Check of Data Types and Values ---")
    print("\nDataFrame Shape:", meta_data_df.shape)
    print("\nColumn Names:")
    for col in meta_data_df.columns:
        print(f"- {col}")  
    print("\n--- Detailed Analysis for Each Column ---")
    for col in meta_data_df.columns:
        print(f"\n{col}:")  
        print(f" Data type: {meta_data_df[col].dtype}")  
        print(f" Missing values: {meta_data_df[col].isna().sum()}")  
        print(f" Unique values: {meta_data_df[col].nunique()}")  
        print(f" First 10 values: {meta_data_df[col].value_counts().head(10).to_dict()}")  
        # Additional statistics for numeric columns
        if pd.api.types.is_numeric_dtype(meta_data_df[col]):  
            print(f" Min: {meta_data_df[col].min()}")  
            print(f" Max: {meta_data_df[col].max()}")  
            print(f" Mean: {meta_data_df[col].mean():.2f}")  

    # Only fill missing values in the 'Class' column with 'nan', leave others untouched
    print("--- Adjusting Only the 'Class' Column ---")
    if 'Class' in meta_data_df.columns:
        print("Adjusting 'Class' column...")
        print(f"Original values: {meta_data_df['Class'].value_counts().head()}")
        # Fill missing values with 'nan'
        meta_data_df['Class'] = meta_data_df['Class'].fillna('nan')
        print(f"Values after adjustment: {meta_data_df['Class'].value_counts().head()}")
        print("Other columns remain unchanged.")

    # Convert the prepared Pandas DataFrame to a Hail Table (without 'types' argument)
    pheno_table = hl.Table.from_pandas(meta_data_df, key='formatted_sample_id')

    print("\n--- Hail Phenotype Table Created and Types Checked (First 5 Rows) ---")
    pheno_table.show(192)

    final_mt_with_pheno = combined_mt.annotate_cols(
        pheno = pheno_table[combined_mt.s_corrected]  # Use corrected IDs
    )

    # assign 'case' and 'control' labels
    final_mt_with_case_control = final_mt_with_pheno.annotate_cols(
        case_control_status = hl.if_else(
            hl.is_defined(final_mt_with_pheno.pheno),
            'case',
            'control'
        )
    )
    return (final_mt_with_case_control,)


@app.cell
def _(final_mt_with_case_control, hl):
    # Since they are in the VEP CSQ, we defined them and saved separately
    csq_fields = [
        "Allele","Consequence","IMPACT","SYMBOL","Gene","Feature_type","Feature","BIOTYPE",
        "EXON","INTRON","HGVSc","HGVSp","cDNA_position","CDS_position","Protein_position",
        "Amino_acids","Codons","Existing_variation","ALLELE_NUM","DISTANCE","STRAND","FLAGS",
        "VARIANT_CLASS","MINIMISED","SYMBOL_SOURCE","HGNC_ID","CANONICAL","MANE","MANE_SELECT",
        "MANE_PLUS_CLINICAL","TSL","APPRIS","CCDS","ENSP","SWISSPROT","TREMBL","UNIPARC",
        "UNIPROT_ISOFORM","GENE_PHENO","SIFT","PolyPhen","DOMAINS","miRNA","HGVS_OFFSET",
        "AF","AFR_AF","AMR_AF","EAS_AF","EUR_AF","SAS_AF","gnomADe_AF","gnomADe_AFR_AF",
        "gnomADe_AMR_AF","gnomADe_ASJ_AF","gnomADe_EAS_AF","gnomADe_FIN_AF","gnomADe_MID_AF",
        "gnomADe_NFE_AF","gnomADe_REMAINING_AF","gnomADe_SAS_AF","gnomADg_AF","gnomADg_AFR_AF",
        "gnomADg_AMI_AF","gnomADg_AMR_AF","gnomADg_ASJ_AF","gnomADg_EAS_AF","gnomADg_FIN_AF",
        "gnomADg_MID_AF","gnomADg_NFE_AF","gnomADg_REMAINING_AF","gnomADg_SAS_AF","MAX_AF",
        "MAX_AF_POPS","CLIN_SIG","SOMATIC","PHENO","PUBMED","MOTIF_NAME","MOTIF_POS",
        "HIGH_INF_POS","MOTIF_SCORE_CHANGE","TRANSCRIPTION_FACTORS","am_class","am_pathogenicity"
    ]

    cleaned_mt_csq = final_mt_with_case_control.annotate_rows(
        CSQ_struct = hl.struct(**{
            field: final_mt_with_case_control.info.CSQ[0].split('\|')[i] for i, field in enumerate(csq_fields)
        })
    )

    cleaned_mt_csq.CSQ_struct.show(5)
    return (cleaned_mt_csq,)


@app.cell
def _(cleaned_mt_csq, hl):
    cleaned_mt_csq_filtered = cleaned_mt_csq.filter_rows(hl.len(cleaned_mt_csq.filters) == 0)
    return (cleaned_mt_csq_filtered,)


@app.cell
def _(hl):
    def filter_genotypes(
        cleaned_mt_csq_filtered: hl.MatrixTable,
        min_vaf: float = 0.15,
        min_dp: int = 15
    ) -> hl.MatrixTable:
        """Filters genotypes in a MatrixTable based on VAF and read depth.

        This function calculates the Variant Allele Frequency (VAF) for each
        genotype and then sets any genotype not meeting the min_vaf and
        min_dp thresholds to missing (NA).

        Args:
            mt: The input Hail MatrixTable.
            min_vaf: The minimum VAF required to keep a genotype.
            min_dp: The minimum read depth (DP) required to keep a genotype.

        Returns:
            A new MatrixTable with low-quality genotypes set to missing.
        """
        # Calculate VAF per genotype, handling cases where DP might be 0.
        mt_with_vaf = cleaned_mt_csq_filtered.annotate_entries(
            VAF=hl.if_else(cleaned_mt_csq_filtered.DP > 0, cleaned_mt_csq_filtered.AD[1] / cleaned_mt_csq_filtered.DP, 0.0)
        )

       # New filtering condition
        condition = (
            (mt_with_vaf.GT.is_hom_ref())  # keep 0/0 always
            | ( (mt_with_vaf.VAF > min_vaf) & (mt_with_vaf.DP > min_dp) )
        )

        filtered_mt = mt_with_vaf.filter_entries(condition)
        return filtered_mt


    def remove_uncalled_variants(cleaned_mt_csq_filtered: hl.MatrixTable) -> hl.MatrixTable:
        """Removes variants that have no called genotypes.

        This should be run after filtering entries to clean up the dataset
        by removing variants that are no longer polymorphic in any sample.

        Args:
            mt: The MatrixTable to clean (e.g., after running filter_genotypes).

        Returns:
            A new MatrixTable with uncalled variants removed.
        """
        # Run variant QC to get stats like the number of called genotypes.
        mt_with_qc = hl.variant_qc(cleaned_mt_csq_filtered)

        # Filter rows (variants) where the number of called genotypes is greater than 0.
        final_mt = mt_with_qc.filter_rows(mt_with_qc.variant_qc.n_called > 0)

        return final_mt
    return filter_genotypes, remove_uncalled_variants


@app.cell
def _(cleaned_mt_csq_filtered, filter_genotypes):
    filtered_mt = filter_genotypes(cleaned_mt_csq_filtered)
    return (filtered_mt,)


@app.cell
def _(hl):
    from bokeh.layouts import gridplot

    def plot_variant_qc_metrics(mt_with_qc: hl.MatrixTable):
        """
        Generates and displays a grid of plots for key variant QC metrics
        using the exact schema provided.

        Args:
            mt_with_qc: A MatrixTable with a `variant_qc` row struct.
        """
        # 1. Create individual histogram plots for each metric
        p1 = hl.plot.histogram(mt_with_qc.variant_qc.call_rate,
                               legend='Call Rate',
                               title='Variant Call Rate')

        p2 = hl.plot.histogram(mt_with_qc.variant_qc.AF[1],
                               legend='Allele Frequency',
                               title='Allele Frequency Spectrum')

        # Corrected path for mean DP
        p3 = hl.plot.histogram(hl.log10(mt_with_qc.variant_qc.dp_stats.mean),
                               legend='Mean Depth (log10)',
                               title='Mean Read Depth')

        # Corrected path for mean GQ
        p4 = hl.plot.histogram(mt_with_qc.variant_qc.gq_stats.mean,
                               legend='Mean GQ',
                               title='Mean Genotype Quality')

        # 2. Arrange the plots in a grid and show them
        plot_grid = gridplot([[p1, p2], [p3, p4]])
        hl.plot.show(plot_grid)
    return (plot_variant_qc_metrics,)


@app.cell
def _(cleaned_mt_csq_filtered, hl, plot_variant_qc_metrics):
    plot_variant_qc_metrics(hl.variant_qc(cleaned_mt_csq_filtered))
    return


@app.cell
def _(filtered_mt, hl, plot_variant_qc_metrics):
    plot_variant_qc_metrics(hl.variant_qc(filtered_mt))
    return


@app.cell
def _(filtered_mt):
    filtered_mt.count_rows()
    return


@app.cell
def _(filtered_mt, remove_uncalled_variants):
    final_mt= remove_uncalled_variants(filtered_mt)
    return (final_mt,)


@app.cell
def _(final_mt):
    final_mt.count_rows()
    return


@app.cell
def _(PCA, final_mt, hl, np, pd, plt, sns):
    # Compute PCA (first 2 components)
    eigenvalues, scores, loadings = hl.hwe_normalized_pca(final_mt.GT, k=2)

    # Add case/control information to the scores Table
    scores = scores.annotate(case_control_status=final_mt.cols()[scores.s_corrected].case_control_status)

    # Convert to Pandas
    scores_pd = scores.to_pandas()

    # Extract PC1 and PC2 columns
    scores_pd['PC1'] = scores_pd['scores'].apply(lambda x: x[0])
    scores_pd['PC2'] = scores_pd['scores'].apply(lambda x: x[1])

    # Scatter plot for PCA
    plt.figure(figsize=(8,6))
    for label, color in [("case", "red"), ("control", "blue")]:
        subset = scores_pd[scores_pd['case_control_status'] == label]
        plt.scatter(subset['PC1'], subset['PC2'], label=label.capitalize(), alpha=0.7, c=color)

    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.legend()
    plt.title("PCA: Case vs Control")
    plt.show()


    # Create a Case-only MatrixTable
    case_mt = final_mt.filter_cols(final_mt.case_control_status == "case")

    # Convert column metadata into a Table
    case_cols_ht = case_mt.cols()

    # Select phenotype fields of interest
    case_cols_ht = case_cols_ht.select(
        **{f: case_cols_ht.pheno[f] for f in [
            'CAA_status',
            'Diagnostic_Age_Status',
            'Sequelae_Status',
            'Family_History_Status',
            'Consanguineous_marriage_status',
            'Degree_of_CM',
            'Family_history_of_CHD_status',
            'KD_in_siblings_status'
        ]}
    )

    # Convert to Pandas DataFrame (GERÇEK VERİ)
    df = case_cols_ht.to_pandas()

    # Feature list
    features = [
        'CAA_status', 'Diagnostic_Age_Status', 'Sequelae_Status',
        'Family_History_Status', 'Consanguineous_marriage_status',
        'Degree_of_CM', 'Family_history_of_CHD_status', 'KD_in_siblings_status'
    ]

    pca_data = df[features].astype(int)

    # PCA
    pca = PCA(n_components=2)
    principal_components = pca.fit_transform(pca_data)

    # PCA result DataFrame
    pca_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'])

    # Unique points
    unique_points = pca_df.drop_duplicates()
    print(f"\nNumber of unique PCA points: {len(unique_points)}")
    print(f"Total number of PCA points: {len(pca_df)}")

    # Jitter adding
    np.random.seed(42)
    jitter_amount = 0.05 
    pca_df['PC1_jittered'] = pca_df['PC1'] + np.random.normal(0, jitter_amount, size=len(pca_df))
    pca_df['PC2_jittered'] = pca_df['PC2'] + np.random.normal(0, jitter_amount, size=len(pca_df))

    # Scatter plot
    plt.figure(figsize=(10, 8), facecolor='#f5f5f5')
    sns.scatterplot(x='PC1_jittered', y='PC2_jittered', data=pca_df, s=70, alpha=0.8)
    plt.title('PCA Distribution of Cases (Points Separated with Jittering)')
    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.show()

    # For every featues scatter
    for feature in features:
        plt.figure(figsize=(10, 8), facecolor='#f5f5f5')
        sns.scatterplot(x='PC1_jittered', y='PC2_jittered', hue=df[feature], data=pca_df, s=70, alpha=0.8, palette='viridis')
        plt.title(f'PCA by {feature} Status (All Patients Visible)', fontsize=16, fontweight='bold', pad=20)
        plt.xlabel('Principal Component 1', fontsize=12)
        plt.ylabel('Principal Component 2', fontsize=12)
        plt.grid(True, linestyle='--', alpha=0.6)
        plt.axhline(0, color='grey', linestyle='--', linewidth=0.8)
        plt.axvline(0, color='grey', linestyle='--', linewidth=0.8)
        plt.box(on=True)
        plt.tight_layout()
        plt.show()

    # Her bir özellik için scatter plot ve kaydetme
    for feature in features:
        plt.figure(figsize=(10, 8), facecolor='#f5f5f5')
        sns.scatterplot(x='PC1_jittered', y='PC2_jittered', hue=df[feature], data=pca_df, s=70, alpha=0.8, palette='viridis')
        plt.title(f'PCA by {feature} Status (All Patients Visible)', fontsize=16, fontweight='bold', pad=20)
        plt.xlabel('Principal Component 1', fontsize=12)
        plt.ylabel('Principal Component 2', fontsize=12)
        plt.grid(True, linestyle='--', alpha=0.6)
        plt.axhline(0, color='grey', linestyle='--', linewidth=0.8)
        plt.axvline(0, color='grey', linestyle='--', linewidth=0.8)
        plt.box(on=True)
        plt.tight_layout()

        # --- Grafik kaydetme satırı burası! ---
        # Dosya adını dinamik olarak oluştur
        file_name = f'pca_plot_{feature}.png'
        plt.savefig(file_name, dpi=300, bbox_inches='tight')
        print(f'Grafik başarıyla kaydedildi: {file_name} 🎉')

        # Grafiği göster
        plt.show()
    return (scores,)


@app.cell
def _(plt, scores):
    scores_correct = scores.to_pandas()

    # Outlier detection thresholds
    pc1_threshold = -0.1
    pc2_lower_threshold = -0.075
    pc2_upper_threshold = 0.075

    # Extract PC1 and PC2 columns from the scores
    scores_correct['PC1'] = scores_correct['scores'].apply(lambda x: x[0])
    scores_correct['PC2'] = scores_correct['scores'].apply(lambda x: x[1])

    # Filter inliers based on the thresholds
    clean_scores_correct = scores_correct[
        (scores_correct['PC1'] > pc1_threshold) & 
        (scores_correct['PC2'] > pc2_lower_threshold) &
        (scores_correct['PC2'] < pc2_upper_threshold)
    ]

    print(f"Original dataset size: {scores_correct.shape[0]}")
    print(f"Cleaned dataset size: {clean_scores_correct.shape[0]}")

    # Find outliers by inverse filtering
    outliers = scores_correct[
        ~((scores_correct['PC1'] > pc1_threshold) &
          (scores_correct['PC2'] > pc2_lower_threshold) &
          (scores_correct['PC2'] < pc2_upper_threshold))
    ]

    # Create figure and axes for visualization
    plt.figure(figsize=(10, 8))

    # Plot inliers 
    plt.scatter(clean_scores_correct['PC1'], clean_scores_correct['PC2'], 
                label='Clean Data (Inliers)', color='skyblue', alpha=0.7)

    # Plot outliers 
    plt.scatter(outliers['PC1'], outliers['PC2'], 
                label='Outliers', color='red', alpha=0.7)

    # Add threshold lines for clarity
    plt.axvline(x=pc1_threshold, color='green', linestyle='--', label='PC1 Threshold')
    plt.axhline(y=pc2_lower_threshold, color='purple', linestyle='--', label='PC2 Lower Threshold')
    plt.axhline(y=pc2_upper_threshold, color='orange', linestyle='--', label='PC2 Upper Threshold')

    # Add title and axis labels for clarity
    plt.title('PCA Scores and Outlier Detection', fontsize=16)
    plt.xlabel('Principal Component 1 (PC1)', fontsize=12)
    plt.ylabel('Principal Component 2 (PC2)', fontsize=12)
    plt.legend()
    plt.grid(True)
    plt.show()
    return


if __name__ == "__main__":
    app.run()
