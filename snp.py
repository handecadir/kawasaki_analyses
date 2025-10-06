import marimo

__generated_with = "0.14.9"
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
    import requests
    import re
    import genebe as gnb
    return hl, mo, pd


@app.cell
def _(hl):
    snp= hl.read_matrix_table('kawasaki_filtered.mt')

    snp = snp.key_cols_by(s_corrected=snp.s.replace('-DNA', ''))
    return (snp,)


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
def _(hl, meta_data_df, pd, snp):
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

    final_mt_with_pheno = snp.annotate_cols(
        pheno = pheno_table[snp.s_corrected]  # Use corrected IDs
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
def _(final_mt_with_case_control):
    final_mt_with_case_control.cols().show(192)
    return


@app.cell
def _(mo):
    mo.md(
        r"""
    --- IDs of Outliers and Mismatched Gender ---
    'M_276', 'M_277', 'M_NG3215-1', 'NG3253-1', M_36 ,'NG2605-1
    """
    )
    return


@app.cell
def _(final_mt_with_case_control, hl):
    # List of outliers and samples with incorrect sex
    outliers_and_wrong_sex = ['M_276', 'M_277', 'M_NG3215-1', 'NG3253-1', 'M_36', 'NG2605-1']

    # Convert the list to a Hail literal set
    remove_set = hl.literal(set(outliers_and_wrong_sex))

    # Remove these samples from final_mt_with_case_control
    cleaned_mt = final_mt_with_case_control.filter_cols(
        ~remove_set.contains(final_mt_with_case_control.s_corrected)
    )

    print("Number of columns before removal:", final_mt_with_case_control.count_cols())
    print("Number of columns after removal:", cleaned_mt.count_cols())

    return (cleaned_mt,)


@app.cell
def _(cleaned_mt):
    cleaned_mt.describe()
    return


@app.cell
def _(cleaned_mt, hl):
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

    snp_result = cleaned_mt.annotate_rows(
        CSQ_struct = hl.struct(**{
            field: cleaned_mt.info.CSQ[0].split('\|')[i] for i, field in enumerate(csq_fields)
        })
    )

    snp_result.CSQ_struct.show(5)
    return (snp_result,)


@app.cell
def _(hl, pd, snp_result):

    import odf
    # .ods dosyasını oku
    target_df = pd.read_excel("variant.ods")

    # Kolonları temizle
    target_df = target_df.rename(columns={
        'Chromosome/scaffold name': 'chrom',
        'Chromosome/scaffold position start (bp)': 'pos_start',
        'Chromosome/scaffold position end (bp)': 'pos_end'
    })

    # Hail Table'a çevir
    ht_targets = hl.Table.from_pandas(target_df)

    # int -> str ve 'chr' ekle
    ht_targets = ht_targets.annotate(
        chrom = 'chr' + hl.str(ht_targets.chrom),  # 16 -> chr16
        locus = hl.locus('chr' + hl.str(ht_targets.chrom), ht_targets.pos_start, reference_genome='GRCh38')
    )



    # Key atama
    ht_targets = ht_targets.key_by('locus')

    # SNP MatrixTable'i locus üzerinden key'le
    snp_final = snp_result.key_rows_by('locus')

    # Filtreleme: sadece Excel'deki loci'leri tut
    filtered_mt = snp_final.filter_rows(
        hl.is_defined(ht_targets[snp_final.locus])
    )

    # Kontrol
    print("Orijinal varyant sayısı:", snp_final.count_rows())
    print("Filtrelenmiş varyant sayısı:", filtered_mt.count_rows())

    filtered_mt.rows().show(10)
    return (filtered_mt,)


@app.cell
def _(filtered_mt):
    filtered_mt.rows().show(10)

    return


@app.cell
def _(filtered_mt):
    filtered_mt.rows().show(5) 
    filtered_mt.cols().show(5) 

    return


@app.cell
def _(filtered_mt):
    # 1️⃣ Entry-level table oluştur
    ht_flat = filtered_mt.entries()

    # 2️⃣ Sample ve meta data ekle (locus zaten key, eklemeye gerek yok)
    ht_flat = ht_flat.annotate(
        Sample_ID = ht_flat.s_corrected,
        case_control = ht_flat.case_control_status,
        Record_No = ht_flat.pheno.Record_No,
        Case_Number = ht_flat.pheno.Case_Number,
        Class = ht_flat.pheno.Class,
        CAA_status = ht_flat.pheno.CAA_status,
        Diagnostic_Age_Status = ht_flat.pheno.Diagnostic_Age_Status,
        Sequelae_Status = ht_flat.pheno.Sequelae_Status,
        Family_History_Status = ht_flat.pheno.Family_History_Status,
        Sex = ht_flat.pheno.Sex,
        Consanguineous_marriage_status = ht_flat.pheno.Consanguineous_marriage_status,
        Degree_of_CM = ht_flat.pheno.Degree_of_CM,
        Age_at_diagnosis_years = ht_flat.pheno.Age_at_diagnosis_years,
        Family_history_of_CHD_status = ht_flat.pheno.Family_history_of_CHD_status,
        KD_in_siblings_status = ht_flat.pheno.KD_in_siblings_status,
        rsid = ht_flat.CSQ_struct.Existing_variation,
        alleles = ht_flat.alleles,
        GT = ht_flat.GT,
        DP = ht_flat.DP,
        AD = ht_flat.AD
    )

    # 3️⃣ Select yaparken locus'u yazmaya gerek yok
    ht_flat = ht_flat.select(
        'Sample_ID', 'case_control', 'GT', 'DP', 'AD', 'rsid', 'alleles',
        'Record_No','Case_Number','Class','CAA_status','Diagnostic_Age_Status',
        'Sequelae_Status','Family_History_Status','Sex','Consanguineous_marriage_status',
        'Degree_of_CM','Age_at_diagnosis_years','Family_history_of_CHD_status',
        'KD_in_siblings_status'
    )

    # 4️⃣ Pandas’a çevir
    df_final = ht_flat.to_pandas()
    print(df_final.head(10))


    return (df_final,)


@app.cell
def _(df_final):
    # df_final Pandas dataframe’i CSV olarak kaydet
    df_final.to_csv('filtered_variants_meta.csv', index=False)

    return


if __name__ == "__main__":
    app.run()
