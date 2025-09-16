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
def _(hl, snp_result):

    rsid_list = ["rs16921209","rs17076896","rs7922552","rs17782904","rs11793049","rs12210919","rs12900413","rs10127456","rs12211370","rs1842579","rs16921209","rs7922552","rs17076896","rs28493229","rs2290692","rs113420705","rs72689236","rs1801274","rs4813003","rs1569723","rs2736340","rs2254546","rs2618476","rs2233152","rs22833188","rs2833195","rs341058","rs767007","rs2517892","rs3118470","rs28493229","rs113420705","rs3755724","rs1800872","rs16944","rs1143627","rs5050","rs7637803","rs1412125","rs899162","rs12041331","rs3789065","rs1568657","rs10129255","rs7604693","rs527409","rs7199343","rs17531088","rs1801274","rs28493229","rs2233152","rs2857151","rs4813003","rs2130392","rs2254546","rs4894410","rs9290065","rs17667932","rs1569723","rs7656244","rs2736340","rs1801274","rs12037447","rs146732504","rs151078858","rs55723436","rs6094136","rs139662037","rs7124405","rs2662865","rs202207863","rs9267431","rs12477499","rs148434007","rs201067154","rs1013532","rs117650853","rs201236531","rs3745213","rs6671847","rs111487401","rs1681087","rs4748329","rs148523506","rs2720378","rs35393613","rs2720378","rs2857602","rs2736340","rs4774175","rs2849322","rs1883832","rs113420414","rs7963257","rs407934","rs1873212","rs1778477","rs1264516","rs3129960","rs7775228","rs2071473","rs79523539","rs3743841","rs73379591","rs10508313"]

    snp_filtered = snp_result.filter_rows(
        hl.any(
            lambda x: hl.literal(rsid_list).contains(x),
            snp_result.CSQ_struct.Existing_variation.split(',')
        )
    )

    return (snp_filtered,)


@app.cell
def _(snp_filtered):
    snp_filtered.rows().show(100)
    return


if __name__ == "__main__":
    app.run()
