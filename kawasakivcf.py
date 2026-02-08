import marimo

__generated_with = "0.17.7"
app = marimo.App(width="full")


@app.cell
def _():
    import marimo as mo
    import statsmodels.api as sm
    import os 
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    import numpy as np
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    from scipy.stats import fisher_exact
    import requests
    import re
    from pprint import pprint
    import hail as hl

    hl.init()
    from hail.plot import show
    return hl, mo, pd


@app.cell
def _(mo):
    mo.md(r"""
    bcftools merge -m none n7_exome_calls.vcf.gz control_exome_calls.vcf.bgz exome_calls.vcf.bgz |
    bcftools norm -m-any -O  z -o merge.vcf.gz

    Lines   total/split/realigned/skipped:	714050/46419/0/0
    """)
    return


@app.cell
def _(hl):
    mt = hl.import_vcf('merge.vep.vcf.gz',reference_genome= "GRCh38", force_bgz=True, array_elements_required=False)
    return (mt,)


@app.cell
def _(mt):
    mt.describe()
    return


@app.cell
def _(mt):
    mt.GT.show(2000)
    return


@app.cell
def _(mt):
    combined_mt = mt.key_cols_by(s_corrected=mt.s.replace('-DNA', ''))
    combined_mt.describe()
    return (combined_mt,)


@app.cell
def _(combined_mt):
    combined_mt.describe()
    return


@app.cell
def _(pd):
    meta_data_df = pd.read_csv("meta.tsv", sep="\t", header=0)
    print(meta_data_df.columns.tolist())

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

    print(meta_data_df.columns.tolist())

    id_column_name = 'Case_Number'
    meta_data_df['formatted_sample_id'] = 'M_' + meta_data_df[id_column_name].astype(str)

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

        if pd.api.types.is_numeric_dtype(meta_data_df[col]):  
            print(f" Min: {meta_data_df[col].min()}")  
            print(f" Max: {meta_data_df[col].max()}")  
            print(f" Mean: {meta_data_df[col].mean():.2f}")  


    print("--- Adjusting Only the 'Class' Column ---")
    if 'Class' in meta_data_df.columns:
        print("Adjusting 'Class' column...")
        print(f"Original values: {meta_data_df['Class'].value_counts().head()}")
        # Fill missing values with 'nan'
        meta_data_df['Class'] = meta_data_df['Class'].fillna('nan')
        print(f"Values after adjustment: {meta_data_df['Class'].value_counts().head()}")
        print("Other columns remain unchanged.")


    pheno_table = hl.Table.from_pandas(meta_data_df, key='formatted_sample_id')

    print("\n--- Hail Phenotype Table Created and Types Checked (First 5 Rows) ---")
    pheno_table.show(192)

    final_mt_with_pheno = combined_mt.annotate_cols(
        pheno = pheno_table[combined_mt.s_corrected]  
    )

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
def _(final_mt_with_case_control):
    final_mt_with_case_control.s_corrected.show(192)
    return


@app.cell
def _(final_mt_with_case_control, hl):
    # VEP CSQ
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
def _(cleaned_mt_csq):
    cleaned_mt_csq.rows().show(200)
    return


@app.cell
def _(cleaned_mt_csq, hl):
    cleaned_mt_csq_filtered = cleaned_mt_csq.filter_rows(hl.len(cleaned_mt_csq.filters) == 0)
    return (cleaned_mt_csq_filtered,)


@app.cell
def _(hl):
    new_genotype_matrix = hl.import_vcf(
        'beni_haila_yukle.split.filtered.vcf.gz',
        reference_genome="GRCh38",
        force_bgz=True,
        array_elements_required=False 
    )
    return (new_genotype_matrix,)


@app.cell
def _(cleaned_mt_csq_filtered, hl, new_genotype_matrix):
    hl.plot.output_notebook()
    old_matrix = cleaned_mt_csq_filtered.unfilter_entries()

    new_renamed = new_genotype_matrix.key_cols_by(s_corrected=new_genotype_matrix.s)

    old_matrix_temp = old_matrix.annotate_entries(
         newGT = new_renamed.index_entries(
            old_matrix.row_key, 
            old_matrix.s_corrected 
        ).GT
    )

    updated_matrix = old_matrix_temp.annotate_entries(
         GT = hl.if_else(
             (old_matrix_temp.case_control_status == "control") & 
             hl.is_defined(old_matrix_temp.newGT),
             old_matrix_temp.newGT,
             old_matrix_temp.GT
         )
    )

    final_matrix = updated_matrix
    return final_matrix, old_matrix


@app.cell
def _(final_matrix):
    final_matrix.entries().show(10)
    return


@app.cell
def _(old_matrix):
    old_matrix.describe()
    return


@app.cell
def _(final_matrix):
    final_matrix.write("kawasaki.mt", overwrite=True)
    return


@app.cell
def _(final_matrix, hl):
    hl.export_vcf(final_matrix, 'kawasaki_final.vcf.bgz') #for relatedness2 
    return


if __name__ == "__main__":
    app.run()
