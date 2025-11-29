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
    import requests
    import re
    import genebe as gnb
    import odf
    return hl, pd


@app.cell
def _(hl):
    snp= hl.read_matrix_table('kawasaki_filtered_outliers.mt')
    return (snp,)


@app.cell
def _(hl, pd, snp):
    target_df = pd.read_excel("variant.ods")

    target_df = target_df.rename(columns={
        'Chromosome/scaffold name': 'chrom',
        'Chromosome/scaffold position start (bp)': 'pos_start',
        'Chromosome/scaffold position end (bp)': 'pos_end'
    })

    ht_targets = hl.Table.from_pandas(target_df)

    ht_targets = ht_targets.annotate(
        chrom = 'chr' + hl.str(ht_targets.chrom),  # 16 -> chr16
        locus = hl.locus('chr' + hl.str(ht_targets.chrom), ht_targets.pos_start, reference_genome='GRCh38')
    )

    # Assign keys
    ht_targets = ht_targets.key_by('locus')

    # Key the SNP MatrixTable by locus
    snp_final = snp.key_rows_by('locus')

    # Filtering: keep only loci present in the Excel file
    filtered_mt = snp_final.filter_rows(
        hl.is_defined(ht_targets[snp_final.locus])
    )

    # Check
    print("Original number of variants:", snp_final.count_rows())
    print("Filtered number of variants:", filtered_mt.count_rows())

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

    ht_flat = filtered_mt.entries()

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

    ht_flat = ht_flat.select(
        'Sample_ID', 'case_control', 'GT', 'DP', 'AD', 'rsid', 'alleles',
        'Record_No','Case_Number','Class','CAA_status','Diagnostic_Age_Status',
        'Sequelae_Status','Family_History_Status','Sex','Consanguineous_marriage_status',
        'Degree_of_CM','Age_at_diagnosis_years','Family_history_of_CHD_status',
        'KD_in_siblings_status'
    )

    df_final = ht_flat.to_pandas()
    print(df_final.head(10))
    return (df_final,)


@app.cell
def _(df_final):

    df_final.to_csv('filtered_variants_meta.csv', index=False)
    return


if __name__ == "__main__":
    app.run()
