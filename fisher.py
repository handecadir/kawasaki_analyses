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
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    from scipy.stats import fisher_exact
    return hl, mo, pd


@app.cell
def _(hl):
    # Load the MatrixTable 
    filtered_missense = hl.read_matrix_table('missense.mt')
    filtered_lof = hl.read_matrix_table('lof_variants.mt')

    print("Missense entry fields:", filtered_missense.entry)
    print("LoF entry fields:", filtered_lof.entry)
    return filtered_lof, filtered_missense


@app.cell
def _(filtered_lof, filtered_missense):
    # Merge the MatrixTables
    combined_mt = filtered_missense.union_rows(filtered_lof)  # Call as a MatrixTable method

    # Get the existing IDs as strings
    # Adjust IDs to match metadata (remove the '-DNA' suffix)
    combined_mt = combined_mt.key_cols_by(s_corrected=combined_mt.s.replace('-DNA', ''))

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
def _(final_mt_with_case_control):
    final_mt_with_case_control.cols().show(192)
    return


@app.cell
def _(mo):
    mo.md(
        r"""
    --- IDs of Outliers and Mismatched Gender ---
    'M_276', 'M_277', 'M_NG3215-1', 'NG3253-1', M_36
    """
    )
    return


@app.cell
def _(final_mt_with_case_control, hl):
    # List of outliers and samples with incorrect sex
    outliers_and_wrong_sex = ['M_276', 'M_277', 'M_NG3215-1', 'NG3253-1', 'M_36']

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
    cleaned_mt.cols().show(187)
    return


@app.cell
def _(cleaned_mt):
    cleaned_mt.describe()
    return


@app.cell
def _(cleaned_mt, hl):
    # Count total cases and controls
    total_cases = cleaned_mt.aggregate_cols(hl.agg.count_where(cleaned_mt.case_control_status == 'case'))
    total_controls = cleaned_mt.aggregate_cols(hl.agg.count_where(cleaned_mt.case_control_status == 'control'))
    print(f"Total number of cases: {total_cases}")
    print(f"Total number of controls: {total_controls}")

    # Step: Count variant numbers per gene for cases and controls 
    # Group rows by gene symbol and aggregate over entries
    entries_table = cleaned_mt.entries().select('gene_symbol', 'case_control_status', 'GT')
    gene_counts_table = entries_table.group_by(entries_table.gene_symbol).aggregate(
        n_cases_variants = hl.agg.count_where(
            (entries_table.case_control_status == 'case') &
            hl.is_defined(entries_table.GT) &
            entries_table.GT.is_non_ref()
        ),
        n_controls_variants = hl.agg.count_where(
            (entries_table.case_control_status == 'control') &
            hl.is_defined(entries_table.GT) &
            entries_table.GT.is_non_ref()
        )
    )

    # Count unique variants per gene
    rows_table = cleaned_mt.rows()
    variants_per_gene = rows_table.group_by(rows_table.gene_symbol).aggregate(
        n_unique_variants = hl.agg.count()
    )

    # Set keys for joining
    variants_per_gene = variants_per_gene.key_by('gene_symbol')
    gene_counts_table = gene_counts_table.key_by('gene_symbol')

    # Add variant counts to gene_counts_table
    gene_counts_table = gene_counts_table.annotate(
        n_unique_variants = variants_per_gene[gene_counts_table.gene_symbol].n_unique_variants
    )

    gene_counts_table.show(5)

    # Perform gene burden analysis using Fisher's exact test
    # Since gene_counts_table is a Table, we can directly use annotate
    gene_results = gene_counts_table.annotate(
        fisher_result = hl.fisher_exact_test(
            hl.int32(gene_counts_table.n_cases_variants),
            hl.int32((total_cases * gene_counts_table.n_unique_variants) - gene_counts_table.n_cases_variants),
            hl.int32(gene_counts_table.n_controls_variants),
            hl.int32((total_controls * gene_counts_table.n_unique_variants) - gene_counts_table.n_controls_variants)
        )
    )

    # Sort results and display 
    gene_burden_results = gene_results.order_by(gene_results.fisher_result.p_value)
    print("\nTop 1000 significant gene burden analysis results:")
    gene_burden_results.show(1000)

    return gene_burden_results, gene_counts_table


@app.cell
def _(gene_burden_results):
    gene_burden_results.show(1000)
    return


@app.cell
def _(gene_counts_table):
    gene_counts_table.show(1000)
    return


@app.cell
def _(gene_burden_results, hl):
    # Bonferroni correction (FWER control): p_adj = min(1, p * m)
    n_tests = gene_burden_results.count()
    adjusted = gene_burden_results.annotate(
        bonferroni_p = hl.min(1.0, gene_burden_results.fisher_result.p_value * n_tests)
    )

    adjusted = adjusted.order_by(adjusted.bonferroni_p)

    print(f"Total number of tests (Bonferroni m): {n_tests}")

    adjusted.show(100)

    # Save the Bonferroni-adjusted results to a TSV file
    adjusted.export('bonferroni_adjusted_results.tsv')

    return


@app.cell
def _(combined_mt, hl):
    et = combined_mt.entries()
    et = et.key_by()
    et = et.select(
            locus = et.locus,
            alleles = et.alleles,
            gene_symbol = et.gene_symbol,
            sample_id = et.s,
            sample_id_corrected = et.s_corrected,
            genotype = hl.or_missing(
                hl.is_defined(et.GT),
                hl.case()
                  .when(et.GT.is_hom_ref(), '0/0')
                  .when(et.GT.is_het(), '0/1')
                  .when(et.GT.is_hom_var(), '1/1')
                  .default('NA')
            )
        )
    output_path = 'combined_genotypes.tsv.bgz'
    et.export(output_path)
    print(f"Exported entries to {output_path}")
    return


@app.cell
def _(gene_burden_results, hl):
    p = hl.plot.qq(gene_burden_results.fisher_result.p_value)
    hl.plot.show(p)
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
