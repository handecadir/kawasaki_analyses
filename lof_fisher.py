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
    return hl, np, plt, sns


@app.cell
def _(hl):
    lof_filtered= hl.read_matrix_table('lof_variants.mt')
    return (lof_filtered,)


@app.cell
def _(lof_filtered):
    lof_filtered.describe()
    return


@app.cell
def _(hl, lof_filtered):
    total_cases = lof_filtered.aggregate_cols(hl.agg.count_where(lof_filtered.case_control_status == 'case'))
    total_controls = lof_filtered.aggregate_cols(hl.agg.count_where(lof_filtered.case_control_status == 'control'))

    print(f"Toplam vaka sayısı: {total_cases}")
    print(f"Toplam kontrol sayısı: {total_controls}")



    gene_counts_mt = lof_filtered.group_rows_by(lof_filtered.gene_symbol).aggregate(
        n_cases_variants=hl.agg.count_where(lof_filtered.case_control_status == 'case'),
        n_controls_variants=hl.agg.count_where(lof_filtered.case_control_status == 'control')
    )


    gene_counts_table = gene_counts_mt.entries()

    gene_counts_table.show(5)


    gene_results = gene_counts_table.annotate(
        fisher_result = hl.fisher_exact_test(
            hl.int32(gene_counts_table.n_cases_variants),
            hl.int32(total_cases - gene_counts_table.n_cases_variants),
            hl.int32(gene_counts_table.n_controls_variants),
            hl.int32(total_controls - gene_counts_table.n_controls_variants)
        )
    )


    gene_burden_results = gene_results.order_by(gene_results.fisher_result.p_value)

    gene_burden_results.show(1000)
    return gene_burden_results, gene_counts_table


@app.cell
def _(gene_burden_results):
    gene_burden_results.show()
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
    adjusted.export('bonferroni_adjusted_results_lof.tsv')
    return


@app.cell
def _(gene_burden_results, np, plt, sns):
    df = gene_burden_results.select(
        gene = gene_burden_results.gene_symbol,
        p_value = gene_burden_results.fisher_result.p_value
    ).to_pandas()

    # Drop missing P-values and sort
    df = df.dropna(subset=['p_value']).sort_values('p_value').reset_index(drop=True)

    # Total number of tests
    n = len(df)

    # Y-Axis: Observed -log10 p-values
    df['observed_logp'] = -np.log10(df['p_value'])

    # X-Axis: Expected -log10 p-values
    df['expected_logp'] = -np.log10((df.index + 1) / (n + 1))

    # Check output
    print(f"Total number of genes: {n}")
    print(df.head())

    plt.figure(figsize=(10, 8), dpi=150)


    sns.scatterplot(
        data=df, 
        x='expected_logp', 
        y='observed_logp', 
        color='navy', 
        alpha=0.6, 
        s=50, 
        edgecolor=None
    )

    max_val = max(df['expected_logp'].max(), df['observed_logp'].max())
    plt.plot([0, max_val], [0, max_val], color='red', linestyle='--', linewidth=1.5, label='Null Hypothesis')

    top_genes = df.head(10)
    for i, row in top_genes.iterrows():
        plt.text(
            x=row['expected_logp'] + 0.2,
            y=row['observed_logp'],        
            s=row['gene'],                 
            fontsize=10,
            fontweight='bold',
            color='black',
            ha='left',                     
            va='center'                    
        )

    plt.title("Q-Q Plot: Gene Burden Analysis LoF", fontsize=14, fontweight='bold')
    plt.xlabel(r"Expected $-log_{10}(p)$", fontsize=12)
    plt.ylabel(r"Observed $-log_{10}(p)$", fontsize=12)
    plt.legend(loc='upper left')
    plt.grid(True, linestyle=':', alpha=0.4)

    plt.savefig('lof_Q-Q_Gene_Burden.png', dpi=300, bbox_inches='tight')

    plt.show()
    return


if __name__ == "__main__":
    app.run()
