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
    rare_filtered_mis = hl.read_matrix_table('missense.mt')
    return (rare_filtered_mis,)


@app.cell
def _(hl, rare_filtered_mis):
    total_cases = rare_filtered_mis.aggregate_cols(hl.agg.count_where(rare_filtered_mis.case_control_status == 'case'))
    total_controls = rare_filtered_mis.aggregate_cols(hl.agg.count_where(rare_filtered_mis.case_control_status == 'control'))

    print(f"Total case number: {total_cases}")
    print(f"Total control number: {total_controls}")

    gene_mt = rare_filtered_mis.group_rows_by(rare_filtered_mis.gene_symbol).aggregate(
        has_variant = hl.agg.any(rare_filtered_mis.GT.is_non_ref())
    )

    gene_counts_table = gene_mt.annotate_rows(
        n_cases_variants = hl.agg.count_where((gene_mt.case_control_status == 'case') & (gene_mt.has_variant)),
        n_controls_variants = hl.agg.count_where((gene_mt.case_control_status == 'control') & (gene_mt.has_variant))
    ).rows()

    gene_results = gene_counts_table.annotate(
        fisher_result = hl.fisher_exact_test(
            hl.int32(gene_counts_table.n_cases_variants),
            hl.int32(total_cases - gene_counts_table.n_cases_variants),
            hl.int32(gene_counts_table.n_controls_variants),
            hl.int32(total_controls - gene_counts_table.n_controls_variants)
        )
    )
    gene_burden_results = gene_results.order_by(gene_results.fisher_result.p_value)

    gene_burden_results.show(100)
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
    adjusted.export('bonferroni_adjusted_results_mis.tsv')
    return (adjusted,)


@app.cell
def _(adjusted):
    adjusted.show(100)
    return


@app.cell
def _():
    from adjustText import adjust_text
    return (adjust_text,)


@app.cell
def _(adjust_text, adjusted, np, plt, sns):
    df = adjusted.select(
        gene = adjusted.gene_symbol,
        p_value = adjusted.fisher_result.p_value,
        bonferroni_p = adjusted.bonferroni_p
    ).to_pandas()

    df = df.dropna(subset=['p_value']).sort_values('p_value').reset_index(drop=True)
    n = len(df)

    df['observed_logp'] = -np.log10(df['p_value'])
    df['expected_logp'] = -np.log10((df.index + 1) / (n + 1))

    print(f"Total number of genes: {n}")

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

    significant_genes = df[df['bonferroni_p'] < 0.05]
    print(f"Number of significant genes to label: {len(significant_genes)}")

    texts = []
    for i, row in significant_genes.iterrows():
        texts.append(plt.text(
            x=row['expected_logp'], 
            y=row['observed_logp'],      
            s=row['gene'],              
            fontsize=9,                  
            fontweight='bold',
            color='black',
            ha='center',                  
            va='center'                  
        ))

    if texts:
        adjust_text(texts, 
                    arrowprops=dict(arrowstyle="-", color='black', lw=0.5), 
                    force_points=(0.2, 0.5), 
                    force_text=(0.5, 0.5))

    plt.title(" Missense Significant Genes (Bonferroni Adjusted)", fontsize=14, fontweight='bold')
    plt.xlabel(r"Expected $-log_{10}(p)$", fontsize=12)
    plt.ylabel(r"Observed $-log_{10}(p)$", fontsize=12)
    plt.legend(loc='lower right')
    plt.grid(True, linestyle=':', alpha=0.4)

    plt.savefig('mis_significant_QQ_Plot.svg', format='svg', bbox_inches='tight')
    plt.show()
    return


if __name__ == "__main__":
    app.run()
