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
    filtered_missense = hl.read_matrix_table('missense.mt')
    filtered_lof = hl.read_matrix_table('lof_variants.mt')

    print("Missense entry fields:", filtered_missense.entry)
    print("LoF entry fields:", filtered_lof.entry)
    return filtered_lof, filtered_missense


@app.cell
def _(filtered_lof, filtered_missense):
    # Merge the MatrixTables
    combined_mt = filtered_missense.union_rows(filtered_lof)  
    return (combined_mt,)


@app.cell
def _(combined_mt):
    combined_mt.describe()
    return


@app.cell
def _(combined_mt, hl):
    # Count total cases and controls
    total_cases = combined_mt.aggregate_cols(hl.agg.count_where(combined_mt.case_control_status == 'case'))
    total_controls = combined_mt.aggregate_cols(hl.agg.count_where(combined_mt.case_control_status == 'control'))
    print(f"Total number of cases: {total_cases}")
    print(f"Total number of controls: {total_controls}")


    entries_table = combined_mt.entries().select('gene_symbol', 'case_control_status', 'GT')
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
    rows_table = combined_mt.rows()
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
def _(gene_burden_results):
    gene_burden_results.describe()
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
    adjusted.export('bonferroni_adjusted_results_lofmis.tsv')
    return


@app.cell
def _():
    from adjustText import adjust_text
    return (adjust_text,)


@app.cell
def _(adjust_text, gene_burden_results, np, plt, sns):
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
    # (df.index + 1) / (n + 1) sıralı P-değerleri için beklenen olasılıktır
    df['expected_logp'] = -np.log10((df.index + 1) / (n + 1))

    # Kontrol çıktısı
    print(f"Total number of genes: {n}")
    print(df.head())

    plt.figure(figsize=(10, 8), dpi=150)

    # --- 1. Dağılım Grafiği (Scatter Plot) ---
    sns.scatterplot(
        data=df, 
        x='expected_logp', 
        y='observed_logp', 
        color='navy', 
        alpha=0.6, 
        s=50, 
        edgecolor=None
    )

    # --- 2. Null Hipotezi Çizgisi ---
    max_val = max(df['expected_logp'].max(), df['observed_logp'].max())
    plt.plot([0, max_val], [0, max_val], color='red', linestyle='--', linewidth=1.5, label='Null Hypothesis')

    # --- 3. Gelişmiş Etiketleme (adjustText Kullanımı) ---

    # Gösterilmesini istediğiniz gen sayısını artırın (Örn: En düşük P-değerine sahip 20 gen)
    num_top_genes = 20 
    top_genes = df.head(num_top_genes)

    texts = []
    for i, row in top_genes.iterrows():
        # Metin nesnelerini oluşturun ancak konumunu otomatik ayarlamayı adjust_text'e bırakın
        texts.append(plt.text(
            x=row['expected_logp'],
            y=row['observed_logp'],      
            s=row['gene'],              
            fontsize=10,
            fontweight='bold',
            color='black',
            ha='center',                 
            va='center'                 
        ))

    # adjust_text'i çağırarak metinlerin çakışmasını engelleyin
    # arrowprops: Metinleri noktalarına bağlamak için oklar kullanır
    adjust_text(texts, 
                arrowprops=dict(arrowstyle="-", color='black', lw=0.5), 
                force_points=(0.2, 0.5), # Noktaları metinden biraz uzak tutmaya zorla
                force_text=(0.5, 0.5))   # Metinleri birbirlerinden uzak tutmaya zorla

    # --- 4. Başlık ve Eksen Etiketleri ---
    plt.title("Q-Q Plot: Gene Burden Analysis LoF", fontsize=14, fontweight='bold')
    plt.xlabel(r"Expected $-log_{10}(p)$", fontsize=12)
    plt.ylabel(r"Observed $-log_{10}(p)$", fontsize=12)
    plt.legend(loc='upper left')
    plt.grid(True, linestyle=':', alpha=0.4)

    # --- 5. SVG Kaydetme ---
    plt.savefig('lof_Q-Q Plot_Gene_Burden.svg', format='svg', bbox_inches='tight')

    plt.show()
    return


if __name__ == "__main__":
    app.run()
