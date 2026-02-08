import marimo

__generated_with = "0.17.7"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    import numpy as np
    import hail as hl
    hl.init()
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    from scipy.stats import chi2
    from matplotlib.patches import Ellipse
    return Ellipse, PCA, chi2, hl, np, pd, plt, sns


@app.cell
def _(hl):
    result = hl.read_matrix_table('kawasaki_filtered.mt')
    return (result,)


@app.cell
def _(result):
    result.describe()
    return


@app.cell
def _(result):
    result.cols().show(192)
    return


@app.cell
def _(PCA, hl, np, pd, plt, result, sns):
    # Compute PCA (first 2 components)
    eigenvalues, scores, loadings = hl.hwe_normalized_pca(result.GT, k=2)

    scores = scores.annotate(case_control_status=result.cols()[scores.s_corrected].case_control_status)

    scores_pd = scores.to_pandas()

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
    file_name_case_control = 'pca_plot_case_control.svg'
    plt.savefig(file_name_case_control, format='svg', bbox_inches='tight')
    plt.show()


    case_mt = result.filter_cols(result.case_control_status == "case")
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

    df = case_cols_ht.to_pandas()

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

    unique_points = pca_df.drop_duplicates()
    print(f"\nNumber of unique PCA points: {len(unique_points)}")
    print(f"Total number of PCA points: {len(pca_df)}")


    np.random.seed(42)
    jitter_amount = 0.05 
    pca_df['PC1_jittered'] = pca_df['PC1'] + np.random.normal(0, jitter_amount, size=len(pca_df))
    pca_df['PC2_jittered'] = pca_df['PC2'] + np.random.normal(0, jitter_amount, size=len(pca_df))


    plt.figure(figsize=(10, 8), facecolor='#f5f5f5')
    sns.scatterplot(x='PC1_jittered', y='PC2_jittered', data=pca_df, s=70, alpha=0.8)
    plt.title('PCA Distribution of Cases (Points Separated with Jittering)')
    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.show()

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


        file_name = f'pca_plot_{feature}.svg'
        plt.savefig(file_name, format='svg', bbox_inches='tight')

        plt.show()
    return scores, scores_pd


@app.cell
def _(np, plt, result, sns):
    case_mt_sex = result.filter_cols(result.case_control_status == 'case')
    df_sex = case_mt_sex.cols().flatten().to_pandas()

    sex_mapping = {'Male': 1, 'Female': 2}
    df_sex['Sex_Numeric'] = df_sex['pheno.Sex'].map(sex_mapping)
    np.random.seed(42)
    jitter_strength = 0.15 

    df_sex['Jitter_X'] = df_sex['Sex_Numeric'] + np.random.normal(0, jitter_strength, size=len(df_sex))
    df_sex['Jitter_Y'] = np.random.normal(0, jitter_strength, size=len(df_sex))

    plt.figure(figsize=(10, 8), facecolor='#f5f5f5')

    sns.scatterplot(
        x='Jitter_X', 
        y='Jitter_Y', 
        hue='pheno.Sex', 
        data=df_sex, 
        s=80, 
        alpha=0.7, 
        palette={'Male': '#1f77b4', 'Female': '#d62728'}, 
        edgecolor='black',
        linewidth=0.5
    )

    plt.title('Distribution of Sex in Cases', fontsize=16, fontweight='bold', pad=20)
    plt.xlabel('Sex Category (1=Male, 2=Female)', fontsize=12)
    plt.ylabel('Random Dispersion', fontsize=12)

    plt.axvline(1, color='grey', linestyle='--', alpha=0.3)
    plt.axvline(2, color='grey', linestyle='--', alpha=0.3)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.box(on=True)
    plt.tight_layout()
    plt.savefig('pca_sex.svg',format='svg', bbox_inches='tight')
    plt.show()
    return


@app.cell
def _(Ellipse, chi2, np, plt, scores_pd):
    def calculate_mahalanobis_genomic(x, data):
        covariance_matrix_g = np.cov(data, rowvar=False)
        covariance_matrix_inv_g = np.linalg.inv(covariance_matrix_g)
        mean_dist_g = x - np.mean(data, axis=0)
        return np.sqrt(np.dot(np.dot(mean_dist_g, covariance_matrix_inv_g), mean_dist_g.T))

    # 2. Calculation
    # scores_pd is already defined in memory; adding a new column to it
    pca_coords_all = scores_pd[['PC1', 'PC2']].values

    # Calculate distances
    scores_pd['mahalanobis_dist_all'] = [calculate_mahalanobis_genomic(x, pca_coords_all) for x in pca_coords_all]

    # 3. Threshold Value (p < 0.01)
    alpha_val = 0.01
    df_deg = 2
    critical_value_g = chi2.ppf((1 - alpha_val), df_deg)
    threshold_dist_g = np.sqrt(critical_value_g)

    scores_pd['is_outlier_genomic'] = scores_pd['mahalanobis_dist_all'] > threshold_dist_g

    # Print statistics
    print(f"Total individuals analyzed: {len(scores_pd)}")
    print(f"Number of outliers detected: {scores_pd['is_outlier_genomic'].sum()}")

    # 4. Visualization
    plt.figure(figsize=(10, 8))

    # Color map
    plot_colors = {"case": "red", "control": "blue"}

    # Renamed loop variables: 'status_label' and 'group_subset'
    for status_label in ["case", "control"]:
        # Normal data (Inlier)
        group_subset = scores_pd[(scores_pd['case_control_status'] == status_label) & (~scores_pd['is_outlier_genomic'])]
        plt.scatter(group_subset['PC1'], group_subset['PC2'], 
                    c=plot_colors[status_label], 
                    label=f'{status_label.capitalize()} (Inlier)', 
                    alpha=0.6, s=30)

        # Outlier data (Outlier) - used 'outlier_subset'
        outlier_subset = scores_pd[(scores_pd['case_control_status'] == status_label) & (scores_pd['is_outlier_genomic'])]
        if not outlier_subset.empty:
            plt.scatter(outlier_subset['PC1'], outlier_subset['PC2'], 
                        c=plot_colors[status_label], 
                        label=f'{status_label.capitalize()} (Outlier)', 
                        marker='x', s=100, linewidth=2)

    # Confidence Ellipse Drawing
    def draw_ellipse_genomic(x, y, ax, **kwargs):
        cov = np.cov(x, y)
        lambda_, v = np.linalg.eig(cov)
        lambda_ = np.sqrt(lambda_)
        scale_factor = np.sqrt(chi2.ppf(1 - alpha_val, 2))

        ell = Ellipse(xy=(np.mean(x), np.mean(y)),
                      width=lambda_[0] * 2 * scale_factor,
                      height=lambda_[1] * 2 * scale_factor,
                      angle=np.rad2deg(np.arccos(v[0, 0])), **kwargs)
        ax.add_patch(ell)

    current_ax = plt.gca()
    draw_ellipse_genomic(scores_pd['PC1'], scores_pd['PC2'], current_ax, 
                         edgecolor='black', linestyle='--', 
                         facecolor='none', linewidth=1.5, 
                         label='99% Confidence Ellipse')

    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.legend(loc='best')
    plt.title("Genomic PCA: Outlier Detection (Unique Vars)")
    plt.grid(True, linestyle='--', alpha=0.3)
    file_name_outlier = 'pca_outlier_mahalanobis.svg'
    plt.savefig(file_name_outlier, format='svg', bbox_inches='tight')
    plt.show()

    final_outliers_df = scores_pd[scores_pd['is_outlier_genomic']]
    return


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
    file_name_threshold = 'pca_outlier_threshold.svg'
    plt.savefig(file_name_threshold, format='svg', bbox_inches='tight')
    plt.show()
    return (outliers,)


@app.cell
def _(hl, result):
    # List of outliers and samples with incorrect sex
    outliers_and_wrong_sex = ['M_NG3215-1', 'NG3253-1', 'M_36', 'M_NG2605-2','M_277']

    # Convert the list to a Hail literal set
    remove_set = hl.literal(set(outliers_and_wrong_sex))

    cleaned_mt = result.filter_cols(
        ~remove_set.contains(result.s_corrected)
    )

    print("Number of columns before removal:", result.count_cols())
    print("Number of columns after removal:", cleaned_mt.count_cols())
    return (cleaned_mt,)


@app.cell
def _(outliers):
    # 1. Outlier örneklerin tam listesini konsola yazdırmak için:
    print("Tespit edilen outlier örnekler:")
    print(outliers)

    # 2. Sadece örnek isimlerini (Sample ID) listelemek için:
    # Not: Hail'den Pandas'a dönüşümde örnek isimleri genellikle 's' sütununda veya index'te tutulur.
    # Eğer isimler index ise:
    print("\nOutlier ID Listesi:")
    print(outliers.index.tolist())
    return


@app.cell
def _(cleaned_mt):
    cleaned_mt.write("kawasaki_filtered_outliers.mt", overwrite=True)
    return


if __name__ == "__main__":
    app.run()
