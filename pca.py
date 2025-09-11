import marimo

__generated_with = "0.14.9"
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
    return PCA, hl, np, pd, plt, sns


@app.cell
def _(hl):
    # We saved the result already, so use this to open it directly
    result = hl.read_matrix_table('kawasaki_filtered.mt')

    # Adjust to match IDs with metadata (remove the '-DNA' suffix)
    result = result.key_cols_by(s_corrected=result.s.replace('-DNA', ''))

    return (result,)


@app.cell
def _(result):
    result.describe()
    return


@app.cell
def _(pd):
    # Read the TSV file
    meta_data_df = pd.read_csv("meta.tsv", sep="\t", header=0)

    print("--- Original Column Names ---")
    print(meta_data_df.columns.tolist())

    # Clean column names to make them more usable
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

    print("\n--- Cleaned Column Names (Final Version) ---")
    print(meta_data_df.columns.tolist())

    # Format the IDs in the 'Case_Number' column by adding the prefix 'M_'
    id_column_name = 'Case_Number'
    meta_data_df['formatted_sample_id'] = 'M_' + meta_data_df[id_column_name].astype(str)

    print("\n--- Metadata IDs formatted with 'M_' prefix (First 5 rows) ---")
    print(meta_data_df[['Case_Number', 'formatted_sample_id']].head())

    return (meta_data_df,)


@app.cell
def _(hl, meta_data_df, pd, result):

    print("--- Data Types and Values Detailed Check ---")
    print("\nDataFrame Shape:", meta_data_df.shape)

    print("\nColumn Names:")
    for col in meta_data_df.columns:
        print(f"- {col}") 

    print("\n--- Detailed Analysis of Each Column ---")
    for col in meta_data_df.columns:
        print(f"\n{col}:")
        print(f" Data type: {meta_data_df[col].dtype}") 
        print(f" Missing values: {meta_data_df[col].isna().sum()}") 
        print(f" Number of unique values: {meta_data_df[col].nunique()}") 
        print(f" First 10 value counts: {meta_data_df[col].value_counts().head(10).to_dict()}") 

        # Additional statistics for numeric columns
        if pd.api.types.is_numeric_dtype(meta_data_df[col]): 
            print(f" Min: {meta_data_df[col].min()}") 
            print(f" Max: {meta_data_df[col].max()}") 
            print(f" Mean: {meta_data_df[col].mean():.2f}") 

    # Only modify missing values in the 'Class' column, leave others unchanged
    print("--- Adjusting Only the Class Column ---")

    if 'Class' in meta_data_df.columns:
        print("Modifying Class column...")
        print(f"Original values: {meta_data_df['Class'].value_counts().head()}")

        # Replace missing values with 'nan'
        meta_data_df['Class'] = meta_data_df['Class'].fillna('nan')

        print(f"Values after modification: {meta_data_df['Class'].value_counts().head()}")
        print("Other columns remain unchanged.")

    # Convert the prepared Pandas DataFrame into a Hail Table
    pheno_table = hl.Table.from_pandas(meta_data_df, key='formatted_sample_id')

    print("\n--- Hail Phenotype Table Created and Types Checked (First 5 Rows) ---")
    pheno_table.show(192)

    # Annotate the MatrixTable with phenotype data
    final_mt_with_pheno = result.annotate_cols(
        pheno = pheno_table[result.s_corrected] 
    )

    # Finally, assign 'case' and 'control' labels
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
def _(PCA, final_mt_with_case_control, hl, np, pd, plt, sns):
    # Compute PCA (first 2 components)
    eigenvalues, scores, loadings = hl.hwe_normalized_pca(final_mt_with_case_control.GT, k=2)

    # Add case/control information to the scores Table
    scores = scores.annotate(case_control_status=final_mt_with_case_control.cols()[scores.s_corrected].case_control_status)

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
    case_mt = final_mt_with_case_control.filter_cols(final_mt_with_case_control.case_control_status == "case")

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

    return (outliers,)


@app.cell
def _(final_mt_with_case_control, hl, outliers):
    # 1. Convert the IDs of outliers (s_corrected) into a list
    outlier_ids = outliers['s_corrected'].tolist()

    # 2. Convert this list into a set that Hail can use
    outlier_set = hl.literal(outlier_ids)

    # 3. Filter out the outliers from the main MatrixTable
    clean_data = final_mt_with_case_control.filter_cols(
        ~outlier_set.contains(final_mt_with_case_control.s_corrected)
    )

    # 4. Check the size of the cleaned dataset
    print("\n--- Dataset Size Check ---")
    print(f"Original dataset size (number of samples): {final_mt_with_case_control.count_cols()}")
    print(f"Cleaned dataset size (number of samples): {clean_data.count_cols()}")

    # See the IDs of the individuals filtered as outliers
    print("\n--- Outlier IDs ---")
    print(outlier_ids)

    return


if __name__ == "__main__":
    app.run()
