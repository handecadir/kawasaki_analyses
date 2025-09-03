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
    return PCA, hl, mo, np, pd, plt, sns


@app.cell
def _(hl):
    nadir_filtered_mis = hl.read_matrix_table('missense.mt')
    lof_filtered= hl.read_matrix_table('lof_variants.mt')
    # 2. Entry şemaları aynı mı kontrol et
    print("Missense entry fields:", nadir_filtered_mis.entry)
    print("LoF entry fields:", lof_filtered.entry)
    return lof_filtered, nadir_filtered_mis


@app.cell
def _(lof_filtered, nadir_filtered_mis):
    # MT’leri birleştir
    combined_mt = nadir_filtered_mis.union_rows(lof_filtered)  # MT methodu olarak çağır
    # Mevcut ID'leri string olarak al
    # ID'leri meta ile eşleştirecek şekilde düzelt (-DNA ekini kaldır)
    combined_mt = combined_mt.key_cols_by(s_corrected=combined_mt.s.replace('-DNA', ''))


    # Kontrol
    print("\n--- Birleştirilmiş MatrixTable ---")
    print("Varyant sayısı, örnek sayısı:", combined_mt.count())
    combined_mt.describe()

    return (combined_mt,)


@app.cell
def _(combined_mt):
    combined_mt.describe()
    return


@app.cell
def _(pd):
    # TSV dosyasını oku
    meta_data_df = pd.read_csv("meta.tsv", sep="\t", header=0)

    print("--- Orijinal Sütun İsimleri ---")
    print(meta_data_df.columns.tolist())

    # 2. Sütun isimlerini temizle ve daha kullanışlı hale getir
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

    print("\n--- Temizlenmiş Sütun İsimleri (Son Hali) ---")
    print(meta_data_df.columns.tolist())

    # 3. 'Case_Number' sütunundaki ID'lerin önüne 'M_' ekleyerek formatla
    id_column_name = 'Case_Number'
    meta_data_df['formatted_sample_id'] = 'M_' + meta_data_df[id_column_name].astype(str)

    print("\n--- Meta Veri ID'leri 'M_' formatına çevrildi (İlk 5 satır) ---")
    print(meta_data_df[['Case_Number', 'formatted_sample_id']].head())

    return (meta_data_df,)


@app.cell
def _(combined_mt, hl, meta_data_df, pd):
    # Doğru kod ✅
    print("--- Veri Tipleri ve Değerler Detaylı Kontrol ---")
    print("\nDataFrame Shape:", meta_data_df.shape)
    print("\nSütun İsimleri:")
    for col in meta_data_df.columns:
        print(f"- {col}") # ✅ Doğru girinti
    print("\n--- Her Sütunun Detaylı Analizi ---")
    for col in meta_data_df.columns:
        print(f"\n{col}:") # ✅ Doğru girinti
        print(f" Veri tipi: {meta_data_df[col].dtype}") # ✅ Doğru girinti
        print(f" Boş değer sayısı: {meta_data_df[col].isna().sum()}") # ✅ Doğru girinti
        print(f" Benzersiz değer sayısı: {meta_data_df[col].nunique()}") # ✅ Doğru girinti
        print(f" İlk 10 değer: {meta_data_df[col].value_counts().head(10).to_dict()}") # ✅ Doğru girinti
        # Sayısal sütunlar için ek istatistikler
        if pd.api.types.is_numeric_dtype(meta_data_df[col]): # ✅ Doğru girinti
            print(f" Min: {meta_data_df[col].min()}") # ✅ Doğru girinti
            print(f" Max: {meta_data_df[col].max()}") # ✅ Doğru girinti
            print(f" Mean: {meta_data_df[col].mean():.2f}") # ✅ Doğru girinti
    # 4. Sadece Class sütunundaki boş değerleri 'nan' yap, diğerlerine dokunma
    print("--- Sadece Class Sütunu Düzenleniyor ---")
    # Class sütunundaki boş değerleri 'nan' yap
    # Doğru kod ✅
    if 'Class' in meta_data_df.columns:
        print("Class sütunu düzenleniyor...")
        print(f"Orijinal değerler: {meta_data_df['Class'].value_counts().head()}")
        # Boş değerleri 'nan' yap
        meta_data_df['Class'] = meta_data_df['Class'].fillna('nan')
        print(f"Düzenleme sonrası değerler: {meta_data_df['Class'].value_counts().head()}")
        print("Diğer sütunlar olduğu gibi bırakıldı.")
    # 5. Hazırlanan Pandas DataFrame'i Hail Table'a dönüştür ('types' argümanı olmadan)
    pheno_table = hl.Table.from_pandas(meta_data_df, key='formatted_sample_id')

    print("\n--- Hail Fenotip Tablosu Oluşturuldu ve Tipleri Kontrol Edildi (İlk 5 satır) ---")
    pheno_table.show(192)

    # BURASI EN KRİTİK KISIM!
    final_mt_with_pheno = combined_mt.annotate_cols(
    pheno = pheno_table[combined_mt.s_corrected] # Düzeltilmiş ID'leri kullan
    )

    # 7. Son olarak, 'case' ve 'control' etiketlerini ata
    # Artık `final_mt_with_pheno` tablosundaki `pheno` alanının dolu olup olmadığını kontrol edebiliriz
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
    # Sadece ilk 100 örneği ve kolonları gösterir
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

    # Sadece bu kolonları al (gerekirse int'e çevir)
    pca_data = df[features].astype(int)

    # PCA
    pca = PCA(n_components=2)
    principal_components = pca.fit_transform(pca_data)

    # PCA sonuç DataFrame
    pca_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'])

    # Benzersiz noktaları kontrol et
    unique_points = pca_df.drop_duplicates()
    print(f"\nNumber of unique PCA points: {len(unique_points)}")
    print(f"Total number of PCA points: {len(pca_df)}")

    # Jitter ekle
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

    # Her feature için ayrı scatter
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
    return


@app.cell
def _(mo):
    mo.md(
        r"""
    --- Aykırı Değerlerin Kimlikleri ---
    ['M_276', 'M_277', 'M_NG3215-1', 'NG3253-1']
    """
    )
    return


if __name__ == "__main__":
    app.run()
