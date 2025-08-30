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
    return PCA, StandardScaler, hl, np, pd, plt, sns


@app.cell
def _(hl):
    #resultu kaydettik ya direkt açmak için bunu kullan artık 
    result = hl.read_matrix_table('kawasaki_filtered.mt')
    # ID'leri meta ile eşleştirecek şekilde düzelt (-DNA ekini kaldır)
    result = result.key_cols_by(s_corrected=result.s.replace('-DNA', ''))
    return (result,)


@app.cell
def _(result):
    result.describe()
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
def _(hl, meta_data_df, pd, result):
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
    final_mt_with_pheno = result.annotate_cols(
    pheno = pheno_table[result.s_corrected] # Düzeltilmiş ID'leri kullan
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
def _(PCA, StandardScaler, final_mt_with_case_control, np, pd, plt, sns):
    # Tüm PCA kodunu tek bir marimo hücresine koy
    # (Gerekiyorsa, `pd`, `np`, `final_mt_with_case_control` gibi değişkenleri de aynı hücrede tanımla)

    print("--- Her Değişken İçin Ayrı Ayrı PCA Analizi ---")

    # Case örneklerini filtrele (sadece case'lerde PCA yap)
    case_samples = final_mt_with_case_control.filter_cols(
        final_mt_with_case_control.case_control_status == 'case'
    )

    case_count = case_samples.count_cols()
    print(f"Case örnek sayısı: {case_count}")
    if case_count == 0:
        print("Hiç case örneği bulunamadı!")
        # return None

    # ... (case_df'yi oluşturan kod burada devam etmeli)
    case_pheno = case_samples.pheno.collect()
    case_df = pd.DataFrame(case_pheno)
    print(f"Case meta data shape: {case_df.shape}")
    print(f"Meta data sütunları: {case_df.columns.tolist()}")

    if case_df is None or case_df.empty:
        print("Case data bulunamadı!")

    # Sayısal sütunları seç (PCA için)
    numeric_columns = case_df.select_dtypes(include=[np.number]).columns.tolist()
    print(f"Sayısal sütunlar: {numeric_columns}")

    pca_results = {}
    for key in numeric_columns:
        print(f"\n--- {key} Değişkeni İçin PCA ---")
        # Veriyi hazırla
        data = case_df[key].dropna()
        if len(data) < 2:
            print(f" {key} için yeterli veri yok (en az 2 değer gerekli)")
            continue
        if data.nunique() < 2:
            print(f" {key} için varyasyon yok (tüm değerler aynı)")
            continue
        # Veriyi 2D array'e çevir (PCA için)
        X = data.values.reshape(-1, 1)
        # Standardize et
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
        # PCA uygula
        pca = PCA(n_components=1)
        X_pca = pca.fit_transform(X_scaled)
        # Sonuçları sakla
        pca_results[key] = {
            'original_data': data,
            'scaled_data': X_scaled,
            'pca_data': X_pca,
            'explained_variance_ratio': pca.explained_variance_ratio_[0],
            'components': pca.components_[0],
            'mean': scaler.mean_[0],
            'scale': scaler.scale_[0]
        }
        print(f" Açıklanan varyans oranı: {pca.explained_variance_ratio_[0]:.4f}")
        print(f" PCA bileşeni: {pca.components_[0]}")
        print(f" Veri sayısı: {len(data)}")
        print(f" Ortalama: {scaler.mean_[0]:.2f}")
        print(f" Standart sapma: {scaler.scale_[0]:.2f}")
    # --- Belirlenen Değişkenler için PCA Analizi ---
    print("--- Genel Veride (Hasta ve Kontrol) PCA Analizi Başlıyor ---")

    # Analiz edilecek 0/1 değişkenlerini seç
    pca_columns = [
        'CAA_status',
        'Sequelae_Status',
        'Family_History_Status',
        'Sex',
        'Consanguineous_marriage_status',
        'Degree_of_CM',
        'Family_history_of_CHD_status',
        'KD_in_siblings_status'
    ]

    # Tüm örnekleri (hem 'case' hem de 'control') içeren bir pandas DataFrame oluştur
    all_samples_pheno = final_mt_with_case_control.cols().collect()
    all_samples_df = pd.DataFrame(all_samples_pheno)
    pca_data = all_samples_df[['pheno'] + ['case_control_status']].dropna()
    pca_df = pd.DataFrame(pca_data['pheno'].tolist())
    pca_df['case_control_status'] = pca_data['case_control_status']

    # 'Sex' sütununu 'Male' -> 0, 'Female' -> 1 olarak dönüştür
    pca_df['Sex'] = pca_df['Sex'].map({'Male': 0, 'Female': 1})

    # Tüm değişkenlerin olduğu DataFrame'i temizle
    pca_df_all = pca_df.dropna(subset=pca_columns)
    X_all = pca_df_all[pca_columns].astype(float)
    y_all_status = pca_df_all['case_control_status']
    y_all_sex = pca_df_all['Sex']

    # Veriyi standartlaştır
    scaler_all = StandardScaler()
    X_scaled_all = scaler_all.fit_transform(X_all)

    # 2 ana bileşen ile PCA'i uygula
    pca_all = PCA(n_components=2)
    X_pca_all = pca_all.fit_transform(X_scaled_all)

    # PCA sonuçlarını bir DataFrame'e koy
    pca_results_df_all = pd.DataFrame(data=X_pca_all, columns=['PC1', 'PC2'])
    pca_results_df_all['status'] = y_all_status.reset_index(drop=True)
    pca_results_df_all['sex'] = y_all_sex.reset_index(drop=True)

    # 'status' sütunundaki değerleri daha anlaşılır etiketlere dönüştür
    pca_results_df_all['status_label'] = pca_results_df_all['status'].map({'case': 'Hasta', 'control': 'Kontrol'})
    # Diğer 0/1 değerlerini de etiketlere dönüştür
    pca_results_df_all['family_history_label'] = pca_df_all['Family_History_Status'].map({0.0: 'Yok', 1.0: 'Var'}).reset_index(drop=True)
    pca_results_df_all['caa_status_label'] = pca_df_all['CAA_status'].map({0.0: 'Yok', 1.0: 'Var'}).reset_index(drop=True)
    pca_results_df_all['sequelae_status_label'] = pca_df_all['Sequelae_Status'].map({0.0: 'Yok', 1.0: 'Var'}).reset_index(drop=True)
    pca_results_df_all['degree_of_cm_label'] = pca_df_all['Degree_of_CM'].map({0.0: 'Yok', 1.0: 'Var'}).reset_index(drop=True)
    pca_results_df_all['kd_in_siblings_label'] = pca_df_all['KD_in_siblings_status'].map({0.0: 'Yok', 1.0: 'Var'}).reset_index(drop=True)
    pca_results_df_all['consanguineous_marriage_label'] = pca_df_all['Consanguineous_marriage_status'].map({0.0: 'Yok', 1.0: 'Var'}).reset_index(drop=True)
    pca_results_df_all['sex_label'] = pca_df_all['Sex'].map({0.0: 'Erkek', 1.0: 'Kadın'}).reset_index(drop=True)

    # Açıklanan varyans oranını yazdır
    print("\n--- PCA Sonuçları ---")
    print(f"PC1'in açıkladığı varyans oranı: {pca_all.explained_variance_ratio_[0]:.2f}")
    print(f"PC2'nin açıkladığı varyans oranı: {pca_all.explained_variance_ratio_[1]:.2f}")
    print(f"Toplam açıklanan varyans: {pca_all.explained_variance_ratio_.sum():.2f}")
    print(f"Analizdeki toplam örnek sayısı: {len(pca_results_df_all)}")

    # Her bir değişken için ayrı bir kümelenme grafiği çiz
    plotting_columns = {
        'status_label': 'Hastalık Durumu',
        'sex_label': 'Cinsiyet',
        'family_history_label': 'Aile Geçmişi',
        'caa_status_label': 'CAA Durumu',
        'sequelae_status_label': 'Sequelae Durumu',
        'degree_of_cm_label': 'Akraba Evliliği Derecesi',
        'kd_in_siblings_label': 'Kardeşlerde Kawasaki',
        'consanguineous_marriage_label': 'Akraba Evliliği'
    }

    for plot_col, title in plotting_columns.items():
        print(f"\n--- {title} Değişkenine Göre Kümelenme ---")
        plt.figure(figsize=(10, 8))
        sns.scatterplot(
            x='PC1', 
            y='PC2', 
            hue=plot_col, 
            data=pca_results_df_all, 
            palette='deep', 
            s=100
        )
        plt.title(f"Kawasaki Hastalığı: {title} Değişkenine Göre PCA Analizi")
        plt.xlabel(f"Principal Component 1 ({pca_all.explained_variance_ratio_[0]*100:.2f}%)")
        plt.ylabel(f"Principal Component 2 ({pca_all.explained_variance_ratio_[1]*100:.2f}%)")
        plt.grid(True)
        plt.show()

    return


if __name__ == "__main__":
    app.run()
