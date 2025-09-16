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
    nadir_filtered_mis = hl.read_matrix_table('missense.mt')
    # 2. Entry şemaları aynı mı kontrol et
    print("Missense entry fields:", nadir_filtered_mis.entry)
    # ID'leri meta ile eşleştirecek şekilde düzelt (-DNA ekini kaldır)
    nadir_filtered_mis  = nadir_filtered_mis .key_cols_by(s_corrected=nadir_filtered_mis .s.replace('-DNA', ''))
    nadir_filtered_mis.describe()
    return (nadir_filtered_mis,)


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
def _(hl, meta_data_df, nadir_filtered_mis, pd):
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
    final_mt_with_pheno = nadir_filtered_mis.annotate_cols(
    pheno = pheno_table[nadir_filtered_mis.s_corrected] # Düzeltilmiş ID'leri kullan
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
def _(mo):
    mo.md(
        r"""
    --- Aykırı Değerlerin ve yanlış cinsiyet eşleşmesinin Kimlikleri ---
    'M_276', 'M_277', 'M_NG3215-1', 'NG3253-1', M_36
    """
    )
    return


@app.cell
def _(final_mt_with_case_control, hl):
    # Çıkarılacak örneklerin listesi
    outliers_and_wrong_sex = ['M_276', 'M_277', 'M_NG3215-1', 'NG3253-1', 'M_36']

    # Listeyi Hail literal set'e dönüştürelim
    remove_set = hl.literal(set(outliers_and_wrong_sex))

    # final_mt_with_case_control'dan bu örnekleri çıkar
    cleaned_mt = final_mt_with_case_control.filter_cols(
        ~remove_set.contains(final_mt_with_case_control.s_corrected)
    )

    print("Önceki kolon sayısı:", final_mt_with_case_control.count_cols())
    print("Temizlenmiş kolon sayısı:", cleaned_mt.count_cols())

    return (cleaned_mt,)


@app.cell
def _(cleaned_mt):
    # Sadece ilk 100 örneği ve kolonları gösterir
    cleaned_mt.cols().show(187)
    return


@app.cell
def _(cleaned_mt):
    cleaned_mt.describe()
    return


@app.cell
def _(cleaned_mt, hl):
    # 1. Adım: Toplam vaka ve kontrol örnek sayısını bulma 🕵️‍♀️
    total_cases = cleaned_mt.aggregate_cols(hl.agg.count_where(cleaned_mt.case_control_status == 'case'))
    total_controls = cleaned_mt.aggregate_cols(hl.agg.count_where(cleaned_mt.case_control_status == 'control'))

    print(f"Toplam vaka sayısı: {total_cases}")
    print(f"Toplam kontrol sayısı: {total_controls}")


    # 2. Adım: Her gen için vaka ve kontrol varyant sayılarını sayma 📊
    # Bu işlem sonucunda oluşan n_cases_variants ve n_controls_variants Entry fields'a eklenir.
    gene_counts_mt = cleaned_mt.group_rows_by(cleaned_mt.gene_symbol).aggregate(
        n_cases_variants=hl.agg.count_where(cleaned_mt.case_control_status == 'case'),
        n_controls_variants=hl.agg.count_where(cleaned_mt.case_control_status == 'control')
    )

    # Entry fields'ı satırlara (rows) dönüştürerek bir Table oluşturuyoruz.
    # Bu Table, artık gen sembolü ile birlikte varyant sayılarını içerecek.
    gene_counts_table = gene_counts_mt.entries()

    gene_counts_table.show(5)


    # 3. Adım: Fisher's exact test ile gen burden analizi yapma 🧬
    # gene_counts_table objesi bir Table olduğu için direkt annotate metodunu kullanıyoruz.
    gene_results = gene_counts_table.annotate(
        fisher_result = hl.fisher_exact_test(
            hl.int32(gene_counts_table.n_cases_variants),
            hl.int32(total_cases - gene_counts_table.n_cases_variants),
            hl.int32(gene_counts_table.n_controls_variants),
            hl.int32(total_controls - gene_counts_table.n_controls_variants)
        )
    )

    # 4. Adım: Sonuçları sıralama ve gösterme 🏆
    gene_burden_results = gene_results.order_by(gene_results.fisher_result.p_value)

    print("\nEn anlamlı 10 gen burden analizi sonucu:")
    gene_burden_results.show(1000)

    return


if __name__ == "__main__":
    app.run()
