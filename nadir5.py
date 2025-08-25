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
    return hl, np, os, pd, plt, sm, sns


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

    # Kontrol
    print("\n--- Birleştirilmiş MatrixTable ---")
    print("Varyant sayısı, örnek sayısı:", combined_mt.count())
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
def _():



    return


@app.cell
def _(hl, meta_data_df, nadir_filtered, pd):
    # 4. Boş değerleri doldur ve VERİ TİPLERİNİ PANDAS'TA KESİN OLARAK BELİRT
    boolean_columns = [
        'CAA_status', 'Diagnostic_Age_Status', 'Sequelae_Status', 'Family_History_Status',
        'Consanguineous_marriage_status', 'Family_history_of_CHD_status',
        'KD_in_siblings_status'
    ]
    for col in boolean_columns:
        if col in meta_data_df.columns:
            meta_data_df[col] = meta_data_df[col].astype(str).str.strip().str.lower().map({'true': 1, 'false': 0}).fillna(0).astype(int)

    numerical_columns_int = ['Record_No', 'Case_Number', 'Degree_of_CM']
    for col in numerical_columns_int:
        if col in meta_data_df.columns:
            meta_data_df[col] = pd.to_numeric(meta_data_df[col], errors='coerce').fillna(0).astype(int)

    numerical_columns_float = ['Age_at_diagnosis_years']
    for col in numerical_columns_float:
        if col in meta_data_df.columns:
            meta_data_df[col] = pd.to_numeric(meta_data_df[col], errors='coerce').fillna(0.0).astype(float)

    string_columns = ['Class', 'Sex']
    for col in string_columns:
        if col in meta_data_df.columns:
            meta_data_df[col] = meta_data_df[col].astype(str).fillna('')

    print("\n--- Pandas DataFrame Sütun Tipleri (Dönüşüm Sonrası) ---")
    print(meta_data_df.dtypes)


    # 5. Hazırlanan Pandas DataFrame'i Hail Table'a dönüştür ('types' argümanı olmadan)
    pheno_table = hl.Table.from_pandas(meta_data_df, key='formatted_sample_id')

    print("\n--- Hail Fenotip Tablosu Oluşturuldu ve Tipleri Kontrol Edildi (İlk 5 satır) ---")
    pheno_table.show(5)

    # 6. Fenotip verilerini MatrixTable'a ekle
    final_mt_with_pheno = nadir_filtered.annotate_cols(pheno = pheno_table[nadir_filtered.s])

    print("\n--- MatrixTable'a Eklenmiş Fenotip Bilgileri (İlk 5 Örnek) ---")
    final_mt_with_pheno.cols().show(5)


    # 7. Eksik Klinik Verisi Olan Örnekleri 'Kontrol' Olarak Etiketleme
    final_mt_with_case_control = final_mt_with_pheno.annotate_cols(
        pheno = final_mt_with_pheno.pheno.annotate(
            case_control_status = hl.if_else(
                hl.is_defined(final_mt_with_pheno.pheno.Case_Number),
                'hasta',
                'kontrol'
            )
        )
    )


    return (final_mt_with_case_control,)


@app.cell
def _(final_mt_with_case_control, hl, np, os, pd):
    print("\n--- MatrixTable'a Eklenmiş Kontrol/Hasta Durumu (İlk 10 Örnek) ---")



    # 'final_mt_with_case_control.cols()' metodunun döndürdüğü bir 'TableExpression' var.

    # Bu TableExpression'ın içinde 'pheno' adında bir alan var.

    # 'select' metodunun içine doğrudan, bu 'TableExpression'ın içindeki 'pheno' alanını kullanarak erişmeliyiz.



    # Bunu yapmanın en kesin yolu, select içine o anki objenin alanını yazmaktır.

    # Yani:

    temp_cols_table = final_mt_with_case_control.cols()



    # Şimdi bu geçici 'temp_cols_table' objesinin içindeki 'pheno'ya erişebiliriz:

    temp_cols_table.select(

        # 'pheno' zaten 'temp_cols_table'ın içinde olduğu için,

        # doğrudan 'pheno.case_control_status' diyebiliriz.

        # Burada 'temp_cols_table.pheno.case_control_status' yazmamız gerekiyor.

        case_control_status = temp_cols_table.pheno.case_control_status,

        Sex = temp_cols_table.pheno.Sex,

        CAA_status = temp_cols_table.pheno.CAA_status

    ).show(200)



    # --- 1. ADIM: Harmonizome Gen Listesini Kaydet ve Yükle ---

    harmonizome_genes = [

        "MLIP", "PTCD1", "CADM1", "EBF2", "CACNA2D4", "TECRL", "PRKAB1", "TANC2", "SNRPA", "TFCP2",

        "ZNF717", "CLEC4A", "RBP2", "ZSCAN25", "ROR2", "HCG23", "LOC100130301", "LOC100287329",

        "PTPN4", "PARK2", "C6ORF10", "IL1RAP", "HS6ST3", "ATP5J2", "CLEC3A", "PROKR2", "RELN",

        "PREX1", "FAM184A", "BMPER", "HLA-DOB", "NUMBL", "RAB4B", "BUD31", "CASK", "LOC100507291",

        "NAALADL2", "TMEM233", "SLC6A6", "TDRD1", "LOC100287225", "ZNF536", "ZNF655", "ZFPM2",

        "MAD1L1", "VIL1", "BLK", "LDLRAD3", "SMC5", "SCGN", "LOC285692", "RIMS2", "BTBD1",

        "MED30", "PDGFD", "KCNN2", "CENPU", "ZFHX3", "CADPS", "NCOA5", "HEG1", "PXDNL", "MYRIP",

        "DNM3", "ITPKC", "FCGR2A", "DAB1", "HEATR5B", "HLA-DRA", "IRF8", "ARPP21", "RAPGEF1",

        "KCNK3", "HLA-DQB2", "SHC4", "CD40", "SLC35F5", "FHIT", "KLC1", "LRRC16A", "DEC1",

        "HRH2", "GABRR1", "ERBB4", "RDH16", "CCDC186", "HLA-DQA2", "EGLN2", "ERAP1", "TTC1",

        "SLC39A13", "DLK1", "MKLN1", "FAM167A", "COL4A2", "ALDH1A2", "PELI1", "TLL1", "GCOM1",

        "MTCH2", "FLT1", "DDX60L", "PRICKLE1", "PPM1L", "CLCN3", "ROBO1", "FSTL4", "EPB41L4B",

        "SGCD", "ERC2", "HNRNPA3", "CSMD1", "NCK2"

    ]



    go_genes_file = 'go_target_genes.txt'



    if not os.path.exists(go_genes_file) or os.stat(go_genes_file).st_size == 0:

        print(f"\n'{go_genes_file}' dosyası bulunamadı veya boş. Genler dosyaya kaydediliyor.")

        with open(go_genes_file, 'w') as f:

            for gene in harmonizome_genes:

                f.write(gene + '\n')

    else:

        print(f"\n'{go_genes_file}' dosyası zaten mevcut. Dosya kaydedilmedi, mevcut dosya kullanılacak.")



    try:

        with open(go_genes_file, 'r') as f:

            target_genes_list = [line.strip() for line in f if line.strip()]

        hl_target_genes = hl.set(target_genes_list)



        print(f"'{go_genes_file}' dosyasından {len(target_genes_list)} hedef gen yüklendi. İlk 10 gen: {target_genes_list[:10]}...")



    except FileNotFoundError:

        print(f"\nUYARI: '{go_genes_file}' dosyası bulunamadı. Lütfen GO gen listesi dosyasını kontrol edin.")

        print("Bu analiz için hedef gen listesi olmadan ilerleyemeyiz.")

        exit()



    # --- 2. ADIM: MatrixTable'ı Hedef Genlere Göre Filtrele ---

    mt_filtered_genes = final_mt_with_case_control.filter_rows(

        hl.is_defined(final_mt_with_case_control.gene_symbol) & \

        hl_target_genes.contains(final_mt_with_case_control.gene_symbol) & \

        hl.is_defined(final_mt_with_case_control.is_canonical) & \

        (final_mt_with_case_control.is_canonical == 1)

    )



    print(f"\nKawasaki hedef genlere göre filtrelenen varyant sayısı: {mt_filtered_genes.count_rows()}")



    # --- 3. ADIM: Nihai Filtrelenmiş MatrixTable ---

    filtered_mt_pato_final = mt_filtered_genes



    print(f"\nNihai filtrelenmiş varyant sayısı (Kawasaki Gen Seti ve Patojenik Filtreler sonrası): {filtered_mt_pato_final.count_rows()}")



    # Eğer hiç varyant bulunamazsa, sonraki adımları atla

    if filtered_mt_pato_final.count_rows() == 0:

        print("\nFiltrelenmiş varyant bulunamadı, varyant yükü hesaplanamaz veya regresyon yapılamaz.")

    else:

        # --- 4. ADIM: Her Birey İçin Varyant Yükünü Hesapla ---

        mt_with_variant_load = filtered_mt_pato_final.annotate_cols(

            variant_load = hl.agg.count_where(

                hl.is_defined(filtered_mt_pato_final.GT) & \

                ~filtered_mt_pato_final.GT.is_hom_ref()

            )

        )



        temp_cols_table = mt_with_variant_load.cols()



        sample_variant_counts_ht = temp_cols_table.select(

            variant_load = temp_cols_table.variant_load

        )



        sample_variant_counts_df = sample_variant_counts_ht.to_pandas()



        print("\nHer birey için varyant yükü hesaplandı (İlk 5):")

        print(sample_variant_counts_df.head())



        # --- 5. ADIM: Varyant Yükünü Fenotip Bilgileriyle Birleştir (YİNE DÜZELTİLDİ! BU SEFER KESİN! ✨) ---

        # MatrixTable'ın sütunlarını bir Table'a dönüştürürken, onun kendi alanlarını kullanmalıyız.

        pheno_table_from_mt_cols = final_mt_with_case_control.cols() # Yeni bir Table objesi oluşturduk



        pheno_data_for_merge_ht = pheno_table_from_mt_cols.select(

            pheno_table_from_mt_cols.pheno.case_control_status, # <-- Düzeltme burada!

            pheno_table_from_mt_cols.pheno.Sex                  # <-- Düzeltme burada!

        )



        pheno_data_for_merge_ht = pheno_data_for_merge_ht.key_by('s') # Zaten anahtar 's' olacak



        final_df = pd.merge(sample_variant_counts_df, pheno_data_for_merge_ht.to_pandas(), on='s', how='left')



        final_df['group_numeric'] = final_df['case_control_status'].apply(

            lambda x: 1 if x == 'hasta' else (0 if x == 'kontrol' else np.nan)

        )



        print("\nVaryant yükü ve grup bilgileri birleştirildi (İlk 5):")

        print(final_df.head())
    return (final_df,)


@app.cell
def _(final_df, np, pd, plt, sm, sns):
    # ... (önceki kodlarının hepsi burada bitiyor, yani 'final_df' burada oluşturulmuş oluyor) ...

    # --- Lojistik Regresyon Öncesi Veri Tipi Kontrol ve Nihai Dönüşüm (ÖNEMLİ!) ---
    # Bu blok, lojistik regresyon adımından HEMEN ÖNCE gelmeli.
    # Amacımız: df_for_regression'ı burada oluşturmak ve sütun tiplerini statsmodels için uygun hale getirmek.

    print("\n--- Lojistik Regresyon Öncesi Veri Tipi Kontrolleri ---")
    print("final_df['variant_load'] tipi:", final_df['variant_load'].dtype)
    print("final_df['group_numeric'] tipi:", final_df['group_numeric'].dtype)

    # df_for_regression'ı burada oluşturuyoruz.
    # final_df'deki NaN değerleri içeren satırları düşürüyoruz ki regresyon modelimiz sorunsuz çalışsın.
    df_for_regression = final_df.dropna(subset=['group_numeric', 'variant_load']).copy()

    # Şimdi df_for_regression içindeki sütunları kesin olarak float tipine dönüştürüyoruz.
    # Bu, statsmodels'ın "object" hatasını gidermede en güvenli yoldur.
    df_for_regression['variant_load'] = df_for_regression['variant_load'].astype(float)
    df_for_regression['group_numeric'] = df_for_regression['group_numeric'].astype(float)

    print("\n--- Dönüşüm Sonrası df_for_regression Sütun Tipleri ---")
    print("df_for_regression['variant_load'] tipi:", df_for_regression['variant_load'].dtype)
    print("df_for_regression['group_numeric'] tipi:", df_for_regression['group_numeric'].dtype)
    print("Dönüşüm sonrası df_for_regression'da kalan satır sayısı:", len(df_for_regression))

    # Ayrıca, regresyon öncesi varyant yükü ve grup dağılımını kontrol edelim.
    # Eğer 'group_numeric' sütununda sadece tek bir değer varsa (örn. sadece 0 veya sadece 1),
    # veya 'variant_load' sütunundaki tüm değerler aynıysa, lojistik regresyon yapılamaz.
    if df_for_regression['group_numeric'].nunique() < 2:
        print("\n❗ YETERLİ GRUP BİLGİSİ YOK: 'group_numeric' sütununda en az iki farklı grup (0 ve 1) bulunamadı. Lojistik regresyon yapılamaz.")
    elif df_for_regression['variant_load'].nunique() <= 1:
        print("\n❗ VARYANT YÜKÜNDE YETERLİ DEĞİŞKENLİK YOK: 'variant_load' sütunundaki tüm değerler aynı veya çok benzer. Lojistik regresyon yapılamaz.")
    else:
        # --- 6. ADIM: Lojistik Regresyon ve İstatistiksel Karşılaştırma ---
        # Artık 'df_for_regression' tanımlı ve veri tipleri doğru.
        X = df_for_regression[['variant_load']]
        y = df_for_regression['group_numeric']

        X = sm.add_constant(X) # Sabit terimi (intercept) modele ekliyoruz.

        try:
            model = sm.Logit(y, X)
            # Modelin yakınsama sorunları yaşaması ihtimaline karşı maxiter'ı artırıyoruz.
            result_logreg = model.fit(disp=False, maxiter=1000)

            print("\n--- Lojistik Regresyon Sonuçları ---")
            print(result_logreg.summary())

            # Olasılık Oranları (Odds Ratios) ve %95 Güven Aralıkları hesaplama
            odds_ratios = pd.DataFrame({
                'OR': result_logreg.params,
                'Lower CI': result_logreg.conf_int()[0],
                'Upper CI': result_logreg.conf_int()[1]
            })
            odds_ratios = odds_ratios.apply(lambda x: np.exp(x))
            print("\n--- Olasılık Oranları (Odds Ratios) ve %95 Güven Aralıkları ---")
            print(odds_ratios)

            # --- 7. ADIM: Görselleştirme ---
            plt.figure(figsize=(8, 6))
            # Kutulu grafik (boxplot) ile varyant yükünün gruplara göre dağılımını görselleştiriyoruz.
            # Burada df_for_regression kullanmak, analiz ettiğimiz temizlenmiş veri setini görselleştirmek için daha tutarlı.
            sns.boxplot(x='group_numeric', y='variant_load', data=df_for_regression)
            plt.title('Varyant Yükünün Gruplara Göre Dağılımı (Kawasaki Gen Seti)')
            plt.xlabel('Grup Durumu (0: Kontrol, 1: Hasta)')
            plt.ylabel('Varyant Yükü')
            plt.xticks(ticks=[0, 1], labels=['Kontrol', 'Hasta']) # Eksen etiketlerini daha açıklayıcı yapıyoruz.
            plt.grid(axis='y', linestyle='--', alpha=0.7)
            plt.show()

        except Exception as e:
            print(f"\nLojistik regresyon sırasında bir hata oluştu: {e}")
            print("Modelin yakınsama sorunları olabilir veya veri setinde sorunlar olabilir (örn. tüm varyant yükleri aynı).")

    return (df_for_regression,)


@app.cell
def _(df_for_regression, np):
    np.asarray(df_for_regression)
    return


@app.cell
def _(final_mt_with_case_control):
    final_mt_with_case_control.describe()
    return


if __name__ == "__main__":
    app.run()
