import marimo

__generated_with = "0.16.5"
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
    import requests
    import re
    return hl, mo, pd, re, requests


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
    'M_276', 'M_277', 'M_NG3215-1', 'NG3253-1', 'M_36', 'NG2605-1'
    """
    )
    return


@app.cell
def _(final_mt_with_case_control, hl):
    # Çıkarılacak örneklerin listesi
    outliers_and_wrong_sex = ['M_276', 'M_277', 'M_NG3215-1', 'NG3253-1', 'M_36','NG2605-1']

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
def _(re, requests):

    def get_pathway_gene_symbols(pathway_id):
        """
        KEGG API'den verilen pathway ID'sine ait gen sembollerini ve pathway adını çeker.
        """
        url = f"http://rest.kegg.jp/get/{pathway_id}"
        try:
            resp = requests.get(url, timeout=10)
            resp.raise_for_status()  # Hata kodlarını (4xx, 5xx) yakalamak için
        except requests.exceptions.RequestException as e:
            print(f"KEGG API'den '{pathway_id}' için veri alınırken hata oluştu: {e}")
            return None, None

        text = resp.text

        pathway_name = ""
        symbols = []
        gene_section = False

        for line in text.split("\n"):
            if line.startswith("NAME"):
                pathway_name = line[len("NAME"):].strip()
            elif line.startswith("GENE"):
                gene_section = True
                line = line[len("GENE"):].strip()
            elif gene_section and re.match(r"^\s+\d+", line):
                line = line.strip()
            else:
                if gene_section:
                    break
                else:
                    continue

            if gene_section and line:
                parts = line.split()
                if len(parts) >= 2:
                    # "GENE" satırından sonraki satırlar birden fazla gen içerebiliyor.
                    # Örneğin: "1234 ABC; 5678 XYZ;"
                    # Bu yüzden her birini ayıklamamız lazım.
                    gene_symbols_on_line = re.findall(r"\s(\S+?);", line)
                    symbols.extend(gene_symbols_on_line)

        # Yinelenen sembolleri kaldır
        symbols = list(set(symbols))
        return symbols, pathway_name

    def main():
        """
        Belirtilen pathway ID'leri için gen sembollerini çeker ve sonuçları yazdırır.
        """
        # Verilen listeyi düzenledim. H ve map ID'lerini "hsa" olarak değiştirdim
        # ve liste elemanlarını ayırdım.
        pathway_ids = ["hsa04370", "hsa04060", "hsa04630", "hsa04620", "hsa04621", "hsa04660", "hsa04657", "hsa04151", "hsa04010", "hsa04350"]

        print("✨ Güncellenen pathway'lerin gen sembolleri ve isimleri: ✨\n")

        for p_id in pathway_ids:
            symbols, name = get_pathway_gene_symbols(p_id)
            if symbols and name:
                print(f"➡️ **{name} ({p_id})**")
                print("   Gen sembolleri:", ", ".join(symbols))
                print("   Toplam gen sayısı:", len(symbols))
                print("-" * 30)
            elif name is None:
                print(f"🤔 {p_id} ID'si için KEGG'den veri alınamadı. Atlanıyor.")
                print("-" * 30)


    # Kodu çalıştırmak için
    if __name__ == "__main__":
        main()
    return


@app.cell
def _():
    vegf= ["PIK3R3", "PPP3R1", "PRKCA", "MAP2K2", "RAC1", "MAPK1", "RAC2", "MAPK13", "AKT2", "HRAS", "VEGFA", "NOS3", "PIK3R2", "HSPB1", "PIK3R1", "PLCG2", "SPHK2", "PTGS2", "PIK3CA", "SHC2", "P3R3URF", "PIK3R3", "CASP9", "MAPK11", "SPHK1", "MAP2K1", "PXN", "PLA2G4E", "PLA2G4F", "PTK2", "AKT1", "PPP3R2", "AKT3", "BAD", "PLCG1", "RAF1", "PRKCG", "MAPKAPK2", "CDC42", "NFATC2", "NRAS", "JMJD7", "PLA2G4B", "PRKCB", "PPP3CB", "PLA2G4D", "PLA2G4B", "MAPK14", "MAPK12", "PPP3CC", "SH2D2A", "PIK3CB", "PIK3CD", "MAPKAPK3", "PLA2G4C", "KDR", "RAC3", "SRC", "PPP3CA", "PLA2G4A", "MAPK3", "KRAS", "RNASE1", "NANOS3", "PXDN","RNASE3", "NCOA3"]

    Cytokine= ["CCL20", "IL1RL1", "IL22", "GHR", "CSF1R", "ACVRL1", "IL23R", "CCL1", "BMP10", "PF4", "CCL28", "CCR9", "GDF10", "IL10RB", "IL4R", "INHBC", "IL18R1", "CCL26", "IFNAR2", "TNFRSF14", "IL6", "LTA", "CXCR4", "GDF2", "IL6ST", "TNFRSF17", "CXCR6", "IL13", "BMPR1B", "TGFB3", "CRLF2", "IFNA16", "TNFSF10", "TNFRSF10A", "IL27", "TNFRSF25", "AMHR2", "BMP8B", "IL13RA2", "OSM", "CNTFR", "CCL3L3", "TNFSF13B", "TNFRSF11B", "TNFRSF11A", "TNFRSF10C", "EPO", "IL24", "CCL8", "ACVR1C", "BMP5", "IFNA8", "IFNK", "IL1A", "IL25", "TNFRSF21", "BMP7", "RELT", "CCR5", "IFNA21", "ACVR2B", "IL31", "IL36G", "IFNW1", "CXCL13", "BMP2", "IFNA7", "IL1B", "IFNA6", "LTB", "GDF5", "IL2RG", "CXCL17", "IFNL3", "GDF9", "IL2", "TSLP", "CD4", "CXCR2", "INHBA", "IL36A", "TNFSF4", "CCL4", "CCL15", "INHBB", "CX3CL1", "IL7R", "CCL7", "TNFSF18", "IL1F10", "AMH", "CCL17", "IL18RAP", "GDF7", "CSH1", "XCL1", "IL34", "IL17RA", "MSTN", "TNFSF15", "FAS", "IFNL1", "IFNGR1", "NGF", "IL5", "RELL2", "IL11RA", "CXCL9", "CCL21", "IL1RL2", "NODAL", "TNFSF11", "IL17A", "CXCR5", "TNF", "CXCL6", "IL2RA", "IL5RA", "IL17RC", "LEPR", "TNFRSF19", "CCL16", "CSF2RA", "IL19", "ACKR4", "GDF11", "IL27RA", "IL10RA", "CTF1", "TNFRSF13B", "TNFRSF12A", "IL17RB", "IL33", "EPOR", "IL36B", "IL1RN", "IL11", "IFNE", "GDF1", "ACVR2A", "IFNA1", "IL1R1", "IL12RB1", "XCR1", "CCL11", "IL36RN", "IFNA2", "IL21", "TNFSF12", "PF4V1", "CXCR3", "TNFRSF13C", "CCL23", "IL32", "BMP3", "IL9R", "IL22RA1", "CLCF1", "IFNG", "GH2", "CSF3R", "IL15RA", "LEP", "CXCL3", "TNFSF13", "IL15", "CCL4L2", "CSHL1", "CXCR1", "INHBE", "BMPR2", "IL17F", "EDAR", "IFNLR1", "CXCL11", "IL1RAP", "IL18", "NGFR", "TNFRSF10B", "CCL5", "OSMR", "GDF15", "IFNL2", "BMP4", "LIFR", "CCL13", "BMP8A", "TNFRSF1B", "IL16", "TNFRSF10D", "IFNA17", "CCL18", "LTBR", "CCL3", "CXCL5", "IL9", "IL17C", "IL23A", "IL6R", "CNTF", "TGFBR2", "IL17RE", "IL3RA", "CXCL2", "IL20RA", "CCR1", "IL12A", "TNFSF8", "ACVR1B", "CSF3", "ACKR3", "IL10", "IFNA13", "CSF1", "CCL4L1", "TNFRSF1A", "CCR7", "PPBP", "TNFSF14", "CD40LG", "CCR6", "IL13RA1", "IFNA4", "BMPR1A", "GDF3", "CCR10", "CCL25", "CCL22", "IL12RB2", "CCR3", "IL26", "CCL2", "MPL", "TGFBR1", "CD40", "PRLR", "IFNA14", "CCR2", "CCR4", "BMP6", "IFNGR2", "CCL27", "BMP15", "IL17D", "CCL24", "IL4", "TGFB2", "CSH2", "GDF6", "TNFSF9", "IFNA10", "THPO", "CCR8", "TNFRSF6B", "ACVR1", "PRL", "CSF2RB", "IL2RB", "IL31RA", "CXCL10", "EDA2R", "CD27", "CCL14", "CXCL14", "IFNAR1", "CXCL16", "IFNB1", "TNFRSF8", "IL21R", "CSF2", "CD70", "IL12B", "RELL1", "TGFB1", "INHA", "FASLG", "XCL2", "GH1", "IL3", "CX3CR1", "IL20RB", "IL37", "IL20", "IFNA5", "CXCL8", "IL1R2", "TNFRSF18", "TNFRSF9", "IL17B", "CXCL1", "CCL19", "EBI3", "TNFRSF4", "LIF", "CXCL12", "CCL3L1", "IL7", "EDA","ACKR2","IL17D","MYDGF","CCM2","BTF3P11","	TIMP1","EPX","MYDGF","IFNA1","ACKR2","CNOT6",]


    toll_like=["MAP2K6", "IKBKE", "TLR7", "IFNAR2", "IRAK1", "IL6", "STAT1", "IFNA16", "LBP", "TIRAP", "TLR6", "CCL3L3", "PIK3CD", "MAPK14", "IFNA8", "TICAM2", "IFNA21", "TAB1", "IFNA7", "IL1B", "IFNA6", "PIK3CB", "CCL4", "JAK1", "CTSK", "MAP2K7", "FADD", "MAPK11", "MAPK1", "CXCL9", "TLR1", "CHUK", "TNF", "AKT2", "MAPK8", "PIK3R1", "TYK2", "SPP1", "MAP3K8", "MAP2K4", "TLR3", "IFNA1", "RIPK1", "FOS", "IFNA2", "TLR9", "MAP3K7", "MAPK12", "IRF5", "MYD88", "MAP2K2", "PIK3R3", "CD86", "CD80", "CCL4L2", "CXCL11", "MAPK10", "TRAF6", "CCL5", "TOLLIP", "TBK1", "JUN", "TICAM1", "IRF3", "IFNA17", "MAP2K1", "NFKBIA", "CCL3", "CASP8", "PIK3R2", "MAPK3", "RAC1", "TRAF3", "IL12A", "IKBKB", "IFNA13", "CCL4L1", "PIK3CA", "IFNA4", "AKT3", "CD40", "STAT2", "IFNA14", "TLR5", "P3R3URF","PIK3R3", "AKT1", "MAPK9", "IKBKG", "TLR8", "CD14", "IFNA10", "TLR4", "RELA", "CXCL10", "IFNAR1", "IFNB1", "IRF7", "TAB2", "LY96", "IL12B", "NFKB1", "TLR2", "IRAK4", "IFNA5", "CXCL8", "MAP2K3", "MAPK13", "IRF9", "CCL3L1","CXXC1","RNASE1","IFNA1"]

    jak_stat=["IL22", "GHR", "IL23R", "IL10RB", "IL4R", "IFNAR2", "IL6", "CDKN1A", "STAM2", "STAT4", "IL6ST", "IL13", "HRAS", "STAT1", "CRLF2", "IFNA16", "IL27", "SOCS5", "IL13RA2", "OSM", "CNTFR", "STAT6", "PIK3CD", "EPO", "IL24", "IFNA8", "IFNK", "STAT5B", "IFNA21", "IL31", "CCND3", "IFNW1", "PIAS3", "IFNA7", "IL22RA2", "IFNA6", "IL2RG", "IFNL3", "SOCS4", "IL2", "TSLP", "EP300", "PIK3CB", "CISH", "SOCS7", "JAK1", "JAK2", "IL7R", "CREBBP", "CSH1", "GFAP", "SOS2", "PIAS2", "IFNL1", "MCL1", "IFNGR1", "IL5", "IL11RA", "AKT2", "IL2RA", "IL5RA", "PIK3R1", "LEPR", "TYK2", "SOCS1", "CSF2RA", "IL19", "IL27RA", "IL10RA", "CTF1", "EPOR", "IL11", "IFNE", "JAK3", "IFNA1", "IL12RB1", "IFNA2", "IL21", "PDGFRA", "IL9R", "PIAS4", "IL22RA1", "SOCS3", "CLCF1", "IFNG", "BCL2L1", "GH2", "CSF3R", "IL15RA", "LEP", "PIK3R3", "IL15", "CSHL1", "IFNLR1", "OSMR", "SOCS6", "EGF", "IFNL2", "LIFR", "SOCS2", "PDGFA", "IFNA17", "PIM1", "IL9", "IL23A", "IL6R", "CNTF", "AOX1", "IL3RA", "PIK3R2", "IL20RA", "IL12A", "PDGFB", "CSF3", "IL10", "CCND1", "IFNA13", "IL13RA1", "IFNA4", "GRB2", "PIK3CA", "PTPN2", "STAM", "AKT3", "IL12RB2", "IL26", "MTOR", "MPL", "STAT2", "PDGFRB", "PRLR", "IFNA14", "STAT3", "PIAS1", "P3R3URF","PIK3R3", "CCND2", "AKT1", "IFNGR2", "STAT5A", "IL4", "CSH2", "IFNA10", "PTPN6", "BCL2", "THPO", "PRL", "CSF2RB", "IL2RB", "IL31RA", "FHL1", "IFNAR1", "IFNB1", "IL21R", "CSF2", "IL12B", "MYC", "GH1", "IL3", "SOS1", "IL20RB", "IL20", "IFNA5", "PTPN11", "EGFR", "IRF9", "LIF", "RAF1","IL7","CDPF1","IL17D","MYDGF","TIMP1","EPX","SOCS4","LONP1","IFNA1","CFH","RNASE3"]

    nod_like=["PKN1", "PLCB4", "ITPR3", "PANX1", "NFKBIB", "IKBKE", "DEFA6", "RNF31", "IFNAR2", "IL6", "BIRC3", "MAVS", "CARD18", "STAT1", "YWHAE", "IFNA16", "TRAF5", "ITPR1", "OAS2", "TRAF2", "DEFB103B", "RHOA", "TRPV2", "MAP1LC3B", "MAPK14", "CASP1", "NLRP1", "VDAC1", "IFNA8", "GSDMD", "IFNA21", "TAB1", "PLCB3", "CASP12", "RIPK2", "IFNA7", "IL1B", "IFNA6", "GABARAPL2", "PYDC1", "ITPR2", "XIAP", "PKN2", "GBP3", "NAMPT", "TXN2", "JAK1", "BRCC3", "FADD", "MAPK11", "MAPK1", "GABARAPL1", "CHUK", "NLRP7", "TNF", "TP53BP1", "MAPK8", "PLCB1", "MFN2", "RBCK1", "PRKCD", "TYK2", "NLRP6", "PYDC2", "TXNIP", "HSP90AB1", "GBP7", "DEFB4A", "GBP1", "IFNA1", "RIPK1", "NLRC4", "GBP5", "IFNA2", "MFN1", "MAP3K7", "DEFB103A", "MAPK12", "RNASEL", "CARD16", "GPRC6A", "TRPM7", "GBP2", "P2RX7", "MYD88", "BCL2L1", "SHARPIN", "CXCL3", "MAP1LC3C", "RIPK3", "NLRX1", "ATG5", "TNFAIP3", "DEFA1", "DHX33", "MAPK10", "IL18", "DEFA1B", "CARD6", "TRAF6", "CCL5", "TBK1", "JUN", "GBP4", "CASP5", "CASR", "PLCB2", "NLRP3", "TICAM1", "CARD8", "TANK", "IFNA17", "IRF3", "VDAC3", "NFKBIA", "TRPM2", "NOD2", "AIM2", "NEK7", "DEFA4", "CASP8", "MAP1LC3B2", "BIRC2", "CXCL2", "PYCARD", "OAS3", "MAPK3", "TRAF3", "TRIP6", "TAB3", "IKBKB", "ANTXR1", "CARD9", "IFNA13", "DEFA3", "IFI16", "MEFV", "GABARAP", "IFNA4", "DEFB4B", "PYDC5", "CASP4", "CCL2", "STAT2", "IFNA14", "NAIP", "MAPK9", "IRGM", "VDAC2", "ATG16L1", "CTSB", "IKBKG", "NLRP12", "MCU", "PSTPIP1", "IFNA10", "OAS1", "BCL2", "TLR4", "TXN", "SUGT1", "RELA", "ANTXR2", "IFNAR1", "DNM1L", "IFNB1", "STING1", "IRF7", "CYBB", "TAB2", "NFKB1", "DEFA5", "IRAK4", "MAP1LC3BP1", "HSP90AA1", "GPSM3", "IFNA5", "CXCL8", "MAPK13", "CAMP", "CXCL1", "IRF9", "CYBA", "NOD1", "MAP1LC3A", "ERBIN", "ATG12","DEPDC1B","DEFA1A3","VDAC1P5","CHAMP1"] 


    t_cell=["NFKBIB", "LAT", "PPP3CB", "HRAS", "PPP3CA", "PRKCQ", "RHOA", "PIK3CD", "PPP2R2B", "LCK", "PPP2R5B", "MAPK14", "PAK3", "ZAP70", "PDCD1", "PAK6", "NFATC2", "NRAS", "PPP2R5E", "IL2", "PPP3R1", "CD4", "MALT1", "PIK3CB", "PPP2R5A", "RASGRP1", "PPP2R2C", "MAP2K7", "GSK3B", "CBLB", "PPP2CB", "CDC42", "MAPK11", "CD3G", "MAPK1", "SOS2", "IL5", "KRAS", "CHUK", "PPP3R2", "TNF", "AKT2", "MAPK8", "PIK3R1", "CD28", "PAK2", "PTPRC", "CD8A", "MAP3K8", "VAV2", "NFKBIE", "PAK4", "FOS", "ITK", "MAP3K7", "VAV3", "MAPK12", "CD3E", "PPP2R5C", "FYN", "IFNG", "MAP2K2", "PIK3R3", "PPP2R1B", "NFATC3", "CD8B2", "MAPK10", "VAV1", "PDPK1", "CD8B", "PAK5", "JUN", "PPP2R3C", "MAP3K14", "PPP2R1A", "MAP2K1", "NFKBIA", "CARD11", "LCP2", "PIK3R2", "MAPK3", "NCK1", "IKBKB", "TEC", "IL10", "CD40LG", "PPP2R3B", "PPP2R5D", "PIK3CA", "GRB2", "PLCG1", "ICOS", "NCK2", "AKT3", "CD3D", "P3R3URF","PIK3R3", "AKT1", "MAPK9", "PPP2R2D", "CTLA4", "PPP2R3A", "PAK1", "IKBKG", "IL4", "BCL10", "CD247", "PTPN6", "RELA", "DLG1", "GRAP2", "NFATC1", "CSF2", "NFKB1", "PPP3CC", "PPP2CA", "CDK4", "SOS1", "PTPN11", "PPP2R2A", "MAPK13", "BUB1B-PAK6", "RAF1","SPNS1", "ARHGEF7","MMAB","PKN1","DLGAP5","RNASE3"]


    il_17=["CCL20", "IKBKE", "IL6", "IL13", "USP25", "TRAF5", "TRAF2", "MUC5B", "MAPK14", "MMP13", "IL25", "S100A7A", "IL1B", "ANAPC5", "MAPK15", "CCL7", "GSK3B", "CCL17", "FADD", "MAPK11", "FOSL1", "IL17RA", "MAPK1", "PTGS2", "IL5", "CHUK", "IL17A", "TNF", "CXCL6", "HSP90B1", "MAPK8", "TRAF4", "IL17RB", "HSP90AB1", "DEFB4A", "CCL11", "FOS", "MAP3K7", "MAPK12", "SRSF1", "IFNG", "CXCL3", "IL17F", "TNFAIP3", "MAPK10", "MMP3", "TRAF6", "TBK1", "JUN", "JUND", "NFKBIA", "MUC5AC", "MAPK6", "CXCL5", "IL17C", "TRAF3IP2", "FOSB", "CASP8", "IL17RE", "CXCL2", "CEBPB", "MAPK3", "TRAF3", "TAB3", "IKBKB", "CSF3", "DEFB4B", "MMP9", "MMP1", "S100A9", "CCL2", "ELAVL1", "LCN2", "MAPK9", "MAPK7", "IL17D", "S100A8", "IKBKG", "CASP3", "IL4", "RELA", "CXCL10", "S100A7", "CSF2", "TAB2", "NFKB1", "MAPK4", "HSP90AA1", "CXCL8", "MAPK13", "IL17B", "CXCL1", "TRADD", "IL17RC","TANK","MYDGF"]



    pi3kakt=["CSF1R", "ITGA9", "GNG11", "COL4A5", "PSPN", "LAMC1", "VEGFD", "THBS1", "CREB3L3", "TNR", "GNGT1", "EPHA2", "IFNA8", "CASP9", "EFNA5", "IFNA6", "EIF4E2", "FGF4", "F2R", "JAK1", "KIT", "MAPK1", "FN1", "RBL2", "PTEN", "ITGA6", "LPAR5", "COMP", "GNG12", "HSP90B1", "PIK3R1", "ITGB5", "JAK3", "PRKCA", "IGH", "IFNA2", "ITGA2B", "PHLPP1", "FGF7", "RET", "PPP2R1B", "THEM4", "FLT3LG", "ANGPT4", "VEGFC", "VEGFA", "FGF22", "LAMB4", "LAMA1", "MAP2K1", "ITGA7", "INS", "PIK3R2", "MAPK3", "CSF1", "PRKAA1", "GRB2", "AKT3", "PRKAA2", "COL6A3", "MLST8", "MTOR", "FGF1", "PRLR", "GNG13", "CCND2", "AKT1", "PPP2R2D", "PPP2R3A", "C8orf44","SGK3", "FLT1", "GYS2", "IL4", "ARTN", "EFNA4","EFNA3", "PRL", "KITLG", "G6PC3", "NFKB1", "ITGA8", "TGFA", "IL3", "SOS1", "PPP2R2A", "TNXB", "PKN1", "EIF4E1B", "ITGB4", "BDNF", "CREB3L4", "COL6A5", "IFNAR2", "COL9A1", "IL6", "GNG2", "LPAR2", "ITGA2", "HRAS", "YWHAE", "OSM", "FGF20", "PIK3CG", "LPAR3", "PPP2R2B", "PPP2R5B", "TSC2", "COL6A1","MDM2", "CCND3", "FGFR1"," NRAS", "GNB3", "PIK3CB", "IGF2", "FGF21", "CCNE2", "CSH1", "NTRK1", "MCL1", "NGF", "RPS6"," GNG10", "GNG5", "ITGA4", "PCK2", "PIK3AP1", "ERBB3", "YWHAQ"," EIF4EBP1", "PPP2R5C", "ERBB4", "ITGA3", "ITGB6", "FGF19", "BRCA1", "NGFR", "OSMR", "PTK2", "KDR", "IL3RA", "FGF18", "IKBKB", "CCND1", "ANGPT1"," COL4A6", "IFNA4", "ERBB2", "COL2A1", "COL4A1", "NTRK2", "GNG4", "BCL2L11", "YWHAZ", "LAMB3", "YWHAH", "YWHAB", "BCL2", "IGF1R", "IL2RB", "DDIT4", "ITGB7", "TLR2", "ITGA1", "ANGPT2", "GH1", "FGF10", "HSP90AA1", "CHRM2", "CRTC2", "PIK3R6", "ATF6B", "RAF1", "PCK1", "IL4R", "IRS1", "MAGI1", "EFNA2", "PDGFC", "COL4A4", "IFNA16", "YWHAG", "LPAR6", "PIK3R5", "TNN", "LAMA5", "EPO", "FGF6", "THBS2", "IL2RG", "IL7R", "CREB3L1", "PPP2CB", "RHEB", "PDGFD", "CREB5", "FGF2", "CHUK", "IL2RA", "CREB1", "COL9A2", "EPOR", "MAGI2", "G6PC2", "LAMA4", "GNB5", "G6PC1", "GH2", "MAP2K2", "COL1A1", "RPS6KB1", "VEGFB", "CSHL1", "FGF23", "RPS6KB2", "IGF1", "PDPK1", "PPP2R3C", "CDKN1B", "IFNA17", "FGFR2", "PPP2R1A", "GNB1", "NRTN", "FLT3", "COL4A3", "NTF4", "CSF3", "SGK1", "FGF8", "INSR", "PPP2R5D", "EREG", "P3R3URF-PIK3R3", "CDK2", "IKBKG", "SYK", "CSH2", "LPAR1", "IFNA10", "ITGB1", "IFNAR1", "IFNB1", "FGFR4", "EIF4B", "COL6A6", "CDK4", "ITGB3", "COL6A2", "IFNA5", "CREB3", "EGFR", "CREB3L2", "GYS1", "IL7", "ITGAV", "FLT4", "GHR", "LAMC3", "GNGT2", "GNG7", "SGK2", "CCNE1", "CDKN1A", "PKN3", "THBS4", "FGF16", "GDNF", "TCL1B", "VTN", "PIK3CD", "TSC1", "CDK6", "LAMA3", "CDC37", "IFNA21", "NR4A1", "IFNA7", "PPP2R5E", "TP53", "IL2", "PKN2", "PPP2R5A", "EFNA1", "JAK2", "PPP2R2C", "GSK3B", "TEK", "NTF3", "FGF5", "SOS2", "FGFR3", "KRAS", "LAMC2", "LPAR4", "AKT2", "RELN", "ITGB8", "SPP1", "COL9A3", "HSP90AB1", "GNG3", "STK11", "IFNA1", "RPTOR", "PDGFRA", "PGF", "BCL2L1", "CSF3R", "PIK3R3", "LAMB2", "CHRM1", "THBS3", "GNB4", "EGF", "FOXO3", "GNG8", "PDGFA", "MTCP1", "COL1A2", "ITGA11", "RXRA", "ATF2", "IL6R", "RAC1", "LAMA2", "PDGFB", "CHAD", "IFNA13", "EFNA3", "PPP2R3B", "TCL1A", "IBSP", "PIK3CA", "SGK3", "ITGA5", "PDGFRB", "IFNA14", "FGF17", "EFNA4", "ITGA10", "LAMB1", "TLR4", "MYB", "RELA", "AREG", "CD19", "HGF", "NOS3", "TNC", "GNB2", "PHLPP2", "ATF4", "MYC", "FASLG", "PPP2CA", "VWF", "MET", "COL4A2", "EIF4E", "FGF3", "BAD", "FGF9","CCM2","TIMP1","EPX","RHEBP1","FGF13","CREB3L4","SGK3","CDPF1","LAMA4","CXXC1","LAMC1","SH2D1A","GDNF","RNASE1","IFNA1","IL6","SOS1","NANOS3","SLTM" ]



    mapk=["RRAS2", "NF1", "RPS6KA2"," MAP3K1", "FLT4", "PPP5C"," DUSP4", "HSPA1B", "CSF1R", "NFKB2", "TAOK1", "BDNF", "CACNA2D3", "CACNA1B", "MAP3K12", "GADD45G", "RRAS", "MAP2K6", "PSPN", "DUSP16", "IRAK1", "CACNA1A", "EFNA2"," MKNK2", "MAP3K4", "DUSP6", "PPP3CB", "HRAS", "PDGFC", "SRF", "MAPKAPK5", "TGFB3", "CACNA1E", "MAPKAPK2", "VEGFD", "PPP3CA", "PTPN5", "MKNK1", "RASGRP3", "FGF16", "GDNF", "TRAF2", "RAPGEF2", "CACNA2D1", "FGF20", "EPHA2", "PRKCG", "CACNG6", "MAPK14", "RASGRF1", "MAP4K3", "PRKACA", "LAMTOR3", "MAP4K1","EFNA5", "IL1A", "FGF6", "NR4A1", "TAB1", "FGFR1", "NRAS", "DUSP1", "IL1B", "TP53", "MAP3K11", "FGF4", "PPP3R1", "CRK", "RASGRF2", "GADD45B", "RASGRP1", "IGF2", "HSPA1L", "CACNA1S", "MAP4K4", "EFNA1", "FGF21", "MAP3K3", "KIT", "MAP2K7", "CACNB3", "PRKCB", "PRKACB", "PLA2G4E", "CDC42", "PPM1A", "MAPK11", "NTF3", "NTRK1", "TEK", "PDGFD", "MAPK1", "FGF5", "FAS", "SOS2", "GADD45A", "FGFR3", "MAP3K10", "NGF", "CACNA1I", "KRAS", "RPS6KA5", "ELK1", "FGF2", "CHUK", "HSPA6", "NLK", "PPP3R2", "GNG12", "DUSP3", "RPS6KA3", "AKT2", "TNF", "MAPK8", "CACNA1C", "CACNB2", "DUSP7", "CACNG1", "HSPA8", "PLA2G4C", "PLA2G4D", "PAK2", "DDIT3", "PPM1B", "MAP3K8", "RPS6KA6"," ERBB3", "VEGFC"," MAP3K5", "MAP2K4", "CACNA1F", "IL1R1", "PRKCA", "FOS", "PDGFRA", "CDC25B", "MAP3K7", "CACNG5", "MAPK12"," MAPK8IP2", "TAOK3", "FGF7", "ERBB4", "MYD88", "PGF", "MAP2K2", "RET","BRAF", "CACNB4", "FGF19", "RAP1B", "VEGFB", "PTPRR", "NFATC3", "FGF23", "ELK4", "FLT3LG", "IL1RAP", "CACNG4", "IGF1", "STMN1", "DUSP5", "NGFR", "PLA2G4F", "MAPK10", "RAC3", "TRAF6", "JUN", "ANGPT4", "PTPN7", "MAPK8IP3", "DUSP10", "EGF", "FGF22", "KDR", "VEGFA", "PLA2G4B", "CACNG2", "PDGFA", "RAP1A", "CACNA1H"," MAP3K14"," MAP3K6", "FGFR2", "JUND", "MAP2K1", "HSPB1", "NRTN", "FLT3", "CACNA1G", "ATF2", "TGFBR2", "MAP3K2", "TRADD", "MAPT", "INS", "MAPK3", "RAC1", "MRAS", "PDGFB", "IKBKB", "CACNG3", "NTF4", "CACNA1D", "FGF18", "CSF1", "EFNA3", "FGF8", "GNA12", "ANGPT1", "INSR", "RASA2", "CACNG8", "TNFRSF1A", "MAP3K13", "CACNA2D2", "ERBB2", "GRB2", "DAXX", "PLA2G4A", "AKT3", "MAX", "TGFBR1", "RPS6KA1", "EREG", "PDGFRB", "FGF1", "FLNA", "MAP3K21", "PRKACG", "RASGRP2", "CRKL", "NTRK2", "MEF2C", "AKT1", "MAP3K9", "MAPK9", "FGF17", "MAPK7", "EFNA4", "ECSIT", "FLT1", "PAK1", "STK4", "CASP3", "RAC2", "IKBKG", "TAOK2", "TGFB2", "RASA1", "CD14", "RPS6KA4", "DUSP9", "HSPA1A", "EFNA4-EFNA3", "IGF1R", "ARTN", "RELA", "KITLG", "ARRB1", "MAPKAPK3", "AREG", "ARAF", "MAPK8IP1", "CACNG7", "JMJD7-PLA2G4B", "HGF", "NFATC1", "MAP2K5", "FGFR4", "ARRB2", "TAB2", "NFKB1", "PPP3CC", "ATF4", "IRAK4", "ANGPT2", "FASLG", "MYC", "TGFA", "DUSP2", "MECOM", "TGFB1", "FGF10", "SOS1", "STK3", "CACNA2D4", "RASGRP4", "RELB", "MET", "CACNB1", "MAP2K3", "MAPK13", "EGFR", "MAP4K2", "FGF3", "RAF1", "DUSP8", "FGF9", "MAP3K20", "DUSP12","FASN","KCNH4","FGF13","GPI","NCOA3","GDNF","RNASE1","PKN1","IL6","SOS1","STK24","SLTM","RNASE3","TANK"]


    tgf_beta=["NRROS", "ID3", "CHRD", "BAMBI", "RBL1", "NOG", "INHBC", "SMURF2", "SIN3A", "E2F4", "CUL1", "TGIF2", "SMAD1", "BMPR1B", "SMAD3", "TGFB3", "NBL1", "THBS1", "AMHR2", "BMP8B", "E2F5", "CDKN2B", "LRRC32", "RHOA", "ACVR2B", "BMP5", "BMP7", "ACVR1C", "BMP2", "GDF5", "HDAC1", "TGIF1", "SKI", "LEFTY1", "DCN", "EP300", "INHBA", "SKP1", "SMAD2", "INHBB", "ZFYVE9", "AMH", "ID4", "PPP2CB", "CREBBP", "GDF7", "LTBP1", "PITX2", "SP1", "MAPK1", "FBN1", "NODAL", "HFE", "TNF", "SMAD6", "EMP3", "RBX1", "GREM1", "ID1", "ACVR2A", "FST", "SMAD9", "RGMA", "HJV", "IFNG", "RGMB", "HAMP", "RPS6KB1", "PPP2R1B", "INHBE", "BMPR2", "RPS6KB2", "ZFYVE16", "ID2", "SMAD4", "SMURF1", "SKIL", "BMP4", "BMP8A", "PPP2R1A", "FMOD", "NEO1", "TGFBR2", "MAPK3", "GREM2", "ACVR1B", "BMPR1A", "TMEM53", "MICOS10-NBL1", "TGFBR1", "TFDP1", "SMAD5", "LEMD3", "BMP6", "NCOR1", "TFRC", "THSD4", "TGFB2", "GDF6", "ACVR1", "IGSF1", "TGFB1", "MYC", "PPP2CA", "SMAD7", "LEFTY2", "ROCK1", "HDAC2", "TF", "TFR2","GARS1","HHAT","F3"]


    lit_genes = ["POTED", "ANKRD26", "PSG4", "GTF3C2", "USH2A", "EOMES", "UGT2B7", "HLA-DQB1", "ABCB5", "ZNF610", "HSD17B4", "TRIM66", "LMO7", "PPP1R13B", "NID2", "HNRNPCL2", "HRNR", "MYO7A", "KIR3DL3", "HLA-B",         
    "TNRC6A", "CEMIP", "PREX2", "GATA5", "EFCC1", "ODF3", "F5", "SORCS2", "RP1L1", "MUC3A", "PLOD3", "FAM240A", "KRTAP10-1", "THSD7A","HLA-DRB1", "IL6ST", "IL17RC", "VEGFB", "ITPKC", "CASP3", "ORAI", "MYH11", "SMAD9", "FCRLA", "PTGER4", "IL17F", "CARD11", "SIGLEC10", "FGFR4", "IL31RA", "FNDC1", "MMP8", "FOXN1", "TTI1", "MYH14", "RBP3", "CD36","CIMAP1A", "MYH7B","E2F1"]
    return (
        Cytokine,
        il_17,
        jak_stat,
        lit_genes,
        mapk,
        nod_like,
        pi3kakt,
        t_cell,
        tgf_beta,
        toll_like,
        vegf,
    )


@app.cell
def _(
    Cytokine,
    il_17,
    jak_stat,
    lit_genes,
    mapk,
    nod_like,
    pi3kakt,
    t_cell,
    tgf_beta,
    toll_like,
    vegf,
):
    kegg_genes = []
    for g in vegf + Cytokine + toll_like + jak_stat +nod_like + t_cell + il_17 + pi3kakt + mapk + tgf_beta + lit_genes :
        kegg_genes.append(g)
    return (kegg_genes,)


@app.cell
def _(cleaned_mt, hl, kegg_genes):
    # Hail'e literal set olarak aktar

    kegg_genes_final_mis = hl.literal(set(kegg_genes))

    # MatrixTable'ı sadece VEGF pathway genleriyle filtrele
    kegg_genes_mis = cleaned_mt.filter_rows(
        kegg_genes_final_mis.contains(cleaned_mt.gene_symbol)
    )

    print("Orijinal varyant sayısı:", cleaned_mt.count_rows())
    print("All pathway varyant sayısı:", kegg_genes_mis.count_rows())

    # İstersen bu yeni MatrixTable'ı kaydet
    kegg_genes_mis.write("kegg_mis.py_variants.mt", overwrite=True)
    return (kegg_genes_mis,)


@app.cell
def _(kegg_genes_mis):
    kegg_genes_mis.entries().show(100)
    return


@app.cell
def _(hl, kegg_genes_mis):
    # Since they are in the VEP CSQ, we defined them and saved separately
    csq_fields = [
        "Allele","Consequence","IMPACT","SYMBOL","Gene","Feature_type","Feature","BIOTYPE",
        "EXON","INTRON","HGVSc","HGVSp","cDNA_position","CDS_position","Protein_position",
        "Amino_acids","Codons","Existing_variation","ALLELE_NUM","DISTANCE","STRAND","FLAGS",
        "VARIANT_CLASS","MINIMISED","SYMBOL_SOURCE","HGNC_ID","CANONICAL","MANE","MANE_SELECT",
        "MANE_PLUS_CLINICAL","TSL","APPRIS","CCDS","ENSP","SWISSPROT","TREMBL","UNIPARC",
        "UNIPROT_ISOFORM","GENE_PHENO","SIFT","PolyPhen","DOMAINS","miRNA","HGVS_OFFSET",
        "AF","AFR_AF","AMR_AF","EAS_AF","EUR_AF","SAS_AF","gnomADe_AF","gnomADe_AFR_AF",
        "gnomADe_AMR_AF","gnomADe_ASJ_AF","gnomADe_EAS_AF","gnomADe_FIN_AF","gnomADe_MID_AF",
        "gnomADe_NFE_AF","gnomADe_REMAINING_AF","gnomADe_SAS_AF","gnomADg_AF","gnomADg_AFR_AF",
        "gnomADg_AMI_AF","gnomADg_AMR_AF","gnomADg_ASJ_AF","gnomADg_EAS_AF","gnomADg_FIN_AF",
        "gnomADg_MID_AF","gnomADg_NFE_AF","gnomADg_REMAINING_AF","gnomADg_SAS_AF","MAX_AF",
        "MAX_AF_POPS","CLIN_SIG","SOMATIC","PHENO","PUBMED","MOTIF_NAME","MOTIF_POS",
        "HIGH_INF_POS","MOTIF_SCORE_CHANGE","TRANSCRIPTION_FACTORS","am_class","am_pathogenicity"
    ]

    kegg_genes_mis_csq = kegg_genes_mis.annotate_rows(
        CSQ_struct = hl.struct(**{
            field: kegg_genes_mis.info.CSQ[0].split('\|')[i] for i, field in enumerate(csq_fields)
        })
    )

    kegg_genes_mis_csq.CSQ_struct.show(5)
    return (kegg_genes_mis_csq,)


@app.cell
def _(kegg_genes_mis_csq):
    kegg_genes_mis_csq.rows().show()
    return


@app.cell
def _(
    Cytokine,
    hl,
    il_17,
    jak_stat,
    kegg_genes_mis_csq,
    lit_genes,
    mapk,
    nod_like,
    pi3kakt,
    t_cell,
    tgf_beta,
    toll_like,
    vegf,
):
    def write_setid_locus(matrix_table, gene_sets_dict, output_file):
        """
        Her gen seti için:
        set adı -> alt alta chr:pos:ref:alt
        """
        rows = matrix_table.rows()

        with open(output_file, "w") as f:
            for set_name, gene_list in gene_sets_dict.items():
                gene_set = set(gene_list)
                filtered = rows.filter(hl.literal(gene_set).contains(rows.gene_symbol))

                # locus ve alleles’i tuple olarak topla
                loci_tuples = filtered.aggregate(
                    hl.agg.collect((filtered.locus.contig, filtered.locus.position,
                                    filtered.alleles[0], filtered.alleles[1]))
                )

                # Python tarafında string oluştur
                for contig, pos, ref, alt in loci_tuples:
                    if None not in (contig, pos, ref, alt):
                        f.write(f"{set_name} {contig}:{pos}:{ref}:{alt}\n")



    # Örnek kullanım
    gene_sets = {
        "VEGF": vegf,
        "TGF_BETA": tgf_beta,
        "PI3KAKT": pi3kakt,
        "CYTOKINE": Cytokine,
        "TOLL_LIKE": toll_like,
        "JAK_STAT": jak_stat,
        "NOD_LIKE": nod_like,
        "T_CELL": t_cell,
        "IL_17": il_17,
        "MAPK": mapk,
        "LIT_GENES": lit_genes
    }

    write_setid_locus(kegg_genes_mis_csq, gene_sets, "plink/kegg_genes_mis.SetID")
    return


@app.cell
def _(kegg_genes_mis_csq):
    kegg_genes_mis_csq.describe()
    return


@app.cell
def _(hl, kegg_genes_mis_csq):
    #familyhistorychd
    hl.export_plink(
        dataset=kegg_genes_mis_csq,
        output='plink/kegg_genes_mis_famhis_chd',          # 1 = male, 2 = female, 0 = unknown
        pheno=kegg_genes_mis_csq.pheno.Family_history_of_CHD_status
    )


    #caa status 
    hl.export_plink(
        dataset=kegg_genes_mis_csq,
        output='plink/kegg_genes_mis_caa',          # 1 = male, 2 = female, 0 = unknown
        pheno=kegg_genes_mis_csq.pheno.CAA_status
    )

    # Diagnostic_Age_Status
    hl.export_plink(
        dataset=kegg_genes_mis_csq,
        output='plink/kegg_genes_mis_diagnostic_age',          # 1 = male, 2 = female, 0 = unknown
        pheno=kegg_genes_mis_csq.pheno.Diagnostic_Age_Status
    )

    #Sequelae_Status
    hl.export_plink(
        dataset=kegg_genes_mis_csq,
        output='plink/kegg_genes_mis_Sequelae_Status',          # 1 = male, 2 = female, 0 = unknown
        pheno=kegg_genes_mis_csq.pheno.Sequelae_Status
    )

    #Family_History_Status
    hl.export_plink(
        dataset=kegg_genes_mis_csq,
        output='plink/kegg_genes_mis_Family_History_Status',          # 1 = male, 2 = female, 0 = unknown
        pheno=kegg_genes_mis_csq.pheno.Family_History_Status
    )


    #Consanguineous_marriage_status
    hl.export_plink(
        dataset=kegg_genes_mis_csq,
        output='plink/kegg_genes_mis_Consanguineous_marriage_status',          # 1 = male, 2 = female, 0 = unknown
        pheno=kegg_genes_mis_csq.pheno.Consanguineous_marriage_status
    )

    #Degree_of_CM
    hl.export_plink(
        dataset=kegg_genes_mis_csq,
        output='plink/kegg_genes_mis_Degree_of_CM',          # 1 = male, 2 = female, 0 = unknown
        pheno=kegg_genes_mis_csq.pheno.Degree_of_CM
    )

    #Age_at_diagnosis_years
    hl.export_plink(
        dataset=kegg_genes_mis_csq,
        output='plink/kegg_genes_mis_Age_at_diagnosis_years',          # 1 = male, 2 = female, 0 = unknown
        pheno=kegg_genes_mis_csq.pheno.Age_at_diagnosis_years
    )

    #KD_in_siblings_status
    hl.export_plink(
        dataset=kegg_genes_mis_csq,
        output='plink/kegg_genes_mis_KD_in_siblings_status',          # 1 = male, 2 = female, 0 = unknown
        pheno=kegg_genes_mis_csq.pheno.KD_in_siblings_status
    )
    return


@app.cell
def _(hl, kegg_genes_mis_csq):
    kegg_genes_mis_csq2 = kegg_genes_mis_csq.annotate_cols(
        case_numeric = hl.case()
            .when(kegg_genes_mis_csq.case_control_status == "case", 2)
            .when(kegg_genes_mis_csq.case_control_status == "control", 1)
            .default(0)
    )

    hl.export_plink(dataset=kegg_genes_mis_csq2,
        output='plink/kegg_genes_mis_case_control',          # 1 = male, 2 = female, 0 = unknown
        pheno=kegg_genes_mis_csq2.case_numeric
    )
    return


@app.cell
def _(hl, kegg_genes_mis_csq):
    kegg_genes_mis_csq3 = kegg_genes_mis_csq.annotate_cols(
        Sex_numeric = hl.case()
            .when(kegg_genes_mis_csq.pheno.Sex == "male", 1)
            .when(kegg_genes_mis_csq.pheno.Sex == "female", 2)
            .default(0)
    )

    hl.export_plink(
        dataset=kegg_genes_mis_csq3,
        output='plink/kegg_genes_mis_Sex',
        pheno=kegg_genes_mis_csq3.Sex_numeric
    )
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
