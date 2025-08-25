import marimo

__generated_with = "0.14.9"
app = marimo.App(width="medium")


@app.cell
def _():
    import  marimo as mo
    return


@app.cell
def _():
    import hail as hl

    hl.init()
    return (hl,)


@app.cell
def _():
    #mt = hl.read_matrix_table('kawasaki_filtered.mt')
    return


@app.cell
def _():
    #mt_nogap = mt.filter_rows(mt.alleles[1] != '*')

    return


@app.cell
def _():
    #result = hl.vep(mt_nogap, 'vep-configuration.json')
    return


@app.cell
def _():
    #result.write('vep_annotated.mt', overwrite=True)
    return


@app.cell
def _():
    #resultu kaydettik ya direkt açmak için bunu kullan artık 
    #result = hl.read_matrix_table('vep_annotated.mt')
    return


@app.cell
def _(hl):
    result= hl.read_matrix_table('vep_annotated.mt')
    return (result,)


@app.cell
def _(result):
    result.describe()
    return


@app.cell
def _(result):
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns

    print("--- En Şiddetli Varyant Etkisi Dağılımı ---")

    # Step 1: Select the 'most_severe_consequence' directly from the MatrixTable's rows
    # and then convert it to a Table.
    # This approach ensures you are working within the correct context.
    consequence_data_table = result.select_rows(
        consequence=result.vep.most_severe_consequence
    ).rows()

    return consequence_data_table, plt, sns


@app.cell
def _(consequence_data_table):
    # Hail Tablosunu Pandas DataFrame'e dönüştürün
    df_consequence = consequence_data_table.to_pandas()

    print("Pandas DataFrame oluşturuldu.")
    return (df_consequence,)


@app.cell
def _(df_consequence):
    # Adım 2: Yüksek Etkili Varyant Tiplerine Göre Filtreleme
    # Yüksek etkili varyant tipleri listesi
    high_impact_consequences = [
        'missense_variant',
        'frameshift_variant',
        'stop_gained',
        'splice_donor_variant',
        'splice_acceptor_variant',
        'start_lost',
        'stop_lost',
        'regulatory_region_variant',
        'intron_variant'
    ]

    # DataFrame'i (bu durumda df_consequence) yüksek etkili varyant tiplerine göre filtrele
    df_high_impact = df_consequence[df_consequence['consequence'].isin(high_impact_consequences)].copy()
    print(f"\nYüksek etkili varyant tiplerine sahip toplam varyant sayısı: {len(df_high_impact)}")
    print("Yüksek etkili varyantların ilk 5 satırı:")
    print(df_high_impact.head())
    # Adım 2: Yüksek Etkili Varyant Tiplerine Göre Filtreleme
    # Yüksek etkili varyant tipleri listesi
    return (df_high_impact,)


@app.cell
def _(df_high_impact, plt, sns):


    # Yüksek etkili varyantların 'consequence' sütununun değer sayımlarını al
    consequence_counts_high_impact = df_high_impact['consequence'].value_counts()

    print("\nYüksek Etkili Varyant Tiplerinin Dağılımı:")
    print(consequence_counts_high_impact)

    # Bar grafiği ile görselleştirme
    plt.figure(figsize=(10, 6))
    sns.barplot(x=consequence_counts_high_impact.index, y=consequence_counts_high_impact.values, palette='viridis')
    plt.title('Yüksek Etkili Varyant Tiplerinin Dağılımı')
    plt.xlabel('Varyant Etkisi (Consequence)')
    plt.ylabel('Varyant Sayısı')
    plt.xticks(rotation=45, ha='right') # Etiketleri daha iyi okumak için döndür
    plt.tight_layout() # Grafiğin düzenini iyileştir
    plt.show()
    return


@app.cell
def _(df_consequence, df_high_impact):
    print("df_consequence DataFrame Sütunları:")
    print(df_consequence.columns.tolist())

    print("\ndf_consequence DataFrame İlk 5 Satır:")
    print(df_consequence.head())

    # Veya sadece df_high_impact için:
    print("\ndf_high_impact DataFrame Sütunları:")
    print(df_high_impact.columns.tolist())

    print("\ndf_high_impact DataFrame İlk 5 Satır:")
    print(df_high_impact.head())
    return


@app.cell
def _(consequence_data_table):
    # Eğer consequence_data_table bir Hail Table ise
    print(consequence_data_table.describe())
    print(consequence_data_table.row.show()) # İlk birkaç satır ve sütunları gösterebilir
    return


@app.cell
def _(df_consequence, df_high_impact, plt, sns):


    def get_variant_type(alleles):
        ref = alleles[0]
        alt = alleles[1]
        if len(ref) == 1 and len(alt) == 1:
            return 'SNP'
        elif len(ref) > len(alt):
            return 'Deletion'
        elif len(alt) > len(ref):
            return 'Insertion'
        return 'Complex' # Veya diğer tipler

    df_consequence['variant_type'] = df_consequence['alleles'].apply(get_variant_type)

    variant_type_counts = df_consequence['variant_type'].value_counts()
    print("\nVaryant Tipi Dağılımı (SNP/InDel):")
    print(variant_type_counts)

    plt.figure(figsize=(7, 5))
    # BURADAKİ HATA DÜZELTİLDİ: sns.pie yerine plt.pie kullanıldı
    plt.pie(variant_type_counts, labels=variant_type_counts.index, autopct='%1.1f%%', startangle=90, colors=sns.color_palette('pastel'))
    plt.title('Varyant Tipi Dağılımı')
    plt.ylabel('') # Pasta grafiğinde Y ekseni etiketi genellikle gizlenir
    plt.tight_layout() # Grafiğin düzenini iyileştirir
    plt.show()

    # Yüksek etkili varyantlar için de benzer analizi yapabiliriz
    df_high_impact['variant_type'] = df_high_impact['alleles'].apply(get_variant_type)
    high_impact_variant_type_counts = df_high_impact['variant_type'].value_counts()
    print("\nYüksek Etkili Varyant Tipleri Dağılımı (SNP/InDel):")
    print(high_impact_variant_type_counts)

    # İsteğe bağlı olarak yüksek etkili varyant tipleri için de pasta grafiği ekleyelim
    plt.figure(figsize=(7, 5))
    plt.pie(high_impact_variant_type_counts, labels=high_impact_variant_type_counts.index, autopct='%1.1f%%', startangle=90, colors=sns.color_palette('dark'))
    plt.title('Yüksek Etkili Varyant Tiplerinin Dağılımı')
    plt.ylabel('')
    plt.tight_layout()
    plt.show()
    return


if __name__ == "__main__":
    app.run()
