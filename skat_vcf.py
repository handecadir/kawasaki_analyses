import marimo

__generated_with = "0.17.7"
app = marimo.App(width="medium")


@app.cell
def _():
    import pandas as pd
    import os
    import csv

    def create_franklin_ready_vcf(file_list, output_filename='franklin_upload_final.vcf'):
        dfs = []

        print("Dosyalar taranıyor...")
        for f in file_list:
            if not os.path.exists(f):
                continue
            try:
                df = pd.read_csv(f)
                df.columns = df.columns.str.strip()
                dfs.append(df)
            except Exception as e:
                print(f"Hata: {f} - {e}")

        if not dfs:
            return

        combined_df = pd.concat(dfs, ignore_index=True)

        # Temel Dönüşümler
        combined_df['Chr'] = combined_df['Chr'].astype(str).apply(lambda x: x if str(x).startswith('chr') else f'chr{x}')
        combined_df['Pos'] = combined_df['Pos'].astype(int)
        combined_df['Ref'] = combined_df['Ref'].astype(str)
        combined_df['Alt'] = combined_df['Alt'].astype(str)

        if 'SNP_ID' in combined_df.columns:
            combined_df['ID'] = combined_df['SNP_ID'].fillna('.')
        else:
            combined_df['ID'] = '.'

        combined_df['QUAL'] = '.'
        combined_df['FILTER'] = 'PASS'

        # INFO Sütunu
        def create_info(row):
            parts = []
            if 'Pathway' in row.index and pd.notna(row['Pathway']):
                val = str(row['Pathway']).replace(" ", "_").replace(";", "_").replace("=", "_")
                parts.append(f"Pathway={val}")
            if 'P.value' in row.index and pd.notna(row['P.value']):
                parts.append(f"Pvalue={row['P.value']}")
            return ";".join(parts) if parts else "."

        combined_df['INFO'] = combined_df.apply(create_info, axis=1)

        # --- YENİ EKLENEN KISIM: FORMAT ve SAMPLE ---
        # GT: Genotype (0/1 - Heterozigot kabul edelim)
        # AD: Allele Depth (Referans 10, Alternatif 10 okuma gibi uyduralım)
        # DP: Total Depth (20 okuma)
        # GQ: Genotype Quality (99 - Yüksek kalite)
        # PL: Phred-scaled likelihoods (0,30,300 gibi standart değerler)

        combined_df['FORMAT'] = 'GT:AD:DP:GQ:PL'
        combined_df['SAMPLE_DATA'] = '0/1:10,10:20:99:0,30,300'

        # Sütun Seçimi
        vcf_df = combined_df[['Chr', 'Pos', 'ID', 'Ref', 'Alt', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE_DATA']].copy()

        # Sütun İsimleri (SAMPLE_DATA yerine örneğe bir isim verelim: 'MySample')
        vcf_df.columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'MySample']

        # Sıralama
        def get_chrom_order(chrom):
            c = chrom.lower().replace('chr', '')
            if c.isdigit(): return int(c)
            if 'x' in c: return 98
            if 'y' in c: return 99
            if 'm' in c: return 100
            return 101

        vcf_df['chrom_order'] = vcf_df['#CHROM'].apply(get_chrom_order)
        vcf_df = vcf_df.sort_values(by=['chrom_order', 'POS'])
        vcf_df = vcf_df.drop(columns=['chrom_order'])

        # Yazma İşlemi
        print(f"VCF oluşturuluyor: {output_filename}")

        with open(output_filename, 'w', newline='') as f:
            # VCF Header
            f.write('##fileformat=VCFv4.2\n')
            f.write('##source=GenotypeMatrix\n')

            # Contig Tanımları
            chroms = [str(i) for i in range(1, 23)] + ['X', 'Y', 'M', 'MT']
            for c in chroms:
                f.write(f'##contig=<ID=chr{c}>\n')

            # INFO Tanımları
            f.write('##INFO=<ID=Pathway,Number=1,Type=String,Description="Pathway info">\n')
            f.write('##INFO=<ID=Pvalue,Number=1,Type=Float,Description="Pvalue info">\n')

            # FORMAT Tanımları (Franklin bunu görmek ister)
            f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            f.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles">\n')
            f.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">\n')
            f.write('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">\n')
            f.write('##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes">\n')

            # Sütun Başlıkları
            f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tMySample\n')

        # Veriyi ekle (Quote None önemli)
        vcf_df.to_csv(output_filename, sep='\t', index=False, mode='a', header=False, quoting=csv.QUOTE_NONE, escapechar='\\')

        print(f"Tamamlandı! '{output_filename}' dosyasını Franklin'e yükleyebilirsin.")

    # Dosya Listesi
    files = [
        "plink/significant_pathways_snps_lof_caa.csv",
        "plink/significant_pathways_snps_lof_case_control.csv",
        "plink/significant_pathways_snps_lof_Consanguineous_marriage.csv",
        "plink/significant_pathways_snps_lof_famhis_chd.csv",
        "plink/significant_pathways_snps_lof_Family_History_Status.csv",
        "plink/significant_pathways_snps_mis_case_control.csv",
        "plink/significant_pathways_snps_mis_Consanguineous_marriage.csv",
        "plink/significant_pathways_snps_mis_famhis_chd.csv",
        "plink/significant_pathways_snps_mis_Family_History_Status.csv",
        "plink/significant_pathways_snps_mis_KD_in_siblings.csv",
        "plink/significant_pathways_snps_mis_Sequelae_Status.csv",
        "plink/significant_pathways_snps_mis_caa.csv"
    ]

    if __name__ == "__main__":
        create_franklin_ready_vcf(files)
    return


if __name__ == "__main__":
    app.run()
