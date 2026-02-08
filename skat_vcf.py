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

        print("Scanning files...")
        for f in file_list:
            if not os.path.exists(f):
                continue
            try:
                df = pd.read_csv(f, sep='\t')
                df.columns = df.columns.str.strip()
                dfs.append(df)
            except Exception as e:
                print(f"Error: {f} - {e}")

        if not dfs:
            print("No data found to process.")
            return

        combined_df = pd.concat(dfs, ignore_index=True)

        combined_df['Chr'] = combined_df['Chr'].astype(str).apply(lambda x: x if str(x).startswith('chr') else f'chr{x}')
        combined_df['Pos'] = combined_df['Pos'].astype(int)
        combined_df['Ref'] = combined_df['Ref'].astype(str).str.upper()
        combined_df['Alt'] = combined_df['Alt'].astype(str).str.upper()

        if 'SNP_ID' in combined_df.columns:
            combined_df['ID'] = combined_df['SNP_ID'].fillna('.')
        else:
            combined_df['ID'] = '.'

        combined_df['QUAL'] = '.'
        combined_df['FILTER'] = 'PASS'
        combined_df['FORMAT'] = 'GT:AD:DP:GQ:PL'
        combined_df['MySample'] = '0/1:10,10:20:99:0,30,300'

        def clean_val(x):
            return str(x).replace(" ", "_").replace(";", "_").replace("=", "_")

        combined_df['Pathway'] = combined_df['Pathway'].apply(lambda x: clean_val(x) if pd.notna(x) else "")
        combined_df['SetID'] = combined_df['SetID'].apply(lambda x: clean_val(x) if pd.notna(x) else "")
        combined_df['adj.P.value'] = combined_df['adj.P.value'].apply(lambda x: str(x) if pd.notna(x) else "")

        group_cols = ['Chr', 'Pos', 'ID', 'Ref', 'Alt', 'QUAL', 'FILTER', 'FORMAT', 'MySample']
    
        vcf_grouped = combined_df.groupby(group_cols).agg({
            'Pathway': lambda x: ",".join(sorted(list(set(filter(None, x))))),
            'SetID': lambda x: ",".join(sorted(list(set(filter(None, x))))),
            'adj.P.value': lambda x: ",".join(list(filter(None, x)))
        }).reset_index()

        def finalize_info(row):
            parts = []
            if row['Pathway']:
                parts.append(f"Pathway={row['Pathway']}")
            if row['SetID']:
                parts.append(f"SetID={row['SetID']}")
            if row['adj.P.value']:
                parts.append(f"adj.P.value={row['adj.P.value']}")
            return ";".join(parts) if parts else "."

        vcf_grouped['INFO'] = vcf_grouped.apply(finalize_info, axis=1)

        final_cols = ['Chr', 'Pos', 'ID', 'Ref', 'Alt', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'MySample']
        vcf_df = vcf_grouped[final_cols].copy()
        vcf_df.columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'MySample']
        def get_chrom_order(chrom):
            c = chrom.lower().replace('chr', '')
            if c.isdigit(): return int(c)
            mapping = {'x': 98, 'y': 99, 'm': 100, 'mt': 100}
            return mapping.get(c, 101)

        vcf_df['chrom_order'] = vcf_df['#CHROM'].apply(get_chrom_order)
        vcf_df = vcf_df.sort_values(by=['chrom_order', 'POS']).drop(columns=['chrom_order'])

        print(f"Generating VCF: {output_filename}")
        with open(output_filename, 'w', newline='', encoding='utf-8') as f:
            f.write('##fileformat=VCFv4.2\n')
            f.write('##reference=GRCh38\n')
            f.write('##source=GenotypeMatrix_Aggregated_v2\n')
        
            for c in [str(i) for i in range(1, 23)] + ['X', 'Y', 'M', 'MT']:
                f.write(f'##contig=<ID=chr{c}>\n')

            # INFO Field Definitions
            f.write('##INFO=<ID=Pathway,Number=.,Type=String,Description="Associated pathways, comma separated">\n')
            f.write('##INFO=<ID=SetID,Number=.,Type=String,Description="Set IDs from original data, comma separated">\n')
            f.write('##INFO=<ID=adj.P.value,Number=.,Type=String,Description="Adjusted P-values, comma separated">\n')

            f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            f.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">\n')
            f.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">\n')
            f.write('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">\n')
            f.write('##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled likelihoods">\n')


            f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tMySample\n')

        vcf_df.to_csv(output_filename, sep='\t', index=False, mode='a', header=False, quoting=csv.QUOTE_NONE, escapechar='\\')
        print(f"Process completed A total of {len(vcf_df)} unique variants saved with SetID information.")

    if __name__ == "__main__":
        files = ["adj.skat.csv"]
        create_franklin_ready_vcf(files)
    return


if __name__ == "__main__":
    app.run()
