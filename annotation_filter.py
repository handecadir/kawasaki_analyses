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
    return (hl,)


@app.cell
def _(hl):
    result = hl.read_matrix_table('kawasaki_filtered_outliers.mt')
    return (result,)


@app.cell
def _(result):
    result.describe()
    return


@app.cell
def _(hl, result):

    all_float_fields = [
        "AF","AFR_AF","AMR_AF","EAS_AF","EUR_AF","SAS_AF",
        "gnomADe_AF","gnomADe_AFR_AF","gnomADe_AMR_AF","gnomADe_ASJ_AF","gnomADe_EAS_AF",
        "gnomADe_FIN_AF","gnomADe_MID_AF","gnomADe_NFE_AF","gnomADe_REMAINING_AF","gnomADe_SAS_AF",
        "gnomADg_AF","gnomADg_AFR_AF","gnomADg_AMI_AF","gnomADg_AMR_AF","gnomADg_ASJ_AF",
        "gnomADg_EAS_AF","gnomADg_FIN_AF","gnomADg_MID_AF","gnomADg_NFE_AF","gnomADg_REMAINING_AF",
        "gnomADg_SAS_AF","MAX_AF", "am_pathogenicity", "SIFT", "PolyPhen"
    ]


    mt = result.annotate_rows(
        CSQ_struct = result.CSQ_struct.annotate(
            **{
                field: hl.if_else(
                    result.CSQ_struct[field] != '',  
                    hl.if_else(
                        # If the field is SIFT or PolyPhen, use the special method
                        (field == "SIFT") | (field == "PolyPhen"),
                        hl.if_else(
                            result.CSQ_struct[field].contains('('),
                            hl.float(result.CSQ_struct[field].split('\\(')[1].split('\\)')[0]),
                            hl.float(result.CSQ_struct[field])
                        ),
                    
                        hl.float(result.CSQ_struct[field])
                    ),
                    hl.missing(hl.tfloat)  
                )
                for field in all_float_fields
            }
        )
    )
    return (mt,)


@app.cell
def _(mt):
    mt.count_rows()
    return


@app.cell
def _(mt):
    mt.CSQ_struct.show(150)
    return


@app.cell
def _(hl, mt):
    csq_set = mt.aggregate_entries(hl.agg.collect_as_set(mt.CSQ_struct.Consequence))
    print(csq_set)
    return


@app.cell
def _(hl, mt):
    mt_filtered = mt.filter_rows(
        (hl.is_missing(hl.float(mt.CSQ_struct.gnomADg_AF))) |
        (hl.float(mt.CSQ_struct.gnomADg_AF) < 0.01)
    )
    return (mt_filtered,)


@app.cell
def _(mt_filtered):
    mt_filtered.rows().show(100)
    return


@app.cell
def _(hl, mt_filtered):
    lof_consequences = hl.set([
        'stop_gained', 
        'splice_donor_variant', 
        'frameshift_variant', 
        'splice_acceptor_variant'
    ])

    lof_mt = mt_filtered.filter_rows(
        (lof_consequences.contains(mt_filtered.CSQ_struct.Consequence)) &
        (mt_filtered.CSQ_struct.CANONICAL == "YES") &  
        (hl.is_missing(mt_filtered.CSQ_struct.gnomADg_AF) | (mt_filtered.CSQ_struct.gnomADg_AF < 0.01))
    )

    lof_mt_annotated = lof_mt.annotate_rows(
        selected_consequence = lof_mt.CSQ_struct.Consequence
    )


    lof_mt_annotated.selected_consequence.show(1000)
    return (lof_mt,)


@app.cell
def _(lof_mt):
    lof_consequences_final_mt = lof_mt.annotate_rows(
        gene_symbol = lof_mt.CSQ_struct.SYMBOL,
        consequence = lof_mt.CSQ_struct.Consequence,
        is_canonical = lof_mt.CSQ_struct.CANONICAL,
        gnomad_genome_af = lof_mt.CSQ_struct.gnomADg_AF,
        gnomad_exome_af  = lof_mt.CSQ_struct.gnomADe_AF,
        polyphen = lof_mt.CSQ_struct.PolyPhen,
        sift = lof_mt.CSQ_struct.SIFT
    )

    lof_consequences_final_mt.show(10)
    return (lof_consequences_final_mt,)


@app.cell
def _(lof_consequences_final_mt):
    lof_consequences_final_mt.count_rows()
    return


@app.cell
def _(lof_consequences_final_mt):
    lof_consequences_final_mt.rows().show(10)
    return


@app.cell
def _(hl, mt_filtered):
    allowed_consequences = hl.set(["missense_variant"])

    # 2. Missense variant filter

    missense_mt = mt_filtered.filter_rows(
        (allowed_consequences.contains(mt_filtered.CSQ_struct.Consequence)) &
        (mt_filtered.CSQ_struct.CANONICAL == "YES") &
        (mt_filtered.CSQ_struct.PolyPhen >= 0.95) &
        (mt_filtered.CSQ_struct.SIFT <= 0.05) &
        (hl.is_missing(mt_filtered.CSQ_struct.gnomADg_AF) | (mt_filtered.CSQ_struct.gnomADg_AF < 0.01))
    )


    missense_mt = missense_mt.annotate_rows(
        gene_symbol = missense_mt.CSQ_struct.SYMBOL,
        consequence = missense_mt.CSQ_struct.Consequence,
        is_canonical = missense_mt.CSQ_struct.CANONICAL,
        gnomad_genome_af = missense_mt.CSQ_struct.gnomADg_AF,
        gnomad_exome_af  = missense_mt.CSQ_struct.gnomADe_AF,
        polyphen = missense_mt.CSQ_struct.PolyPhen,
        sift = missense_mt.CSQ_struct.SIFT
    )


    missense_mt.show(10)
    return (missense_mt,)


@app.cell
def _(missense_mt):
    missense_mt.rows().show(10)
    return


@app.cell
def _(missense_mt):
    missense_mt.count_rows()
    return


@app.cell
def _(lof_consequences_final_mt, missense_mt):
    missense_mt.write("missense.mt", overwrite=True)
    lof_consequences_final_mt.write("lof_variants.mt", overwrite=True)
    return


@app.cell
def _(lof_consequences_final_mt, missense_mt):
    # Export aligned subset of columns while retaining full schemas in memory
    miss_rows = missense_mt.rows()
    miss_rows = miss_rows.select(
        miss_rows.gene_symbol,
        miss_rows.consequence,
        miss_rows.is_canonical,
        miss_rows.gnomad_genome_af,
        miss_rows.gnomad_exome_af,
        miss_rows.polyphen,
        miss_rows.sift
    )
    miss_rows.export('rare_filtered_variants.csv')

    lof_rows = lof_consequences_final_mt.rows()
    lof_rows = lof_rows.select(
        lof_rows.gene_symbol,
        lof_rows.consequence,
        lof_rows.is_canonical,
        lof_rows.gnomad_genome_af,
        lof_rows.gnomad_exome_af,
        lof_rows.polyphen,
        lof_rows.sift
    )
    lof_rows.export('rare_filtered_lof_variants.csv')
    return


if __name__ == "__main__":
    app.run()
