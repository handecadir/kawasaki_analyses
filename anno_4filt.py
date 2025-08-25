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
    return (mo,)


@app.cell
def _():
    import hail as hl

    hl.init()
    return (hl,)


@app.cell
def _(hl):
    #resultu kaydettik ya direkt açmak için bunu kullan artık 
    result = hl.read_matrix_table('kawasaki_filtered.mt')
    return (result,)


@app.cell
def _(result):
    result.describe()
    return


@app.cell
def _(hl, result):
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

    result2 = result.annotate_rows(
        CSQ_struct = hl.struct(**{
            field: result.info.CSQ[0].split('\|')[i] for i, field in enumerate(csq_fields)
        })
    )

    result2.CSQ_struct.show(5)

    return (result2,)


@app.cell
def _(hl, result2):
    # 1️⃣ Collect all fields to be converted in a single list
    all_float_fields = [
        "AF","AFR_AF","AMR_AF","EAS_AF","EUR_AF","SAS_AF",
        "gnomADe_AF","gnomADe_AFR_AF","gnomADe_AMR_AF","gnomADe_ASJ_AF","gnomADe_EAS_AF",
        "gnomADe_FIN_AF","gnomADe_MID_AF","gnomADe_NFE_AF","gnomADe_REMAINING_AF","gnomADe_SAS_AF",
        "gnomADg_AF","gnomADg_AFR_AF","gnomADg_AMI_AF","gnomADg_AMR_AF","gnomADg_ASJ_AF",
        "gnomADg_EAS_AF","gnomADg_FIN_AF","gnomADg_MID_AF","gnomADg_NFE_AF","gnomADg_REMAINING_AF",
        "gnomADg_SAS_AF","MAX_AF", "am_pathogenicity", "SIFT", "PolyPhen"
    ]


    mt = result2.annotate_rows(
        CSQ_struct = result2.CSQ_struct.annotate(
            **{
                field: hl.if_else(
                    result2.CSQ_struct[field] != '',  # Proceed if the value is not empty
                    hl.if_else(
                        # If the field is SIFT or PolyPhen, use the special method
                        (field == "SIFT") | (field == "PolyPhen"),
                        hl.if_else(
                            result2.CSQ_struct[field].contains('('),
                            hl.float(result2.CSQ_struct[field].split('\\(')[1].split('\\)')[0]),
                            hl.float(result2.CSQ_struct[field])
                        ),
                        # Otherwise (for AF/gnomAD fields), convert directly to float
                        hl.float(result2.CSQ_struct[field])
                    ),
                    hl.missing(hl.tfloat)  # If the value is empty, set as missing
                )
                for field in all_float_fields
            }
        )
    )

    return (mt,)


@app.cell
def _(mo):
    mo.md(
        r"""
    # 2️⃣ AF ve gnomAD alanlarını float’a çevir
    af_fields = [
        "AF","AFR_AF","AMR_AF","EAS_AF","EUR_AF","SAS_AF",
        "gnomADe_AF","gnomADe_AFR_AF","gnomADe_AMR_AF","gnomADe_ASJ_AF","gnomADe_EAS_AF",
        "gnomADe_FIN_AF","gnomADe_MID_AF","gnomADe_NFE_AF","gnomADe_REMAINING_AF","gnomADe_SAS_AF",
        "gnomADg_AF","gnomADg_AFR_AF","gnomADg_AMI_AF","gnomADg_AMR_AF","gnomADg_ASJ_AF",
        "gnomADg_EAS_AF","gnomADg_FIN_AF","gnomADg_MID_AF","gnomADg_NFE_AF","gnomADg_REMAINING_AF",
        "gnomADg_SAS_AF","MAX_AF","am_pathogenicity,"SIFT","PolyPhen"
    ]


    # AF alanlarını float'a çevir
    mt = result2.annotate_rows(
        CSQ_struct = result2.CSQ_struct.annotate(
            **{ 
                field: hl.if_else(
                    result2.CSQ_struct[field] != '', 
                    hl.float(result2.CSQ_struct[field]), 
                    hl.missing(hl.tfloat)
                ) 
                for field in af_fields
            }
        )
    )
    """
    )
    return


@app.cell
def _(mo):
    mo.md(
        r"""
    poly_sift_fields = ["SIFT", "PolyPhen"]

    # Sadece PolyPhen ve SIFT için özel bir annotasyon yapalım
    mt_poly = result2.annotate_rows(
        CSQ_struct = result2.CSQ_struct.annotate(
            **{
                field: hl.if_else(
                    result2.CSQ_struct[field] != '',
                    hl.if_else(
                        # İlk parantezin açılıp açılmadığını kontrol et
                        result2.CSQ_struct[field].contains('('),
                        # Eğer parantez varsa, içindeki değeri al
                        hl.float(result2.CSQ_struct[field].split('\\(')[1].split('\\)')[0]),
                        # Yoksa, doğrudan float'a çevir
                        hl.float(result2.CSQ_struct[field])
                    ),
                    # Eğer değer boşsa, missing (eksik) olarak ayarla
                    hl.missing(hl.tfloat)
                )
                for field in poly_sift_fields
            }
        )
    )

    mt_poly.CSQ_struct.show(750)
    """
    )
    return


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
        consequence_terms = lof_mt.CSQ_struct.Consequence,
        polyphen_score = lof_mt.CSQ_struct.PolyPhen,   
        sift_score = lof_mt.CSQ_struct.SIFT,           
        is_canonical = lof_mt.CSQ_struct.CANONICAL,    
        gnomad_genome_AF = lof_mt.CSQ_struct.gnomADg_AF,
        gnomad_exome_AF  = lof_mt.CSQ_struct.gnomADe_AF
    )

    lof_consequences_final_mt = lof_consequences_final_mt.select_rows(
        lof_consequences_final_mt.gene_symbol,
        lof_consequences_final_mt.consequence_terms,
        lof_consequences_final_mt.is_canonical,
        lof_consequences_final_mt.gnomad_genome_AF,
        lof_consequences_final_mt.gnomad_exome_AF
    )

    lof_consequences_final_mt.show(10)


    return (lof_consequences_final_mt,)


@app.cell
def _(lof_consequences_final_mt):
    lof_consequences_final_mt.rows().show(1000)
    return


@app.cell
def _(hl, mt_filtered):

    # (missense) variant
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
        polyphen = missense_mt.CSQ_struct.PolyPhen,
        sift = missense_mt.CSQ_struct.SIFT,
        gnomadg_af = missense_mt.CSQ_struct.gnomADg_AF,
        gnomade_af = missense_mt.CSQ_struct.gnomADe_AF
    )


    missense_mt.show(10)
    return (missense_mt,)


@app.cell
def _(missense_mt):
    missense_mt.rows().show(200)
    return


@app.cell
def _(lof_consequences_final_mt, missense_mt):

    missense_mt.write("missense.mt", overwrite=True)
    lof_consequences_final_mt.write("lof_variants.mt", overwrite=True)

    return


@app.cell
def _(lof_consequences_final_mt, missense_mt):
    missense_mt.rows().export('nadir_filtered_variants.csv')
    lof_consequences_final_mt.rows().export('nadir_filtered_lof_variants.csv')
    return


if __name__ == "__main__":
    app.run()
