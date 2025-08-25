import marimo

__generated_with = "0.14.9"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    return


@app.cell
def _():
    import hail as hl
    return (hl,)


@app.cell
def _(hl):
    mt = hl.import_vcf('test2.vep.vcf.gz',reference_genome= "GRCh38", force_bgz=True, array_elements_required=False)


    return (mt,)


@app.cell
def _(mt):
    mt.info.CSQ.show(5)
    return


@app.cell
def _(hl, mt):
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


    mt2 = mt.annotate_rows(
        CSQ_struct = hl.struct(**{
            field: mt.info.CSQ[0].split('\|')[i] for i, field in enumerate(csq_fields)
        })
    )

    mt2.CSQ_struct.show(5)
    return (mt2,)


@app.cell
def _(hl, mt2):
    # 2️⃣ AF ve gnomAD alanlarını float’a çevir
    af_fields = [
        "AF","AFR_AF","AMR_AF","EAS_AF","EUR_AF","SAS_AF",
        "gnomADe_AF","gnomADe_AFR_AF","gnomADe_AMR_AF","gnomADe_ASJ_AF","gnomADe_EAS_AF",
        "gnomADe_FIN_AF","gnomADe_MID_AF","gnomADe_NFE_AF","gnomADe_REMAINING_AF","gnomADe_SAS_AF",
        "gnomADg_AF","gnomADg_AFR_AF","gnomADg_AMI_AF","gnomADg_AMR_AF","gnomADg_ASJ_AF",
        "gnomADg_EAS_AF","gnomADg_FIN_AF","gnomADg_MID_AF","gnomADg_NFE_AF","gnomADg_REMAINING_AF",
        "gnomADg_SAS_AF","MAX_AF"
    ]

    # her AF alanını float’a cast et
    mt3 = mt2.annotate_rows(
        CSQ_struct = mt2.CSQ_struct.annotate(**{
            field: hl.if_else(mt2.CSQ_struct[field] != '', hl.float(mt2.CSQ_struct[field]), hl.missing(hl.tfloat))
            for field in af_fields
        })
    )

    mt3.CSQ_struct.show(5)

    return


app._unparsable_cell(
    r"""
    mt2.CSQ_struct.
    """,
    name="_"
)


@app.cell
def _(mt2):
    mt2.describe(
    )
    return


@app.cell
def _(mt):
    mt_nogap = mt.filter_rows(mt.alleles[1] != '*')

    return (mt_nogap,)


@app.cell
def _(hl, mt_nogap):
    result = hl.vep(mt_nogap, 'vep-configuration.json', csq=True)
    return (result,)


@app.cell
def _(result):
    result.describe()
    return


@app.cell
def _(result):
    gnomad_af_table = result.select_rows(
        colocated_gnomad_afs = result.vep.colocated_variants.frequencies
    )
    gnomad_af_table.rows().show(1000, truncate=False)

    return


if __name__ == "__main__":
    app.run()
