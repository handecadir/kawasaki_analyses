import marimo

__generated_with = "0.14.9"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo 
    return (mo,)


@app.cell
def _(mo):
    mo.md(
        r"""
    ```
    bcftools merge -m none n7_exome_calls.vcf.gz control_exome_calls.vcf.bgz exome_calls.vcf.bgz |
    bcftools norm -m-any -O  z -o merge.vcf.gz
    ```

    Lines   total/split/realigned/skipped:	714050/46419/0/0
    """
    )
    return


@app.cell
def _():
    import hail as hl

    hl.init()
    return (hl,)


@app.cell
def _(hl):
    from pprint import pprint

    from hail.plot import show

    hl.plot.output_notebook()
    return


@app.cell
def _(hl):
    mt = hl.import_vcf('merge.vep.vcf.gz',reference_genome= "GRCh38", force_bgz=True, array_elements_required=False)


    return (mt,)


@app.cell
def _(mt):
    mt.rows().select().show(5)
    return


@app.cell
def _(mt):
    mt.s.show(5)
    return


@app.cell
def _(mt):
    mt.entry.take(100)
    return


@app.cell
def _(mt):
    mt.GT.show()

    return


@app.cell
def _(mt):
    mt.AD.show()
    return


@app.cell
def _(mt):
    mt.describe()
    return


@app.cell
def _(mt):
    mt.GT.show(2000)
    return


@app.cell
def _(hl, mt):
    locus = hl.locus("chr1", 17358, reference_genome="GRCh38")
    mt_locus = mt.filter_rows(mt.locus == locus)
    mt_locus.entries().show(n=2000)   # ilk 20 satırı göster

    return


@app.cell
def _(mt):
    # MatrixTable'i kaydet
    mt.write("kawasaki.mt", overwrite=True)
    return


if __name__ == "__main__":
    app.run()
