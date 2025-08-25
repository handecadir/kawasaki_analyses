import marimo

__generated_with = "0.14.9"
app = marimo.App(width="medium")


@app.cell
def _():
    import hail as hl

    hl.init()
    return (hl,)


@app.cell
def _(hl):
    hl.import_vcf('exome_calls.vcf.gz').write('exome_calls.mt', overwrite=True)

    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
