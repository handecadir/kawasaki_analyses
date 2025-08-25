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
    from pprint import pprint

    from hail.plot import show

    hl.plot.output_notebook()
    return


@app.cell
def _(hl):
    hl.utils.get_1kg('data/')

    return


@app.cell
def _(hl):
    hl.import_vcf('data/1kg.vcf.bgz').write('data/1kg.mt', overwrite=True)

    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
