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

    hl.init()
    return (hl,)


@app.cell
def _(hl):
    mt = hl.read_matrix_table('kawasaki.mt')
    return (mt,)


@app.cell
def _(hl):
    def filter_genotypes(
        mt: hl.MatrixTable,
        min_vaf: float = 0.15,
        min_dp: int = 15
    ) -> hl.MatrixTable:
        """Filters genotypes in a MatrixTable based on VAF and read depth.

        This function calculates the Variant Allele Frequency (VAF) for each
        genotype and then sets any genotype not meeting the min_vaf and
        min_dp thresholds to missing (NA).

        Args:
            mt: The input Hail MatrixTable.
            min_vaf: The minimum VAF required to keep a genotype.
            min_dp: The minimum read depth (DP) required to keep a genotype.

        Returns:
            A new MatrixTable with low-quality genotypes set to missing.
        """
        # Calculate VAF per genotype, handling cases where DP might be 0.
        mt_with_vaf = mt.annotate_entries(
            VAF=hl.if_else(mt.DP > 0, mt.AD[1] / mt.DP, 0.0)
        )

        # Filter entries, turning genotypes that fail the check into missing values.
        filtered_mt = mt_with_vaf.filter_entries(
            (mt_with_vaf.VAF > min_vaf) & (mt_with_vaf.DP > min_dp)
        )
        return filtered_mt




    def remove_uncalled_variants(mt: hl.MatrixTable) -> hl.MatrixTable:
        """Removes variants that have no called genotypes.

        This should be run after filtering entries to clean up the dataset
        by removing variants that are no longer polymorphic in any sample.

        Args:
            mt: The MatrixTable to clean (e.g., after running filter_genotypes).

        Returns:
            A new MatrixTable with uncalled variants removed.
        """
        # Run variant QC to get stats like the number of called genotypes.
        mt_with_qc = hl.variant_qc(mt)

        # Filter rows (variants) where the number of called genotypes is greater than 0.
        final_mt = mt_with_qc.filter_rows(mt_with_qc.variant_qc.n_called > 0)

        return final_mt
    return filter_genotypes, remove_uncalled_variants


@app.cell
def _(filter_genotypes, mt):
    filtered_mt = filter_genotypes(mt)
    return (filtered_mt,)


@app.cell
def _(filtered_mt, hl):
    hl.variant_qc(filtered_mt).row.describe()
    return


@app.cell
def _(hl):
    from bokeh.layouts import gridplot

    def plot_variant_qc_metrics(mt_with_qc: hl.MatrixTable):
        """
        Generates and displays a grid of plots for key variant QC metrics
        using the exact schema provided.

        Args:
            mt_with_qc: A MatrixTable with a `variant_qc` row struct.
        """
        # 1. Create individual histogram plots for each metric
        p1 = hl.plot.histogram(mt_with_qc.variant_qc.call_rate,
                               legend='Call Rate',
                               title='Variant Call Rate')

        p2 = hl.plot.histogram(mt_with_qc.variant_qc.AF[1],
                               legend='Allele Frequency',
                               title='Allele Frequency Spectrum')

        # Corrected path for mean DP
        p3 = hl.plot.histogram(hl.log10(mt_with_qc.variant_qc.dp_stats.mean),
                               legend='Mean Depth (log10)',
                               title='Mean Read Depth')

        # Corrected path for mean GQ
        p4 = hl.plot.histogram(mt_with_qc.variant_qc.gq_stats.mean,
                               legend='Mean GQ',
                               title='Mean Genotype Quality')

        # 2. Arrange the plots in a grid and show them
        plot_grid = gridplot([[p1, p2], [p3, p4]])
        hl.plot.show(plot_grid)
    return (plot_variant_qc_metrics,)


@app.cell
def _(hl, mt, plot_variant_qc_metrics):
    plot_variant_qc_metrics(hl.variant_qc(mt))
    return


@app.cell
def _(hl, mt, plot_variant_qc_metrics):
    plot_variant_qc_metrics(hl.variant_qc(mt))
    return


@app.cell
def _(filtered_mt, hl, plot_variant_qc_metrics):
    plot_variant_qc_metrics(hl.variant_qc(filtered_mt))
    return


@app.cell
def _(filtered_mt):
    filtered_mt.count_rows()
    return


@app.cell
def _(filtered_mt, remove_uncalled_variants):
    final_mt= remove_uncalled_variants(filtered_mt)
    return (final_mt,)


@app.cell
def _(final_mt):
    final_mt.count_rows()
    return


@app.cell
def _(final_mt):
    final_mt.write("kawasaki_filtered.mt")
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
