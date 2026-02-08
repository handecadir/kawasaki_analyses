import marimo

__generated_with = "0.17.7"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import hail as hl

    hl.init()
    from bokeh.layouts import gridplot
    return gridplot, hl


@app.cell
def _(hl):
    mt = hl.read_matrix_table('kawasaki.mt')
    return (mt,)


@app.cell
def _(mt):
    mt.entries().show(10)
    return


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


        condition = (
            (mt_with_vaf.GT.is_hom_ref())  # keep 0/0 always
            | ( (mt_with_vaf.VAF > min_vaf) & (mt_with_vaf.DP > min_dp) )
        )

        filtered_mt = mt_with_vaf.filter_entries(condition)
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
def _(gridplot, hl):
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
def _(filtered_mt, hl, plot_variant_qc_metrics):
    plot_variant_qc_metrics(hl.variant_qc(filtered_mt))
    return


@app.cell
def _(filtered_mt):
    filtered_mt.count_rows()
    return


@app.cell
def _(filtered_mt, hl, plot_variant_qc_metrics, remove_uncalled_variants):
    final_mt= remove_uncalled_variants(filtered_mt)
    plot_variant_qc_metrics(hl.variant_qc(final_mt))
    return (final_mt,)


@app.cell
def _(final_mt):
    final_mt.count_rows()
    return


@app.cell
def _(final_mt):
    final_mt.describe()
    return


@app.cell
def _(final_mt, gridplot, hl):
    # 1. Calculate the Call Rate for Case and Control groups separately and add it to the variants (rows)
    mt_check = final_mt.annotate_rows(
        call_rate_case = hl.agg.filter(
            final_mt.case_control_status == 'case', 
            hl.agg.fraction(hl.is_defined(final_mt.GT))
        ),
        call_rate_control = hl.agg.filter(
            final_mt.case_control_status == 'control', 
            hl.agg.fraction(hl.is_defined(final_mt.GT))
        )
    )

    #2. Plot the graphs (Side-by-side comparison)
    p_case = hl.plot.histogram(mt_check.call_rate_case, 
                               title='Case Group Call Rate', 
                               legend='Case CR', 
                               range=(0,1)) 

    p_control = hl.plot.histogram(mt_check.call_rate_control, 
                                  title='Control Group Call Rate', 
                                  legend='Control CR', 
                                  range=(0,1))


    hl.plot.show(gridplot([[p_case, p_control]]))
    return


@app.cell
def _(gridplot, hl, mt):
    mt_raw_check = mt.annotate_rows(
        cr_case = hl.agg.filter(
            mt.case_control_status == 'case', 
            hl.agg.fraction(hl.is_defined(mt.GT))
        ),
        cr_control = hl.agg.filter(
            mt.case_control_status == 'control', 
            hl.agg.fraction(hl.is_defined(mt.GT))
        )
    )

    # Calculate the difference between two groups
    mt_raw_check = mt_raw_check.annotate_rows(
        diff = hl.abs(mt_raw_check.cr_case - mt_raw_check.cr_control)
    )

    # Filter out variants with high difference (e.g. > 0.2)
    problematic_variants = mt_raw_check.filter_rows(mt_raw_check.diff > 0.2)

    count = problematic_variants.count()
    print(f"Number of variants with more than 20% call difference between Case and Control: {count}")

    problematic_variants.rows().select('cr_case', 'cr_control', 'diff').show(10)

    # Plot Call Rate histograms for Case and Control groups
    p_case_raw = hl.plot.histogram(mt_raw_check.cr_case, 
                                   title='Case Group Call Rate (RAW Data)', 
                                   legend='Case CR', 
                                   range=(0,1)) 

    p_control_raw = hl.plot.histogram(mt_raw_check.cr_control, 
                                      title='Control Group Call Rate (RAW Data)', 
                                      legend='Control CR', 
                                      range=(0,1))

    hl.plot.show(gridplot([[p_case_raw, p_control_raw]]))
    return (count,)


@app.cell
def _(count, filtered_mt, gridplot, hl):
    mt_raw_check2 = filtered_mt.annotate_rows(
        cr_case2 = hl.agg.filter(
            filtered_mt.case_control_status == 'case',  
            hl.agg.fraction(hl.is_defined(filtered_mt.GT))
        ),
        cr_control2 = hl.agg.filter(
            filtered_mt.case_control_status == 'control', 
            hl.agg.fraction(hl.is_defined(filtered_mt.GT))
        )
    )

    mt_raw_check2 = mt_raw_check2.annotate_rows(
        diff2 = hl.abs(mt_raw_check2.cr_case2 - mt_raw_check2.cr_control2)
    )

    problematic_variants2 = mt_raw_check2.filter_rows(mt_raw_check2.diff2 > 0.2)

    count2 = problematic_variants2.count()
    print(f"Number of variants with more than 20% call difference between Case and Control: {count}")


    problematic_variants2.rows().select('cr_case2', 'cr_control2', 'diff2').show(10)

    p_case_raw2 = hl.plot.histogram(mt_raw_check2.cr_case2, 
                                   title='Case Group Call Rate (First Filtered Data)', 
                                   legend='Case CR', 
                                   range=(0,1)) 

    p_control_raw2 = hl.plot.histogram(mt_raw_check2.cr_control2, 
                                      title='Control Group Call Rate (First Filtered Data)', 
                                      legend='Control CR', 
                                      range=(0,1))

    hl.plot.show(gridplot([[p_case_raw2, p_control_raw2]]))
    return


@app.cell
def _(final_mt):
    final_mt.count_rows()
    return


@app.cell
def _(final_mt):
    final_mt.write("kawasaki_filtered.mt", overwrite=True)
    return


if __name__ == "__main__":
    app.run()
