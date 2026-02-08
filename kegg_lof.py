import marimo

__generated_with = "0.17.7"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import statsmodels.api as sm
    import os 
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    import numpy as np
    import hail as hl

    hl.init()
    import requests
    import re
    return hl, re, requests


@app.cell
def _(hl):
    lof_filtered= hl.read_matrix_table('lof_variants.mt')
    return (lof_filtered,)


@app.cell
def _(lof_filtered):
    lof_filtered.describe()
    return


@app.cell
def _(re, requests):

    def get_pathway_gene_symbols(pathway_id):
        """
        Fetches gene symbols and the pathway name for a given pathway ID from the KEGG API.
        """
        url = f"http://rest.kegg.jp/get/{pathway_id}"
        try:
            resp = requests.get(url, timeout=10)
            resp.raise_for_status()  # To catch error codes (4xx, 5xx)
        except requests.exceptions.RequestException as e:
            print(f"Error occurred while fetching data for '{pathway_id}' from KEGG API: {e}")
            return None, None

        text = resp.text

        pathway_name = ""
        symbols = []
        gene_section = False

        for line in text.split("\n"):
            if line.startswith("NAME"):
                pathway_name = line[len("NAME"):].strip()
            elif line.startswith("GENE"):
                gene_section = True
                line = line[len("GENE"):].strip()
            elif gene_section and re.match(r"^\s+\d+", line):
                line = line.strip()
            else:
                if gene_section:
                    break
                else:
                    continue

            if gene_section and line:
                parts = line.split()
                if len(parts) >= 2:

                    gene_symbols_on_line = re.findall(r"\s(\S+?);", line)
                    symbols.extend(gene_symbols_on_line)


        symbols = list(set(symbols))
        return symbols, pathway_name

    def main():
        """
        Fetches gene symbols for the specified pathway IDs and prints the results.
        """

        pathway_ids = ["hsa04370", "hsa04060", "hsa04630","hsa04660", "hsa04350"]

        print(" Gene symbols and names of the updated pathways: \n")

        for p_id in pathway_ids:
            symbols, name = get_pathway_gene_symbols(p_id)
            if symbols and name:
                print(f" **{name} ({p_id})**")
                print("   Gene symbols:", ", ".join(symbols))
                print("   Total number of genes:", len(symbols))
                print("-" * 30)
            elif name is None:
                print(f" Data could not be retrieved from KEGG for ID {p_id}. Skipping.")
                print("-" * 30)

    if __name__ == "__main__":
        main()
    return


@app.cell
def _():
    vegf= ["PIK3R3", "PPP3R1", "PRKCA", "MAP2K2", "RAC1", "MAPK1", "RAC2", "MAPK13", "AKT2", "HRAS", "VEGFA", "NOS3", "PIK3R2", "HSPB1", "PIK3R1", "PLCG2", "SPHK2", "PTGS2", "PIK3CA", "SHC2", "P3R3URF", "CASP9", "MAPK11", "SPHK1", "MAP2K1", "PXN", "PLA2G4E", "PLA2G4F", "PTK2", "AKT1", "PPP3R2", "AKT3", "BAD", "PLCG1", "RAF1", "PRKCG", "MAPKAPK2", "CDC42", "NFATC2", "NRAS", "JMJD7", "PRKCB", "PPP3CB", "PLA2G4D", "PLA2G4B", "MAPK14", "MAPK12", "PPP3CC", "SH2D2A", "PIK3CB", "PIK3CD", "MAPKAPK3", "PLA2G4C", "KDR", "RAC3", "SRC", "PPP3CA", "PLA2G4A", "MAPK3", "KRAS", "RNASE1", "NANOS3", "PXDN","RNASE3", "NCOA3"]

    Cytokine= ["CCL20", "IL1RL1", "IL22", "GHR", "CSF1R", "ACVRL1", "IL23R", "CCL1", "BMP10", "PF4", "CCL28", "CCR9", "GDF10", "IL10RB", "IL4R", "INHBC", "IL18R1", "CCL26", "IFNAR2", "TNFRSF14", "IL6", "LTA", "CXCR4", "GDF2", "IL6ST", "TNFRSF17", "CXCR6", "IL13", "BMPR1B", "TGFB3", "CRLF2", "IFNA16", "TNFSF10", "TNFRSF10A", "IL27", "TNFRSF25", "AMHR2", "BMP8B", "IL13RA2", "OSM", "CNTFR", "CCL3L3", "TNFSF13B", "TNFRSF11B", "TNFRSF11A", "TNFRSF10C", "EPO", "IL24", "CCL8", "ACVR1C", "BMP5", "IFNA8", "IFNK", "IL1A", "IL25", "TNFRSF21", "BMP7", "RELT", "CCR5", "IFNA21", "ACVR2B", "IL31", "IL36G", "IFNW1", "CXCL13", "BMP2", "IFNA7", "IL1B", "IFNA6", "LTB", "GDF5", "IL2RG", "CXCL17", "IFNL3", "GDF9", "IL2", "TSLP", "CD4", "CXCR2", "INHBA", "IL36A", "TNFSF4", "CCL4", "CCL15", "INHBB", "CX3CL1", "IL7R", "CCL7", "TNFSF18", "IL1F10", "AMH", "CCL17", "IL18RAP", "GDF7", "CSH1", "XCL1", "IL34", "IL17RA", "MSTN", "TNFSF15", "FAS", "IFNL1", "IFNGR1", "NGF", "IL5", "RELL2", "IL11RA", "CXCL9", "CCL21", "IL1RL2", "NODAL", "TNFSF11", "IL17A", "CXCR5", "TNF", "CXCL6", "IL2RA", "IL5RA", "IL17RC", "LEPR", "TNFRSF19", "CCL16", "CSF2RA", "IL19", "ACKR4", "GDF11", "IL27RA", "IL10RA", "CTF1", "TNFRSF13B", "TNFRSF12A", "IL17RB", "IL33", "EPOR", "IL36B", "IL1RN", "IL11", "IFNE", "GDF1", "ACVR2A", "IFNA1", "IL1R1", "IL12RB1", "XCR1", "CCL11", "IL36RN", "IFNA2", "IL21", "TNFSF12", "PF4V1", "CXCR3", "TNFRSF13C", "CCL23", "IL32", "BMP3", "IL9R", "IL22RA1", "CLCF1", "IFNG", "GH2", "CSF3R", "IL15RA", "LEP", "CXCL3", "TNFSF13", "IL15", "CCL4L2", "CSHL1", "CXCR1", "INHBE", "BMPR2", "IL17F", "EDAR", "IFNLR1", "CXCL11", "IL1RAP", "IL18", "NGFR", "TNFRSF10B", "CCL5", "OSMR", "GDF15", "IFNL2", "BMP4", "LIFR", "CCL13", "BMP8A", "TNFRSF1B", "IL16", "TNFRSF10D", "IFNA17", "CCL18", "LTBR", "CCL3", "CXCL5", "IL9", "IL17C", "IL23A", "IL6R", "CNTF", "TGFBR2", "IL17RE", "IL3RA", "CXCL2", "IL20RA", "CCR1", "IL12A", "TNFSF8", "ACVR1B", "CSF3", "ACKR3", "IL10", "IFNA13", "CSF1", "CCL4L1", "TNFRSF1A", "CCR7", "PPBP", "TNFSF14", "CD40LG", "CCR6", "IL13RA1", "IFNA4", "BMPR1A", "GDF3", "CCR10", "CCL25", "CCL22", "IL12RB2", "CCR3", "IL26", "CCL2", "MPL", "TGFBR1", "CD40", "PRLR", "IFNA14", "CCR2", "CCR4", "BMP6", "IFNGR2", "CCL27", "BMP15", "IL17D", "CCL24", "IL4", "TGFB2", "CSH2", "GDF6", "TNFSF9", "IFNA10", "THPO", "CCR8", "TNFRSF6B", "ACVR1", "PRL", "CSF2RB", "IL2RB", "IL31RA", "CXCL10", "EDA2R", "CD27", "CCL14", "CXCL14", "IFNAR1", "CXCL16", "IFNB1", "TNFRSF8", "IL21R", "CSF2", "CD70", "IL12B", "RELL1", "TGFB1", "INHA", "FASLG", "XCL2", "GH1", "IL3", "CX3CR1", "IL20RB", "IL37", "IL20", "IFNA5", "CXCL8", "IL1R2", "TNFRSF18", "TNFRSF9", "IL17B", "CXCL1", "CCL19", "EBI3", "TNFRSF4", "LIF", "CXCL12", "CCL3L1", "IL7", "EDA","ACKR2","CCM2","BTF3P11","TIMP1","EPX","MYDGF","CNOT6"]

    jak_stat=["IL22", "GHR", "IL23R", "IL10RB", "IL4R", "IFNAR2", "IL6", "CDKN1A", "STAM2", "STAT4", "IL6ST", "IL13", "HRAS", "STAT1", "CRLF2", "IFNA16", "IL27", "SOCS5", "IL13RA2", "OSM", "CNTFR", "STAT6", "PIK3CD", "EPO", "IL24", "IFNA8", "IFNK", "STAT5B", "IFNA21", "IL31", "CCND3", "IFNW1", "PIAS3", "IFNA7", "IL22RA2", "IFNA6", "IL2RG", "IFNL3", "SOCS4", "IL2", "TSLP", "EP300", "PIK3CB", "CISH", "SOCS7", "JAK1", "JAK2", "IL7R", "CREBBP", "CSH1", "GFAP", "SOS2", "PIAS2", "IFNL1", "MCL1", "IFNGR1", "IL5", "IL11RA", "AKT2", "IL2RA", "IL5RA", "PIK3R1", "LEPR", "TYK2", "SOCS1", "CSF2RA", "IL19", "IL27RA", "IL10RA", "CTF1", "EPOR", "IL11", "IFNE", "JAK3", "IFNA1", "IL12RB1", "IFNA2", "IL21", "PDGFRA", "IL9R", "PIAS4", "IL22RA1", "SOCS3", "CLCF1", "IFNG", "BCL2L1", "GH2", "CSF3R", "IL15RA", "LEP", "PIK3R3", "IL15", "CSHL1", "IFNLR1", "OSMR", "SOCS6", "EGF", "IFNL2", "LIFR", "SOCS2", "PDGFA", "IFNA17", "PIM1", "IL9", "IL23A", "IL6R", "CNTF", "AOX1", "IL3RA", "PIK3R2", "IL20RA", "IL12A", "PDGFB", "CSF3", "IL10", "CCND1", "IFNA13", "IL13RA1", "IFNA4", "GRB2", "PIK3CA", "PTPN2", "STAM", "AKT3", "IL12RB2", "IL26", "MTOR", "MPL", "STAT2", "PDGFRB", "PRLR", "IFNA14", "STAT3", "PIAS1", "P3R3URF","CCND2", "AKT1", "IFNGR2", "STAT5A", "IL4", "CSH2", "IFNA10", "PTPN6", "BCL2", "THPO", "PRL", "CSF2RB", "IL2RB", "IL31RA", "FHL1", "IFNAR1", "IFNB1", "IL21R", "CSF2", "IL12B", "MYC", "GH1", "IL3", "SOS1", "IL20RB", "IL20", "IFNA5", "PTPN11", "EGFR", "IRF9", "LIF", "RAF1","IL7","CDPF1","IL17D","MYDGF","TIMP1","EPX","LONP1","CFH","RNASE3"]

    t_cell=["NFKBIB", "LAT", "PPP3CB", "HRAS", "PPP3CA", "PRKCQ", "RHOA", "PIK3CD", "PPP2R2B", "LCK", "PPP2R5B", "MAPK14", "PAK3", "ZAP70", "PDCD1", "PAK6", "NFATC2", "NRAS", "PPP2R5E", "IL2", "PPP3R1", "CD4", "MALT1", "PIK3CB", "PPP2R5A", "RASGRP1", "PPP2R2C", "MAP2K7", "GSK3B", "CBLB", "PPP2CB", "CDC42", "MAPK11", "CD3G", "MAPK1", "SOS2", "IL5", "KRAS", "CHUK", "PPP3R2", "TNF", "AKT2", "MAPK8", "PIK3R1", "CD28", "PAK2", "PTPRC", "CD8A", "MAP3K8", "VAV2", "NFKBIE", "PAK4", "FOS", "ITK", "MAP3K7", "VAV3", "MAPK12", "CD3E", "PPP2R5C", "FYN", "IFNG", "MAP2K2", "PIK3R3", "PPP2R1B", "NFATC3", "CD8B2", "MAPK10", "VAV1", "PDPK1", "CD8B", "PAK5", "JUN", "PPP2R3C", "MAP3K14", "PPP2R1A", "MAP2K1", "NFKBIA", "CARD11", "LCP2", "PIK3R2", "MAPK3", "NCK1", "IKBKB", "TEC", "IL10", "CD40LG", "PPP2R3B", "PPP2R5D", "PIK3CA", "GRB2", "PLCG1", "ICOS", "NCK2", "AKT3", "CD3D", "P3R3URF", "AKT1", "MAPK9", "PPP2R2D", "CTLA4", "PPP2R3A", "PAK1", "IKBKG", "IL4", "BCL10", "CD247", "PTPN6", "RELA", "DLG1", "GRAP2", "NFATC1", "CSF2", "NFKB1", "PPP3CC", "PPP2CA", "CDK4", "SOS1", "PTPN11", "PPP2R2A", "MAPK13", "BUB1B-PAK6", "RAF1","SPNS1", "ARHGEF7","MMAB","PKN1","DLGAP5","RNASE3"]

    tgf_beta=["NRROS", "ID3", "CHRD", "BAMBI", "RBL1", "NOG", "INHBC", "SMURF2", "SIN3A", "E2F4", "CUL1", "TGIF2", "SMAD1", "BMPR1B", "SMAD3", "TGFB3", "NBL1", "THBS1", "AMHR2", "BMP8B", "E2F5", "CDKN2B", "LRRC32", "RHOA", "ACVR2B", "BMP5", "BMP7", "ACVR1C", "BMP2", "GDF5", "HDAC1", "TGIF1", "SKI", "LEFTY1", "DCN", "EP300", "INHBA", "SKP1", "SMAD2", "INHBB", "ZFYVE9", "AMH", "ID4", "PPP2CB", "CREBBP", "GDF7", "LTBP1", "PITX2", "SP1", "MAPK1", "FBN1", "NODAL", "HFE", "TNF", "SMAD6", "EMP3", "RBX1", "GREM1", "ID1", "ACVR2A", "FST", "SMAD9", "RGMA", "HJV", "IFNG", "RGMB", "HAMP", "RPS6KB1", "PPP2R1B", "INHBE", "BMPR2", "RPS6KB2", "ZFYVE16", "ID2", "SMAD4", "SMURF1", "SKIL", "BMP4", "BMP8A", "PPP2R1A", "FMOD", "NEO1", "TGFBR2", "MAPK3", "GREM2", "ACVR1B", "BMPR1A", "TMEM53", "MICOS10-NBL1", "TGFBR1", "TFDP1", "SMAD5", "LEMD3", "BMP6", "NCOR1", "TFRC", "THSD4", "TGFB2", "GDF6", "ACVR1", "IGSF1", "TGFB1", "MYC", "PPP2CA", "SMAD7", "LEFTY2", "ROCK1", "HDAC2", "TF", "TFR2","GARS1","HHAT","F3"]

    lit_genes = [
        "ADM", "BLK", "CACNA1E", "CARD11", "CASP3", "CD36", "CD40", "CD40L", 
        "CXCL16", "CXCL8", "FCRLA", "FCGR2A", "FGA", "FGFR4", "FNDC1", "FOXN1", 
        "GYG1", "HIST2H2AC", "HLA-DRB1", "IL17F", "IL17RC", "IL1B", "IL1RN", 
        "IL31RA", "IL6ST", "ITPKC", "MMP17", "MMP25", "MMP8", "MMP9", "MYH11", 
        "MYH14", "ORAI1", "PGS1", "PTGER4", "RAB32", "RBP3", "SH3GLB1", 
        "SIGLEC10", "SLC8A1", "SMAD3", "SMAD9", "TLR2", "TLR4", "TNFSF13B", 
        "TNFSF14", "VEGFA", "VEGFB", "VPS9D1", "SPI1", "HCK", "LILRB2", 
        "CD40LG", "H2AC20", "MYH7B", "TPGS1", "E2F1"]
    return Cytokine, jak_stat, lit_genes, t_cell, tgf_beta, vegf


@app.cell
def _(Cytokine, jak_stat, lit_genes, t_cell, tgf_beta, vegf):
    kegg_genes = []
    for g in vegf + Cytokine + jak_stat + t_cell +  tgf_beta + lit_genes :
        kegg_genes.append(g)
    return (kegg_genes,)


@app.cell
def _(hl, kegg_genes, lof_filtered):
    kegg_genes_final = hl.literal(set(kegg_genes))


    kegg_genes_lof = lof_filtered.filter_rows(
        kegg_genes_final.contains(lof_filtered.gene_symbol)
    )

    print("Original variant number:", lof_filtered.count_rows())
    print("Total pathway variant number:", kegg_genes_lof.count_rows())

    kegg_genes_lof.write("kegg_lof_variants.mt", overwrite=True)
    return (kegg_genes_lof,)


@app.cell
def _(Cytokine, jak_stat, lit_genes, t_cell, tgf_beta, vegf):
    gene_lists = {
        "vegf": vegf,
        "Cytokine": Cytokine,
        "jak_stat": jak_stat,
        "t_cell": t_cell,
        "tgf_beta": tgf_beta,
        "lit_genes": lit_genes
    }

    for list_name, genes in gene_lists.items():
        print(list_name)
        print("*" * 100)
        for gene in genes:
            print(gene)
        print("#" * 100)
    return


@app.cell
def _(kegg_genes_lof):
    kegg_genes_lof.entries().show(100)
    return


@app.cell
def _(
    Cytokine,
    hl,
    jak_stat,
    kegg_genes_lof,
    lit_genes,
    t_cell,
    tgf_beta,
    vegf,
):
    def write_setid_locus(matrix_table, gene_sets_dict, output_file):
        """
        for each gene set:
        set name -> chr:pos:ref:alt
        """
        rows = matrix_table.rows()

        with open(output_file, "w") as f:
            for set_name, gene_list in gene_sets_dict.items():
                gene_set = set(gene_list)
                filtered = rows.filter(hl.literal(gene_set).contains(rows.gene_symbol))

                loci_tuples = filtered.aggregate(
                    hl.agg.collect((filtered.locus.contig, filtered.locus.position,
                                    filtered.alleles[0], filtered.alleles[1]))
                )

                for contig, pos, ref, alt in loci_tuples:
                    if None not in (contig, pos, ref, alt):
                        f.write(f"{set_name} {contig}:{pos}:{ref}:{alt}\n")


    gene_sets = {
        "VEGF": vegf,
        "TGF_BETA": tgf_beta,
        "Cytokine": Cytokine,
        "JAK/STAT": jak_stat,
        "T-Cell": t_cell,
        "Literature_Genes": lit_genes}
    write_setid_locus(kegg_genes_lof, gene_sets, "plink/kegg_genes_lof.SetID")
    return


@app.cell
def _(hl, kegg_genes_lof):
    #familyhistorychd
    hl.export_plink(
        dataset=kegg_genes_lof,
        output='plink/kegg_genes_lof_famhis_chd',          
        pheno=kegg_genes_lof.pheno.Family_history_of_CHD_status
    )

    #caa status 
    hl.export_plink(
        dataset=kegg_genes_lof,
        output='plink/kegg_genes_lof_caa',          
        pheno=kegg_genes_lof.pheno.CAA_status
    )

    #Sequelae_Status
    hl.export_plink(
        dataset=kegg_genes_lof,
        output='plink/kegg_genes_lof_Sequelae_Status',          
        pheno=kegg_genes_lof.pheno.Sequelae_Status
    )

    #Family_History_Status
    hl.export_plink(
        dataset=kegg_genes_lof,
        output='plink/kegg_genes_lof_Family_History_Status',          
        pheno=kegg_genes_lof.pheno.Family_History_Status
    )


    #Consanguineous_marriage_status
    hl.export_plink(
        dataset=kegg_genes_lof,
        output='plink/kegg_genes_lof_Consanguineous_marriage_status',           
        pheno=kegg_genes_lof.pheno.Consanguineous_marriage_status
    )

    #KD_in_siblings_status
    hl.export_plink(
        dataset=kegg_genes_lof,
        output='plink/kegg_genes_lof_KD_in_siblings_status',         
        pheno=kegg_genes_lof.pheno.KD_in_siblings_status
    )
    return


@app.cell
def _(hl, kegg_genes_lof):
    kegg_genes_lof_csq2 = kegg_genes_lof.annotate_cols(
        case_numeric = hl.case()
            .when(kegg_genes_lof.case_control_status == "case", 1)
            .when(kegg_genes_lof.case_control_status == "control", 0)
            .or_missing()  
    )

    hl.export_plink(
        dataset=kegg_genes_lof_csq2,
        output='plink/kegg_genes_lof_case_control',
        pheno=kegg_genes_lof_csq2.case_numeric
    )
    return


if __name__ == "__main__":
    app.run()
