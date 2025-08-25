# Kawasaki Analyses

We downloaded the hail image as our base image.
```
docker.io/hailgenetics/hail:0.2.135-py3.11
```
We later installed marimo and vep on the container.

Run the container
```
podman run -it --memory=22g -e JAVA_TOOL_OPTIONS="-Xmx18g" -v /home/handecadir/Desktop/PhD_folders/projeler/kawasaki/:/home/ubuntu/:z -w /home/ubuntu -p 2718:2718 localhost/hail_marimo_vep:latest bash
```

We have our cases in two VCF files, exome\_calls.vcf.bgz contains 85 cases. n7\_exome\_calls.vcf.gz contains 7 cases. control\_exome\_calls.vcf.gz contains 100 controls. We have merged them all together using bcftools. While merging we haven't merged any of the alleles using `-m` flag. We later split any merged alleles using bcftools norm.

```
bcftools merge -m none exome_calls.vcf.bgz n7_exome_calls.vcf.gz control_exome_calls.vcf.bgz |
    bcftools norm -m-any -O  z -o merge.vcf.gz
```

```
Lines   total/split/realigned/skipped:	714050/46419/0/0
```

```
apt install tabix
/ensembl-vep/vep -i merge.vcf.gz -o merge.vep.vcf.gz --force_overwrite --cache --offline --assembly GRCh38 --allele_number --no_stats --minimal --pick --vcf --plugin AlphaMissense,file=AlphaMissense_hg38.tsv.gz,transcript_match=1 --compress_output bgzip --everything
```

## Initial analysis
This file reads in the merged vcf files and writes them as kawasaki.mt in hail format.

 __file: initial\_exploration.py__

## Quality filtering
This file reads in the kawasaki.mt. We filter for VAF>0.15 and DP>15. Resulting variants are saved as kawasaki\_filtered.mt.

 __file: quality\_filtering.py__

## Sex estimation
This file reads kawasaki\_filtered.mt and meta.csv. Calculates the sex with hail `impute\_sex` function and compares them to given sex from the metadata.

 __file: sex\_check.py__

## Annotation
This file reads in the kawasaki\_filtered.mt and we are annotating our variants with using hail vep function. We are writing vep\_annotated.mt.
 __file: annotation.py__

## Consequence and frequency filtering
We are reading vep\_annotated.mt. Filtering for missinse and nonsense variants. In order to filter high impact missense variants we are using polyphen >= 0.95 and SIFT <= 0.5. For nonsense variants we are using stop\_gained, splice\_donor\_variant, frameshift\_variant, splice\_acceptor\_variant consequences. For both missense and nonsense variants we are filtering rare variants using gnomadg\_af < 0.01.

 __file: variant_filtering.py__


## Association

 __file: association.py__
