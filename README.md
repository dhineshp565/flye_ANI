# flye_ANI
Pipeline for Whole genome Assembly using Dragonflye and ANI analysis using fastANI and ANIclustermap

Usage

```
nextflow main.nf --input samplelist.csv --outdir test14 -profile docker --gsize 2.0M  --reference /data/Dhinesh/flye_ANI/PDS23214_flye.fasta
```
```
options
--input     csv file with headers SampleName,SamplePath
--outdir    output directory
--profile   docker or conda
--gsize     genome size (2.0M for genome size 2.0 Mbases)
--reference path for reference fasta file (.fa,.fasta)
```