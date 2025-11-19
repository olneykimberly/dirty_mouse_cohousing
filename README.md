# dirty_mouse_cohousing
8 wk old mice co-housed for 8 weeks with dirty mice from a local pet store. 
RNA extracted from brain samples. 
Groups are clean, bedding, or cohoused (n=4 per sex per group, n= 24 total).


| Group                     | Count   | Sex (F, M)|
| ------------------------- |:-------:|:-------:|
| Clean                     | 8    | (4,4)   |
| Bedding                   | 8    | (4,4)   |
| Cohoused (CH)             | 8    | (4,4)   |

This git repo contains scripts for the following:
-   Metadata analysis
-   Processing of bulk RNA-sequencing data
-   Generation of manuscript figures 
-   Generation of shiny app for exploration of the results, view app [here](https://fryerlab.shinyapps.io/dirty_mice_cohousing/)


*Characterize cell-type transcriptional alterations across neuropathologies:*


## Set up conda environment
This workflow uses conda. For information on how to install conda [here](https://docs.conda.io/projects/conda/en/latest/user-guide/index.html)

To create the conda environment:
```
conda env create -n dirty_mice --file dirty_mice.yml

# To activate this environment, use
#
#     $ conda activate dirty_mice
#
# To deactivate an active environment, use
#
#     $ conda deactivate dirty_mice
```

## Reference genome
Reference genome and annotation were downloaded prior to running snakemake. 
Ensembl refdata-gex-GRCm39-2024-A fasta

```
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCm39-2024-A.tar.gz
STAR --runMode genomeGenerate --runThreadN 8 --genomeDir refdata-gex-GRCm39-2024-A_STARv2.7.11_150sjdb  --genomeFastaFiles genome.fa --sjdbGTFfile genes.gtf --sjdbOverhang 150
```

## Snakemake for trimming reads and alignment to the reference genome
See the snakefile for specific details for each step. 
The config file was generated using the 01 and 02 scripts for obtaining sample information and creating the config file. 
```
Snakemake -s Snakefile 
```
Output includes STAR reads per gene, which can then be read into R for differential expression.

## Variance assessment and differential expression
```
R 04_differential_expression.Rmd
```
Output is differentially expressed genes for all pairwise comparisons.

## References
All packages used in this workflow are publicly available. If you use this workflow please cite the packages used. 
If you use the data in this workflow cite the following:
[et al. 2026]()

## Contacts

| Contact | Email |
| --- | --- |
| Donna Roscoe | droscoe@tgen.org |
| Kimberly Olney, PhD | kolney@tgen.org |
| John Fryer, PhD | jfryer@tgen.org |
