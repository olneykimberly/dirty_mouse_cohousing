
pathToRef <- "/tgen_labs/jfryer/projects/references/mouse/"
# Read in annotation file
gtf.file <- paste0("/tgen_labs/jfryer/projects/references/mouse/refdata-gex-GRCm39-2024-A/genes/genes.gtf") #  refdata-gex-GRCm39-2024-A/genes/genes.gtf
gtf.gr <- rtracklayer::import(gtf.file)
# save gtf as data frame
gtf.df <- as.data.frame(gtf.gr)
# get gene id, transcript id, gene name, seqname which is chromosome, and biotype from gtf
genes <-
  gtf.df[, c("seqnames",
             "width",
             "gene_id",
             "gene_name",
             "gene_type", 
             "type")]
# Up date naming these columns using the correct column information.
names(genes)[names(genes) == "seqnames"] <- "Chr"
names(genes)[names(genes) == "lenght"] <- "width"
# keep gene_id to merge with counts data
genes$GENEID <- genes$gene_id
genes <- subset(genes, type == "gene")
#protein_coding_genes <- subset(genes, genes$gene_biotype == "protein_coding")

#protein_coding_genes <- protein_coding_genes %>% mutate(gene_name = coalesce(gene_name,gene_id))
# saveRDS(protein_coding_genes$gene_name, file = paste0("../rObjects/gene_options.rds"))
write.table(
  genes,
  "/tgen_labs/jfryer/kolney/dirty_mice/dirty_mouse_cohousing/scripts/ensembl_mouse_genes.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
# make tx2gene
# txdb <-
#   makeTxDbFromGFF(paste0(pathToRef, "Sus_scrofa.Sscrofa11.1.107.gtf"), format = "gtf")
# txdb_keys <- keys(txdb, keytype = "TXNAME")
# keytypes(txdb) # list of the different key types
# tx2gene <-
#   AnnotationDbi::select(txdb, txdb_keys, "GENEID", "TXNAME")
# 
# # clean up
# rm(txdb, pathToRef, gtf.file, gtf.gr, gtf.df, protein_coding_genes, txdb_keys)