annots = builtin_annotations() %>% grep("mm10", ., value = T) %>% grep("cpg", ., value = T, invert = T)%>% grep("lnc", ., value = T, invert = T)%>% grep("bound", ., value = T, invert = T)%>% grep("intron", ., value = T, invert = T)
annotations = build_annotations(genome = 'mm10', annotations = annots)
annotations = annotations[annotations$type != "mm10_genes_introns"]

## LIFTOVER TO MM39
library(rtracklayer)
ch = import.chain("mm10ToMm39.over.chain")
seqlevelsStyle(annotations)
annotations_lifted = liftOver(annotations, ch)
annotations_lifted = unlist(annotations_lifted)
annotations = annotations_lifted

gencode_annots = fread("~/Downloads/gencode.vM29.annotation.gtf.gz", data.table = F)
colnames(gencode_annots) = c("chrom", "source", "type", "start", "end", "score", "strand", "frame", "attributes")
gencode.exons = gencode_annots[gencode_annots$type == "exon",]
gencode.exons$gene_id = strv_split2(gencode.exons$attributes, 'gene_id "', '";')
gencode.exons$gene_name = strv_split2(gencode.exons$attributes, 'gene_name "', '";')
gencode.exons$transcript_id = strv_split2(gencode.exons$attributes, 'transcript_id "', '";')
gencode.exons$exon_id = strv_split2(gencode.exons$attributes, 'exon_id "', '";')
gencode.exons$exon_number = as.numeric(strv_split2(gencode.exons$attributes, 'exon_number ', ';'))

n_cores=1
save(annotations, gencode.exons, n_cores, file="mouse_annotation_data.rdata")
