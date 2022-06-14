all_enhancer = NULL
for(cell_type in c(paste0("E0", c("07","09", 67:74, 81,82)))){
  core_marks_dense = fread(paste0("~/Downloads/roadmap_chromhmm/", cell_type, "_15_coreMarks_hg38lift_dense.bed.gz"), header=F, data.table = F)
  core_marks_dense = core_marks_dense[,c(1,2,3,4)]
  colnames(core_marks_dense) = c("seqnames", "start", "end", "type")
  core_marks_dense = core_marks_dense[core_marks_dense$type == "7_Enh",]
  core_marks_dense$cell_type = cell_type
  all_enhancer = rbind(all_enhancer, core_marks_dense)
}

enhancers.reduced = GenomicRanges::reduce(makeGRangesFromDataFrame(all_enhancer))
ol_to_orig = as.data.frame(findOverlaps(enhancers.reduced, makeGRangesFromDataFrame(all_enhancer)))
ol_to_orig.counts = table(ol_to_orig$queryHits) %>% as.data.frame()

enhancers.reduced.common = enhancers.reduced[as.numeric(as.character(ol_to_orig.counts$Var1[ol_to_orig.counts$Freq > 1]))]
enhancers.reduced.common = enhancers.reduced

enhancers.reduced.common$id = paste0("brain_enhancer", ":", c(1:length(enhancers.reduced.common)))
enhancers.reduced.common$tx_id = NA
enhancers.reduced.common$gene_id = NA
enhancers.reduced.common$symbol = NA
enhancers.reduced.common$type = "brain_enhancer"

annots = builtin_annotations() %>% grep("hg38", ., value = T) %>% grep("cpg", ., value = T, invert = T)%>% grep("lnc", ., value = T, invert = T)%>% grep("bound", ., value = T, invert = T)%>% grep("intron", ., value = T, invert = T)
annotations = build_annotations(genome = 'hg38', annotations = annots)
annotations = annotations[annotations$type != "hg38_genes_introns"]
annotations.atr.brainenh = c(annotations, enhancers.reduced.common)

gencode_annots = fread("~/Downloads/gencode.v39.annotation.gtf.gz", data.table = F)
colnames(gencode_annots) = c("chrom", "source", "type", "start", "end", "score", "strand", "frame", "attributes")
gencode.exons = gencode_annots[gencode_annots$type == "exon",]
gencode.exons$gene_id = strv_split2(gencode.exons$attributes, 'gene_id "', '";')
gencode.exons$gene_name = strv_split2(gencode.exons$attributes, 'gene_name "', '";')
gencode.exons$transcript_id = strv_split2(gencode.exons$attributes, 'transcript_id "', '";')
gencode.exons$exon_id = strv_split2(gencode.exons$attributes, 'exon_id "', '";')
gencode.exons$exon_number = as.numeric(strv_split2(gencode.exons$attributes, 'exon_number ', ';'))

save(annotations.atr.brainenh, gencode.exons, n_cores, annotate_ranges, file="annotation_data.rdata")