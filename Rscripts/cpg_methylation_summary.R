library(data.table)
library(tidyverse)
library(halpme)
library(methylKit)
library(BSgenome.Mmusculus.UCSC.mm39)
library(aaRon)
library(Repitools)
setwd("~/Projects/UTAS/mouse_heart_nome/")

source("~/Projects/UTAS/NOMe_analysis/Rscripts/QC_functions.R")
source("~/Projects/UTAS/NOMe_analysis/Rscripts/methylation_summary_functions.R")

sample_info = fread("sample_info_mouse.txt", data.table = F, header = T)
sample_info$heart_chamber = factor(sample_info$heart_chamber, levels = c("RA", "LA", "RV", "LV"))
sample_info$mouse = factor(sample_info$mouse, levels = c("m1to5", "m6to10", "m11to15"))

gene_annot = rtracklayer::import("~/Downloads/gencode.vM29.annotation.gtf.gz")
prom_locs = promoters(gene_annot[gene_annot$type == "transcript" & gene_annot$transcript_type == "protein_coding" & gene_annot$transcript_support_level %in% c(1,2)], upstream = 1000, downstream = 500)


nome_cpg = fread("biscuit/NOMe.dedup.hcg.bed.gz", data.table = F)
nome_cpg$seqnames[is.na(nome_cpg$seqnames)] = "lambda"


chrom_cov = nome_cpg  %>% dplyr::select(ends_with(".cov") | starts_with("seqnames")) %>% group_by(seqnames) %>% summarise_all(sum)
missing_cov = colnames(chrom_cov[,which(apply(chrom_cov[,-1],2, function(x) any(x==0)))+1])
if(length(missing_cov) > 0){
  message("missing all chrom data for: ", paste(missing_cov, collapse = "; "))
}

### check lambda conversion
nome_lambda = nome_cpg[nome_cpg$seqnames == "lambda" | nome_cpg$seqnames =="J02459.1",]
lambda_conversion = colSums(nome_lambda[,grep(".C", colnames(nome_lambda))]) /
  colSums(nome_lambda[,grep(".cov", colnames(nome_lambda))])


nome_cpg = nome_cpg[nome_cpg$seqnames != "lambda" & nome_cpg$seqnames !="J02459.1" ,]
cpg.grange = makeGRangesFromDataFrame(nome_cpg)
ol.prom = findOverlaps(prom_locs, cpg.grange)

#cpg_windows = total_meth_inoverlap(nome_cpg, prom_locs, min_window_cov = 10)
cpg_windows = convert_biscuit2methwindows(nome_cpg, window_width = 1000, window_spacing = 100, min_window_cov = 10)
#cpg_windows = windows_reduced
cpg_windows.LETTER = cpg_windows[(seqnames(cpg_windows) %in% c("chrX", "chrY", "chrM"))]
cpg_windows = cpg_windows[!(seqnames(cpg_windows) %in% c("chrX", "chrY", "chrM"))]
## Plot PCA




#ol = findOverlaps(cpg_windows, c(prom_locs, gene_annot[gene_annot$type == "transcript" & gene_annot$transcript_type == "protein_coding" & gene_annot$transcript_support_level %in% c(1,2)]))
#ol = findOverlaps(cpg_windows, prom_locs)


#meth_percent = as.data.frame(elementMetadata(cpg_windows[unique(ol@from)]))
meth_percent = as.data.frame(elementMetadata(cpg_windows))
colnames(meth_percent) = colnames(elementMetadata(cpg_windows))
meth_percent = meth_percent[,grep(".meth_percent", colnames(meth_percent))]

dim(meth_percent)

min = apply(meth_percent, 1 ,min)
max = apply(meth_percent, 1 ,max)

summary(min)
summary(max)


var = apply(meth_percent, 1 ,var)
length(var)
pca = prcomp(t(meth_percent[which(var >= sort(var, decreasing = T)[10000]),]))
#pca = prcomp(t(meth_percent))
#pca = prcomp(t(meth_percent[which(var >= sort(var, decreasing = T)[500] & (min<=0.3)),]))

PCA = as.data.frame(pca$x)
PCA$sample_id = gsub(".meth_percent", "", colnames(meth_percent))
PCA = left_join(PCA, sample_info, by="sample_id")
ggplot(PCA, aes(x=PC1, y=PC2, col=heart_chamber, shape = mouse,label=sample_id))  + geom_point(size=10) + theme_classic() +scale_color_brewer(palette = "Paired")
ggplot(PCA, aes(x=PC1, y=PC2, shape=heart_chamber, col = mouse,label=sample_id))  + geom_point(size=10) + theme_classic() +scale_shape_manual(values = c(15:18))


## Plot PCA
rot = as.data.frame(pca$rotation)
rot$rowN = rownames(rot)

rot = arrange(rot, abs(PC1), desc=T)
top = meth_percent[as.numeric(rot$rowN[1:100]),]
#top = gpc_windows_all[as.numeric(rot$rowN[1:100]),]

#colnames(top_data) = c("seqnames", "start", "end", "strand", paste0(c("cov_", "C_", "T_"), rep(top_samples, each=3)))
#methyl_percent = top_data[,-c(1:4)][,c(seq(2,(ncol(top_data)-4), by=3))] / top_data[,-c(1:4)][,c(seq(1,(ncol(top_data)-4), by=3))]
#colnames(methyl_percent) = top_samples
pheatmap::pheatmap(top)

pheatmap::pheatmap(top_data[,-c(1:4)][,c(seq(1,(ncol(top_data)-4), by=3))])



#####
gpc_files = list.files("biscuit/NOME_bychr/", full.names = T)

for(i in seq_along(gpc_files)){
  if(i == 1) first_window = T
  nome_gpc = fread(gpc_files[i], data.table = F)

  try({
  gpc_windows = convert_biscuit2methwindows(nome_gpc, window_width = 200, window_spacing = 20, min_window_cov = 1)

  if(first_window == T){
    gpc_windows_all = gpc_windows
    first_window = F
  }else{
    suppressWarnings({
      gpc_windows_all = c(gpc_windows_all, gpc_windows)
    })
  }
  })
}


## Plot PCA
gpc_windows_all[which(!is.na(meth_percent$NOMe.meth_percent))]

meth_percent = as.data.frame(elementMetadata(gpc_windows_all))
colnames(meth_percent) = colnames(elementMetadata(gpc_windows_all))
meth_percent = meth_percent[,apply(meth_percent, 2, function(x) length(which(is.na(x))) <= length(x)*0.1)]
var = apply(meth_percent, 1 ,var)
pca = prcomp(t(meth_percent[which(var >= sort(var, decreasing = T)[500]),]))
pca = prcomp(t(meth_percent))

screeplot(pca)
pca$sdev[1:2]


PCA = as.data.frame(pca$x)
PCA$sample_id = gsub(".meth_percent", "", colnames(elementMetadata(gpc_windows)))
PCA = left_join(PCA, sample_to_treat, by="sample_id")
ggplot(PCA, aes(x=PC1, y=PC2, col=ad_pathology, label=sample_id))  + geom_text()

rot = as.data.frame(pca$rotation)
rot$rowN = rownames(rot)

rot = arrange(rot, abs(PC2), desc=T)
top = meth_percent[as.numeric(rot$rowN[1:100]),]
top = gpc_windows_all[as.numeric(rot$rowN[1:100]),]

#colnames(top_data) = c("seqnames", "start", "end", "strand", paste0(c("cov_", "C_", "T_"), rep(top_samples, each=3)))
#methyl_percent = top_data[,-c(1:4)][,c(seq(2,(ncol(top_data)-4), by=3))] / top_data[,-c(1:4)][,c(seq(1,(ncol(top_data)-4), by=3))]
#colnames(methyl_percent) = top_samples
pheatmap::pheatmap(top)

pheatmap::pheatmap(top_data[,-c(1:4)][,c(seq(1,(ncol(top_data)-4), by=3))])


1


























cov_cols = which(str_sub(colnames(nome_cpg), -4,-1) == ".cov")
C_cols = which(str_sub(colnames(nome_cpg), -2,-1) == ".C")
high_cov = which(rowSums(nome_cpg[,cov_cols] > 1) == length(cov_cols))



meth_percent = nome_cpg[high_cov,C_cols] / nome_cpg[high_cov,cov_cols]

pca = prcomp(t(meth_percent))
PCA = as.data.frame(pca$x)
PCA$sample_id = gsub(".C", "", colnames(meth_percent))

PCA = left_join(PCA, sample_to_treat, by="sample_id")

ggplot(PCA, aes(x=PC1, y=PC2, col=ad_pathology, label=sample_id))  + geom_text()





out_folder = "biscuit/methylkit"

convert_biscuit2methylkit(nome_cpg, min_cov = 3, out_folder = "biscuit/methylkit")

directory = out_folder
sample_to_treat = sample_info[,c("sample_id", "ad_pathology")]
sample_to_treat$sample_id = gsub("_RPT", "", sample_to_treat$sample_id)

read_in_methylKitObj = function(directory, sample_to_treat = sample_info[,c("sample_id", "ad_pathology")]){

  meth_files = list.files(directory, full.names = T)

  file.list = list()
  for(filen in meth_files){
    file.list = append(file.list, filen)
  }

  sample.list = list()
  for(sam in basename(meth_files) %>% gsub("_methylkit.txt", "", .) %>% gsub(".gz", "",.)){
    sample.list = append(sample.list, sam)
  }

  treats = sample_to_treat[,2][match(unlist(sample.list), sample_to_treat[,1])]
  treats = as.factor(treats)
  #as.numeric(treats)-1

  meth_obj = methRead(file.list,
                      assembly = "hg38",
                      treatment = as.numeric(treats)-1,
                      sample.id = sample.list,
                      context = "CpG", mincov = 1)

  meth = unite(meth_obj)
  return(meth)

}
meth@treatment
meth@sample.ids

clusterSamples(meth, dist = "correlation", method="ward", plot = TRUE)
meth_data = as.data.frame(meth@.Data)

PCASamples(meth, screeplot = F, sd.filter = T, sd.threshold = 0.2)

pca = PCASamples(meth,obj.return=T, sd.filter = T, sd.threshold = 0.2)
PCA = as.data.frame(pca$x)
PCA$sample = rownames(PCA)

rot = as.data.frame(pca$rotation)
rot$rowN = rownames(rot)

rot = arrange(rot, abs(PC1), desc=T)
top = meth[as.numeric(rot$rowN[1:100]),]
top_samples = top@sample.ids %>% gsub("High_|Low_|Intermediate_", "", .)

top_data = as.data.frame(top@.Data)
colnames(top_data) = c("seqnames", "start", "end", "strand", paste0(c("cov_", "C_", "T_"), rep(top_samples, each=3)))
methyl_percent = top_data[,-c(1:4)][,c(seq(2,(ncol(top_data)-4), by=3))] / top_data[,-c(1:4)][,c(seq(1,(ncol(top_data)-4), by=3))]
colnames(methyl_percent) = top_samples
pheatmap::pheatmap(methyl_percent)

pheatmap::pheatmap(top_data[,-c(1:4)][,c(seq(1,(ncol(top_data)-4), by=3))])

PCA$treat = halpme::strv_split(PCA$sample, "_", 1)
PCA$sample_id = str_remove_all(PCA$sample, "High_") %>% str_remove_all("Intermediate_") %>% str_remove_all("Low_")

PCA = left_join(PCA, sample_to_treat, by="sample_id")

ggplot(PCA, aes(x=PC1, y=PC2, col=ad_pathology, label=sample_id))  + geom_text()
