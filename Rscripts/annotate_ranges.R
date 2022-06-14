library(annotatr)

annotate_ranges = function(ranges, annotations, output="both", tryUnique=TRUE){
  names(annotations) = NULL
  suppressWarnings({
    ranges.annot = ChIPpeakAnno::annotatePeakInBatch(ranges, AnnotationData=annotations, output = output)
  })
  ranges.annot$annot_index = as.numeric(str_remove_all(ranges.annot$feature, "ann"))
  ranges.annot$feature_type = annotations$type[ranges.annot$annot_index]
  ranges.annot$new_dist = ranges.annot$shortestDistance
  ranges.annot$new_dist[ranges.annot$fromOverlappingOrNearest == "Overlapping"] = 0
  ranges.annot$insideFeature = factor(ranges.annot$insideFeature, levels = c("includeFeature", "inside", "overlapStart", "overlapEnd", "upstream", "downstream"))
  names(ranges.annot) = NULL
  ranges.annot = arrange(as.data.frame(ranges.annot), peak, new_dist, insideFeature, shortestDistance)
  if(tryUnique){
    ranges.annot.singlepeak = ranges.annot[!duplicated(paste0(ranges.annot$peak, ranges.annot$feature_strand)),]
    annot_name = paste(ranges.annot$peak, ranges.annot$new_dist, ranges.annot$insideFeature, ranges.annot$feature_strand, sep="_")
    annot_name.singlepeak = paste(ranges.annot.singlepeak$peak, ranges.annot.singlepeak$new_dist, ranges.annot.singlepeak$insideFeature,  ranges.annot.singlepeak$feature_strand, sep="_")
    ranges.annot = ranges.annot[annot_name %in% annot_name.singlepeak,] 
  }
  return(ranges.annot)
}
