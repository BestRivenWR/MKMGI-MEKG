suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))
stopifnot(nzchar(Sys.getenv('MEKG_WORKSPACE')))
root <- file.path(Sys.getenv('MEKG_WORKSPACE'),'v31')
obj <- readRDS(file.path(root,'derived/GSE270866_merged_seurat.rds'))
x <- GetAssayData(obj,assay='RNA',slot='counts')
meta <- obj@meta.data[colnames(x),,drop=FALSE]
cat('RNA_COUNTS_DIM',dim(x),'\n')
stopifnot(ncol(x)==31506, nrow(x)>=28137, all(x@x>=0), all(x@x==floor(x@x)))
stopifnot(identical(rownames(meta),colnames(x)))
stopifnot(!anyNA(meta$annotation.predicted), !anyNA(meta$seurat_clusters))
q <- data.frame(metric=c('nuclei','native_RNA_features','RNA_counts_total',
  'nCount_RNA_total','cells_with_RNA_count_discrepancy','from_filtered','from_raw',
  'non_Singlet_DFPrediction_retained'),value=c(ncol(x),nrow(x),sum(x),sum(meta$nCount_RNA),
  sum(Matrix::colSums(x)!=meta$nCount_RNA),sum(meta$source_matrix=='filtered'),
  sum(meta$source_matrix=='raw'),sum(meta$DFPrediction!='Singlet')))
write.table(q,file.path(root,'derived/seurat_export_qc.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
group_stats <- list()
out <- gzfile(file.path(root,'derived/seurat_group_expression.tsv.gz'),'wt')
writeLines('native_gene_id\tgroup_id\tUMI_sum\tdetected_nuclei',out)
for(level in c('cell_type','cluster')) {
  lab <- if(level=='cell_type') as.character(meta$annotation.predicted) else as.character(meta$seurat_clusters)
  keys <- paste(meta$orig.ident,lab,sep='|')
  levels_order <- sort(unique(keys))
  m <- sparseMatrix(i=seq_along(keys),j=match(keys,levels_order),x=1,
    dims=c(ncol(x),length(levels_order)))
  sums <- as.matrix(x %*% m)
  detected <- as.matrix((x>0) %*% m)
  for(j in seq_along(levels_order)) {
    idx <- which(keys==levels_order[j]);gid <- paste(level,levels_order[j],sep='|')
    total <- sum(sums[,j]);stopifnot(total>0,length(idx)>0)
    st <- data.frame(group_id=gid,group_level=level,replicate=as.character(meta$orig.ident[idx[1]]),
      label=lab[idx[1]],nuclei=length(idx),UMI_total=total,
      from_raw=sum(meta$source_matrix[idx]=='raw'),from_filtered=sum(meta$source_matrix[idx]=='filtered'),
      non_Singlet_DFPrediction=sum(meta$DFPrediction[idx]!='Singlet'),
      mean_author_annotation_score=mean(meta$annotation.predicted.score[idx]),
      not_spatial_coordinates=TRUE,stringsAsFactors=FALSE)
    group_stats[[length(group_stats)+1]] <- st
    d <- data.frame(native_gene_id=rownames(x),group_id=gid,UMI_sum=sums[,j],detected_nuclei=detected[,j])
    write.table(d,out,sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)
  }
  cat('EXPORTED',level,length(levels_order),'groups\n')
}
close(out)
write.table(do.call(rbind,group_stats),file.path(root,'derived/seurat_group_stats.tsv'),
  sep='\t',quote=FALSE,row.names=FALSE)
emb <- Embeddings(obj,'SCT_umap')
write.table(cbind(barcode=rownames(emb),emb),gzfile(file.path(root,'derived/seurat_umap.tsv.gz')),
  sep='\t',quote=FALSE,row.names=FALSE)
capture.output(sessionInfo(),file=file.path(root,'derived/seurat_R_session.txt'))
print(q)
file.create(file.path(root,'derived/seurat_export_COMPLETE'))
