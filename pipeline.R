library(genomicsutil)
library(ioutil)
library(miscutil)
library(SingleCellExperiment)
library(scater)
library(RColorBrewer)

count_fn = "/work-zfs/abattle4/ashis/progdata/scnet/ye_lab_data/umi_counts_by_cell_type/umi_counts_CD8_T_cells.feather"
metadata_fn = "/work-zfs/abattle4/ashis/progdata/scnet/ye_lab_data/umi_counts_by_cell_type/metadata_CD8_T_cells.txt_per_cell.txt"
# TODO: edit gene annotation to have expression gene ids as rownames in the annotation file
#gene_annot_fn = "/work-zfs/abattle4/lab_data/annotation/gencode.v26/gencode.v26.annotation.gene.txt" 
gene_annot_fn = "/work-zfs/abattle4/ashis/prog/scproc/results/gene_annot_CD8.txt"
out_pfx = "results/processed"

pca_n_features = 1000  # the number of most variable features to use for PCA.
pca_n_components = 2   # the number of principal components to obtain.
pca_outlier_n_mad = 5  # the number of mad to declare outlier in PCA

qc_outlier_n_mad = 3   # the number of mad to declare outlier in QC metrics

### auto output paths
qc_plot_fn = sprintf('%s_qc_plot.pdf', out_pfx)
outlier_plot_fn = sprintf('%s_outlier_plot.pdf', out_pfx)

### read count data
if(endsWith(x = count_fn, suffix = '.feather')){
  count_df = read_feather_df(fn = count_fn, rownames.col = 1)
} else {
  count_df = read_df(fn = count_fn)
}

sce <- SingleCellExperiment(list(counts=as.matrix(count_df)))

### read cell-level metadata and add to sce
metadata_df = read_df(metadata_fn, stringsAsFactors = T)
stopifnot(all(colnames(count_df) %in% rownames(metadata_df)))  # all cells must have metadata
metadata_df = metadata_df[colnames(count_df),,drop=F]          # identical order
for(metacol in colnames(metadata_df))
  colData(sce)[,metacol] = metadata_df[,metacol]


### add gene-level metadata to sce
gene_annot_df = read_df(gene_annot_fn) # rownames correspond to gene ids

# gene_annot_df0 = gene_annot_df
# gene_annot_df2 = gene_annot_df
# gene_annot_df2 = gene_annot_df2[!duplicated(gene_annot_df2$gene_name),]
# rownames(gene_annot_df2) = gene_annot_df2$gene_name
# genes = rownames(sce)
# matched_idx = match(genes, rownames(gene_annot_df2))
# unmatched_idx = setdiff(1:nrow(gene_annot_df2), matched_idx)
# matched_idx[is.na(matched_idx)] = unmatched_idx[1:sum(is.na(matched_idx))]
# gene_annot_df2 = gene_annot_df2[matched_idx, ]
# rownames(gene_annot_df2) = rownames(sce)
# gene_annot_df = gene_annot_df2

gene_annot_df = gene_annot_df[rownames(sce), ,drop = F]
stopifnot(all(rownames(sce) %in% rownames(gene_annot_df)))
stopifnot(all(c('gene_id', 'gene_name', 'chr', 'start_pos', 'end_pos', 'strand') %in% colnames(gene_annot_df)))
for(annotcol in colnames(gene_annot_df))
  rowData(sce)[,annotcol] = gene_annot_df[,annotcol]

rowData(sce)$chr = extend_chr(rowData(sce)$chr)
rowData(sce)$gene_length = rowData(sce)$end_pos - rowData(sce)$start_pos + 1
rowData(sce)$tss = mapply(function(st, en, strand) ifelse(strand == "+", st, en), rowData(sce)$start_pos, rowData(sce)$end_pos, rowData(sce)$strand)

### compute QC metrics
mito <- which(rowData(sce)$chr=="chrM")
sce <- calculateQCMetrics(sce, feature_controls=list(Mt=mito))

### plot QC
pdf(qc_plot_fn)

hist(sce$total_counts/1e3, xlab="Library sizes (thousands)", main="", breaks=20, col="grey80", ylab="Number of cells")
hist(sce$total_features_by_counts, xlab="Number of expressed genes", main="", breaks=20, col="grey80", ylab="Number of cells")
hist(sce$pct_counts_Mt, xlab="Mitochondrial proportion (%)", ylab="Number of cells", breaks=20, main="", col="grey80")

plotColData(sce, y="total_counts", x="batch_cov")
plotColData(sce, y="total_counts", x="disease_cov")
plotColData(sce, y="total_features_by_counts", x="batch_cov")
plotColData(sce, y="total_features_by_counts", x="disease_cov")
plotColData(sce, y="pct_counts_Mt", x="batch_cov")
plotColData(sce, y="pct_counts_Mt", x="disease_cov")

plot(sce$total_features_by_counts, sce$total_counts/1e6, xlab="Number of expressed genes", ylab="Library size (millions)")
plot(sce$total_features_by_counts, sce$total_counts/1e6, xlab="Number of expressed genes", ylab="Library size (millions)")
plot(sce$total_features_by_counts, sce$pct_counts_Mt, xlab="Number of expressed genes", ylab="Mitochondrial proportion (%)")

dev.off()

### drop cells based on outliers
low.libsize.drop <- isOutlier(sce$total_counts, nmads=qc_outlier_n_mad, type="lower", log=TRUE)
high.libsize.drop <- isOutlier(sce$total_counts, nmads=qc_outlier_n_mad, type="higher", log=TRUE)
low.feature.drop <- isOutlier(sce$total_features_by_counts, nmads=qc_outlier_n_mad, type="lower", log=TRUE)
high.feature.drop <- isOutlier(sce$total_features_by_counts, nmads=qc_outlier_n_mad, type="higher", log=TRUE)
high.mito.drop <- isOutlier(sce$pct_counts_Mt, nmads=qc_outlier_n_mad, type="higher")
qc.drop = low.libsize.drop | high.libsize.drop | low.feature.drop | high.feature.drop | high.mito.drop

sce <- runPCA(sce, use_coldata=TRUE, ncomponents = pca_n_components, ntop = min(c(pca_n_features, dim(sce))), detect_outliers=F, exprs_values = 'counts')
pcs <- reducedDim(sce)
pca_outliers <- apply(pcs, MARGIN = 2, FUN = function(vals){
  isOutlier(vals, nmads=pca_outlier_n_mad, type="both", log=F)
})
pca.drop <- rowSums(pca_outliers) > 0

colData(sce)$low.libsize.drop = low.libsize.drop
colData(sce)$high.libsize.drop = high.libsize.drop
colData(sce)$low.feature.drop = low.feature.drop
colData(sce)$high.feature.drop = high.feature.drop
colData(sce)$high.mito.drop = high.mito.drop
colData(sce)$qc.drop = qc.drop
colData(sce)$pca.drop = pca.drop

### plot PCs colored by outlier status and known covarites
pdf(outlier_plot_fn)
outlier_colors = c(rgb(0,1,0,0.5), rgb(1,0,0,0.5)) # non outlier first, the outlier
outlier_pch = c(20, 18) # non outlier first, the outlier
tmp <- lapply(c('low.libsize.drop', 'high.libsize.drop', 'low.feature.drop', 'high.feature.drop', 'high.mito.drop', 'qc.drop', 'pca.drop'), 
              FUN = function(metric){
                outlier_colors = c(rgb(0,1,0,0.5), rgb(1,0,0,0.5)) # non outlier first, the outlier
                outlier_pch = c(20, 18) # non outlier first, the outlier
                plot(pcs[,1], pcs[,2], 
                     xlab = 'PC1', ylab = 'PC2', 
                     col = outlier_colors[as.integer(colData(sce)[,metric]) + 1], 
                     pch = outlier_pch[as.integer(colData(sce)[,metric]) + 1],
                     main = metric)
                
                legend('topright', 
                       legend = c('Not outlier', 'Outlier'),
                       col = c(outlier_colors),
                       pch = c(outlier_pch),
                       bg = rgb(1,1,1,0.3))
                invisible(NULL)
              })

point_colors <- brewer.pal(3, "RdGy")
point_palette <- colorRampPalette(point_colors)
tmp <- lapply(c(colnames(metadata_df), c('total_counts', 'total_features_by_counts', 'pct_counts_Mt')),
              FUN = function(cur_cov){
                cov_values = colData(sce)[,cur_cov]
                if(is.numeric(cov_values)){
                  cov_order = findInterval(cov_values, sort(cov_values))
                  plot(x = pcs[,1], y = pcs[,2], xlab = 'PC1', ylab = 'PC2', col = point_palette(length(cov_values))[cov_order], pch = 20, main = cur_cov)
                  legend("topright", col=point_palette(2), pch=20, legend=range(cov_values))
                } else {
                  fac_values = factor(cov_values)
                  plot(x = pcs[,1], y = pcs[,2], xlab = 'PC1', ylab = 'PC2', col = fac_values, pch = 20, main = cur_cov)
                  legend("topright", col=1:length(levels(fac_values)), pch=20, legend=levels(fac_values))
                }
                
                invisible(NULL)
              })

dev.off()

keep <- !(low.libsize.drop | high.libsize.drop | low.feature.drop | high.feature.drop | high.mito.drop | pca.drop)
drop_count = c(low.libsize.drop=sum(low.libsize.drop), 
               high.libsize.drop=sum(high.libsize.drop),
               low.feature.drop=sum(low.feature.drop),
               high.feature.drop=sum(high.feature.drop),
               high.mito.drop=sum(high.mito.drop),
               pca.drop = sum(pca.drop),
               Remaining=sum(keep))
verbose_print('Number of dropped cells')
print(drop_count)

### normalize data

### remove doublets (using doubletCells())

### cell cylcle identification and plot
# http://master.bioconductor.org/packages/release/workflows/vignettes/simpleSingleCell/inst/doc/umis.html#4_cell_cycle_classification

### Correct cell cycle, known batches, covariates and top PCs (?)

### examine most expressed genes (scater: plotHighestExprs)

### examine sparsity (scater:plotExprsFreqVsMean(sce))

### plotExplanatoryVariables

### find highly variant genes

