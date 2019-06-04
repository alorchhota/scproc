library(genomicsutil)
library(ioutil)
library(miscutil)
library(SingleCellExperiment)
library(scater)
library(RColorBrewer)
library(scran)
library(future)
library(sctransform)
library(umap)
library(dplyr)
library(ggplot2)
library(limma)
library(argparser)

args <- arg_parser("program");
args <- add_argument(args, "-count", help="umi count data", default="/work-zfs/abattle4/ashis/progdata/scnet/ye_lab_data/umi_counts_by_cell_type/umi_counts_Dendritic_cells.feather")
args <- add_argument(args, "-meta", help="cell meata data or covariate data", default="/work-zfs/abattle4/ashis/progdata/scnet/ye_lab_data/umi_counts_by_cell_type/metadata_Dendritic_cells.txt_per_cell.txt")
args <- add_argument(args, "-annot", help="gene annotation file", default="/work-zfs/abattle4/ashis/progdata/scproc/ye_annot_hg19/gene_annot_hg19_ye.txt")
args <- add_argument(args, "-core", help="number of cores", default=5)
args <- add_argument(args, "-o", help="output prefix", default="results/dendrytic/dendrytic")

argv = parse_args(args)
count_fn = argv$count
metadata_fn = argv$meta
gene_annot_fn = argv$annot
n_cores = argv$core
out_pfx = argv$o

# count_fn = "/work-zfs/abattle4/ashis/progdata/scnet/ye_lab_data/umi_counts_by_cell_type/umi_counts_CD8_T_cells.feather"
# metadata_fn = "/work-zfs/abattle4/ashis/progdata/scnet/ye_lab_data/umi_counts_by_cell_type/metadata_CD8_T_cells.txt_per_cell.txt"
# # TODO: edit gene annotation to have expression gene ids as rownames in the annotation file
# # gene_annot_fn = "/work-zfs/abattle4/lab_data/annotation/gencode.v26/gencode.v26.annotation.gene.txt" 
# # gene_annot_fn = "/work-zfs/abattle4/lab_data/annotation/gencode.v19/gencode.v19.annotation.gene.txt"
# gene_annot_fn = "/work-zfs/abattle4/ashis/progdata/scproc/ye_annot_hg19/gene_annot_hg19_ye.txt"
# out_pfx = "results/cd8_processed"

# count_fn = "/work-zfs/abattle4/ashis/progdata/scnet/ye_lab_data/umi_counts_by_cell_type/umi_counts_Dendritic_cells.feather"
# metadata_fn = "/work-zfs/abattle4/ashis/progdata/scnet/ye_lab_data/umi_counts_by_cell_type/metadata_Dendritic_cells.txt_per_cell.txt"
# gene_annot_fn = "/work-zfs/abattle4/ashis/progdata/scproc/ye_annot_hg19/gene_annot_hg19_ye.txt"
# out_pfx = "results/dendrytic/dendrytic_"


pca_n_features = 2000  # the number of most variable features to use for PCA.
pca_n_components = 4   # the number of principal components to obtain.
pca_outlier_n_mad = 5  # the number of mad to declare outlier in PCA
umap_n_reduced_dim = 50


qc_outlier_n_mad = 5   # the number of mad to declare outlier in QC metrics

individual_cov = "ind_cov"

max_mem = 30 * 1024 ^ 3 # 30GB

min_cell_per_gene = 10
min_individual_per_gene = 5
min_cells_for_itra_individual_pca_outlier = 20


vst_n_genes = 2000 # Number of genes to use when estimating parameters, default 2000
vst_batch_var = "batch_cov"
vst_latent_var = c('log10_total_counts')
#vst_latent_var_nonreg = NULL # c('disease_cov', 'pop_cov', 'well', 'total_features_by_counts', 'pct_counts_Mt', 'phase')
#vst_latent_var_nonreg = c("G1_score", "G2M_score", 'total_features_by_counts')
#vst_latent_var_nonreg = c("G1_score", "G2M_score", 'well')
vst_latent_var_nonreg = c("G1_score", "G2M_score")
pseudobulk_vst_latent_var_nonreg = NULL
vst_n_cells = NULL # Number of cells to use when estimating VST parameters, NULL means all.

max_doublet_score = 1.5
doublet_irlba_nv = 20
n_pcs_to_remove = 5

n_top_variable_genes = 500
hvg_min_count = 2
hvg_min_sample_per_gene = 5  # each gene must have count >= hvg_min_count in >= hvg_min_sample_per_gene samples.
hvg_min_total_variance_quantile = 0.1
fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))

### save settings
settings_df = t(data.frame(count_fn = count_fn,
                           metadata_fn = metadata_fn,
                           gene_annot_fn = gene_annot_fn,
                           out_pfx = out_pfx,
                           pca_n_features = pca_n_features,
                           pca_n_components = pca_n_components,
                           pca_outlier_n_mad = pca_outlier_n_mad,
                           umap_n_reduced_dim = umap_n_reduced_dim,
                           qc_outlier_n_mad = qc_outlier_n_mad,
                           individual_cov = individual_cov,
                           n_cores = n_cores,
                           max_mem = max_mem,
                           min_cell_per_gene = min_cell_per_gene,
                           min_individual_per_gene = min_individual_per_gene,
                           min_cells_for_itra_individual_pca_outlier = min_cells_for_itra_individual_pca_outlier,
                           vst_n_genes = vst_n_genes,
                           vst_batch_var = vst_batch_var,
                           vst_latent_var = paste(vst_latent_var, sep = ';', collapse = ';'),
                           vst_latent_var_nonreg = paste(vst_latent_var_nonreg, sep = ';', collapse = ';'),
                           pseudobulk_vst_latent_var_nonreg = paste(pseudobulk_vst_latent_var_nonreg, sep = ';', collapse = ';'),
                           vst_n_cells = paste(vst_n_cells, sep = ';', collapse = ';'), 
                           max_doublet_score = max_doublet_score,
                           doublet_irlba_nv = doublet_irlba_nv,
                           n_pcs_to_remove = n_pcs_to_remove,
                           n_top_variable_genes = n_top_variable_genes,
                           hvg_min_total_variance_quantile = hvg_min_total_variance_quantile,
                           hvg_min_count = hvg_min_count,
                           hvg_min_sample_per_gene = hvg_min_sample_per_gene,
                           stringsAsFactors = F))
settings_fn = sprintf('%s_settings.txt', out_pfx)
write_df(settings_df, file = settings_fn, row.names = T, col.names = F)
settings_data_fn = sprintf('%s_settings.RData', out_pfx)
save(list = rownames(settings_df), file = settings_data_fn)


### plot counter
plot_count = 0
get_plot_count <- function(){
  plot_count <<- plot_count + 1
  return(plot_count)
}


### variable checking
stopifnot(pca_n_components >= 2)

### read count data
verbose_print('reading count data.')
if(endsWith(x = count_fn, suffix = '.feather')){
  count_df = read_feather_df(fn = count_fn, rownames.col = 1)
} else {
  count_df = read_df(fn = count_fn)
}

sce <- SingleCellExperiment(list(counts=as.matrix(count_df)))

# clear memory
rm(count_df)
gc(reset = T)


### get genes expressed in at least n cells
verbose_print(sprintf('filtering out genes expressed in < %d cells.', min_cell_per_gene))
n_cells_per_gene = rowSums(counts(sce)>0)
sce = sce[rownames(sce)[n_cells_per_gene >= min_cell_per_gene],]
logcounts(sce) <- log2(counts(sce) + 1)


### read cell-level metadata and add to sce
verbose_print('adding cell metadata.')
metadata_df = read_df(metadata_fn, stringsAsFactors = T)
stopifnot(all(colnames(sce) %in% rownames(metadata_df)))  # all cells must have metadata
metadata_df = metadata_df[colnames(sce),,drop=F]          # identical order
for(metacol in colnames(metadata_df))
  colData(sce)[,metacol] = metadata_df[,metacol]

observed_covariates = c(colnames(metadata_df), c('total_counts', 'log10_total_counts', 'total_features_by_counts', 'pct_counts_Mt', 'phase', 'G1_score', 'G2M_score', 'S_score', 'doublet_score'))

### add gene-level metadata to sce
verbose_print('adding gene-level metadata.')
gene_annot_df = read_df(gene_annot_fn) # rownames correspond to gene ids

# gene_annot_df = read_df(gene_annot_fn, row.names = F) # rownames correspond to gene ids
# gene_annot_df0 = gene_annot_df
# gene_annot_df2 = gene_annot_df
# gene_annot_df2 = gene_annot_df2[!duplicated(gene_annot_df2$gene_name),]
# rownames(gene_annot_df2) = gene_annot_df2$gene_name
# genes = rownames(sce)
# matched_idx = match(genes, rownames(gene_annot_df2))
# gene_annot_df = gene_annot_df2[matched_idx, ,drop=F]
# rownames(gene_annot_df) = genes
# gene_annot_df$gene_name = genes
# annot_names = colnames(gene_annot_df)
# annot_names[annot_names=='strand'] = 'gene_strand'
# colnames(gene_annot_df) = annot_names
# gene_annot_df[is.na(gene_annot_df$chr), 'chr'] = 'chrNA'


gene_annot_df = gene_annot_df[rownames(sce), ,drop = F]
stopifnot(!'strand' %in% colnames(gene_annot_df))
stopifnot(all(rownames(sce) %in% rownames(gene_annot_df)))
stopifnot(all(c('gene_id', 'gene_name', 'chr', 'start_pos', 'end_pos', 'gene_strand') %in% colnames(gene_annot_df)))
for(annotcol in colnames(gene_annot_df))
  rowData(sce)[,annotcol] = gene_annot_df[,annotcol]

rowData(sce)$chr = extend_chr(rowData(sce)$chr)
rowData(sce)$gene_length = rowData(sce)$end_pos - rowData(sce)$start_pos + 1
rowData(sce)$tss = mapply(function(st, en, strand) ifelse(strand == "+", st, en), rowData(sce)$start_pos, rowData(sce)$end_pos, rowData(sce)$gene_strand)

### get genes expressed in at least n individuals
verbose_print(sprintf('filtering out genes expressed in <%d individuals.', min_individual_per_gene))
counts_mat = counts(sce)
indiv_values = colData(sce)[,individual_cov]
n_ind_per_gene = apply(counts_mat, MARGIN = 1, function(g_expr){
  n_ind = length(unique(indiv_values[g_expr>0]))
  return(n_ind)
})
genes_expressed_in_min_indiv = names(n_ind_per_gene[n_ind_per_gene>=min_individual_per_gene])
sce <- sce[genes_expressed_in_min_indiv,]

# clear memory
rm(counts_mat)
gc(reset = T)

### compute QC metrics
verbose_print('computing cell-level quality metrices.')
mito <- which(rowData(sce)$chr=="chrM")
sce <- calculateQCMetrics(sce, feature_controls=list(Mt=mito), exprs_values = "counts")

### compute mean-log-expression and variance-log-expression
verbose_print('adding mean and variance of log expression of genes across cells.')
log_expr = logcounts(sce)
mean_log_expr = rowMeans(log_expr)
variance_log_expr = apply(log_expr, 1, var)
rowData(sce)$mean_log2_expr = mean_log_expr
rowData(sce)$variance_log2_expr = variance_log_expr
# clear memory
rm(log_expr)
gc(reset = T)

### cell cylcle identification and plot
verbose_print('identifying cell cycles.')
# http://master.bioconductor.org/packages/release/workflows/vignettes/simpleSingleCell/inst/doc/umis.html#4_cell_cycle_classification
set.seed(100)
hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
ensembl_gene_ids = as.character(sapply(rowData(sce)$gene_id, function(s) strsplit(s, split = '\\.')[[1]][1]))
bpp = MulticoreParam(workers = n_cores)
assignments <- cyclone(sce, hs.pairs, gene.names=ensembl_gene_ids, BPPARAM=bpp, assay.type="counts")
colData(sce)$phase = factor(assignments$phases)
colData(sce)$G1_score = assignments$score$G1
colData(sce)$S_score = assignments$score$S
colData(sce)$G2M_score = assignments$score$G2M

### function to observe relationship between expression and covariates
plot_expr_by_cov <- function(sceobj, ncomponents = 2, plt_fn, observed_covariates, plot_func = plotUMAP){
  pdf(plt_fn)
  for(covariate in intersect(observed_covariates, colnames(colData(sceobj)))){
    print(plot_func(sceobj, ncomponents = ncomponents, colour_by = covariate))
  }
  dev.off()
}


### plot QC
plot_qc <- function(sceobj, plt_fn, observed_covariates){
  pdf(plt_fn)
  
  hist(sceobj$total_counts/1e3, xlab="Library sizes (thousands)", main="", breaks=20, col="grey80", ylab="Number of cells")
  hist(sceobj$total_features_by_counts, xlab="Number of expressed genes", main="", breaks=20, col="grey80", ylab="Number of cells")
  hist(sceobj$pct_counts_Mt, xlab="Mitochondrial proportion (%)", ylab="Number of cells", breaks=20, main="", col="grey80")
  
  for(yvar in c('total_counts', 'total_features_by_counts', 'pct_counts_Mt')){
    for(xvar in intersect(observed_covariates, colnames(colData(sceobj))))
      print(plotColData(sceobj, y=yvar, x=xvar))
  }
  
  plot(sceobj$total_features_by_counts, sceobj$total_counts/1e6, xlab="Number of expressed genes", ylab="Library size (millions)")
  plot(sceobj$total_features_by_counts, sceobj$total_counts/1e6, xlab="Number of expressed genes", ylab="Library size (millions)")
  plot(sceobj$total_features_by_counts, sceobj$pct_counts_Mt, xlab="Number of expressed genes", ylab="Mitochondrial proportion (%)")
  
  # gene counts vs gene length
  plotRowData(sceobj, x = 'gene_length', y = 'log10_total_counts', colour_by = 'pct_dropout_by_counts' )
  
  # mean-log-expression vs mean-log-variance plot
  print(plotRowData(sceobj, x='mean_log2_expr', y='variance_log2_expr', colour_by = 'log10_total_counts'))
  
  # examine most expressed genes (scater: plotHighestExprs)
  print(plotHighestExprs(sceobj, n=50) + fontsize)
  
  try(print(plotExprsFreqVsMean(sceobj)))  # it may fail [singular gradient in nls()]
  
  # print(plotScater(sceobj, colour_by = "phase", nfeatures = 300, exprs_values = "counts"))
  
  # print(plotExplanatoryVariables(sceobj))
  exp_variables = intersect(c("total_counts", "total_features_by_counts", "pct_counts_Mt", "phase"), colnames(colData(sceobj)))
  print(plotExplanatoryVariables(sceobj, variables = exp_variables))
  # print(plotExplanatoryVariables(sceobj, variables = c("total_counts", "total_features_by_counts", "pct_counts_Mt")))
  
  dev.off()  
}

###
verbose_print('plotting QC metrices.')
qc_plot_fn = sprintf('%s_%d_qc_plot.pdf', out_pfx, get_plot_count())
plot_qc(sce, plt_fn = qc_plot_fn, observed_covariates = observed_covariates)

### plot expr-cov relationship
verbose_print('plotting expression-covariate relationships in raw counts.')
set.seed(100)
sce <- runPCA(sce, ncomponents = umap_n_reduced_dim, ntop = min(c(pca_n_features, dim(sce))), detect_outliers=F, exprs_values = 'counts', method = 'irlba')
sce <- runUMAP(sce, ncomponents = pca_n_components, use_dimred = 'PCA',  detect_outliers=F)
expr_cov_plt_fn = sprintf('%s_%d_expr_cov_raw_counts_pca.pdf', out_pfx, get_plot_count())
plot_expr_by_cov(sceobj = sce, ncomponents = pca_n_components, plt_fn = expr_cov_plt_fn, observed_covariates = observed_covariates, plot_func = plotPCA)
expr_cov_plt_fn = sprintf('%s_%d_expr_cov_raw_counts_umap.pdf', out_pfx, get_plot_count())
plot_expr_by_cov(sceobj = sce, ncomponents = pca_n_components, plt_fn = expr_cov_plt_fn, observed_covariates = observed_covariates, plot_func = plotUMAP)

### drop cells based on outliers
verbose_print('identifying outliers.')
low.libsize.drop <- isOutlier(sce$total_counts, nmads=qc_outlier_n_mad, type="lower", log=TRUE)
high.libsize.drop <- isOutlier(sce$total_counts, nmads=qc_outlier_n_mad, type="higher", log=TRUE)
low.feature.drop <- isOutlier(sce$total_features_by_counts, nmads=qc_outlier_n_mad, type="lower", log=TRUE)
#high.feature.drop <- isOutlier(sce$total_features_by_counts, nmads=qc_outlier_n_mad, type="higher", log=TRUE)
high.feature.drop <- rep(FALSE, ncol(sce))
high.mito.drop <- isOutlier(sce$pct_counts_Mt, nmads=qc_outlier_n_mad, type="higher")
qc.drop = low.libsize.drop | high.libsize.drop | low.feature.drop | high.feature.drop | high.mito.drop

sce <- runPCA(sce, use_coldata=TRUE, ncomponents = pca_n_components, ntop = min(c(pca_n_features, dim(sce))), detect_outliers=F, exprs_values = 'counts')
pcs <- reducedDims(sce)$PCA_coldata
pca_outliers <- apply(pcs, MARGIN = 2, FUN = function(vals){
  isOutlier(vals, nmads=pca_outlier_n_mad, type="both", log=F)
})
pca.drop <- rowSums(pca_outliers) > 0

# drop outlier cells from same individual based on PCs
ind_pca.drop = rep(F, ncol(sce))
names(ind_pca.drop) = colnames(sce)
for (ind in unique(colData(sce)[,individual_cov])){ 
  ind_cells <- rownames(colData(sce))[colData(sce)[,individual_cov] == ind] # cell names associted with ind
  if(length(ind_cells) <= min_cells_for_itra_individual_pca_outlier)
    next
  ind_sce <- sce[, ind_cells]
  
  ind_sce <- runPCA(ind_sce, use_coldata=TRUE, ncomponents = pca_n_components, ntop = min(c(pca_n_features, dim(ind_sce))), detect_outliers=F, exprs_values = 'counts')
  ind_pcs <- reducedDims(ind_sce)$PCA_coldata
  ind_pca_outliers <- apply(ind_pcs, MARGIN = 2, FUN = function(vals){
    isOutlier(vals, nmads=pca_outlier_n_mad, type="both", log=F)
  })
  ind_pca.drop[ind_cells] = rowSums(ind_pca_outliers) > 0
  
  rm(ind_sce)
}
gc(reset = T)

colData(sce)$low.libsize.drop = low.libsize.drop
colData(sce)$high.libsize.drop = high.libsize.drop
colData(sce)$low.feature.drop = low.feature.drop
colData(sce)$high.feature.drop = high.feature.drop
colData(sce)$high.mito.drop = high.mito.drop
colData(sce)$qc.drop = qc.drop
colData(sce)$pca.drop = pca.drop
colData(sce)$ind_pca.drop = ind_pca.drop
colData(sce)$final.drop = low.libsize.drop | high.libsize.drop | low.feature.drop | high.feature.drop | high.mito.drop | pca.drop | ind_pca.drop


### plot PCs colored by outlier status and known covarites
verbose_print('plotting outliers.')
outlier_plot_fn = sprintf('%s_%d_outlier_plot.pdf', out_pfx, get_plot_count())
pdf(outlier_plot_fn)
outlier_colors = c(rgb(0,1,0,0.5), rgb(1,0,0,0.5)) # non outlier first, the outlier
outlier_pch = c(20, 18) # non outlier first, the outlier
tmp <- lapply(c('low.libsize.drop', 'high.libsize.drop', 'low.feature.drop', 'high.feature.drop', 'high.mito.drop', 'qc.drop', 'pca.drop', 'ind_pca.drop'), 
              FUN = function(metric){
                outlier_colors = c(rgb(0,1,0,0.5), rgb(1,0,0,0.5)) # non outlier first, the outlier
                outlier_pch = c(20, 18) # non outlier first, the outlier
                plot(pcs[,1], pcs[,2], 
                     xlab = 'QC PC1', ylab = 'QC PC2', 
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
tmp <- lapply(intersect(observed_covariates, colnames(colData(sce))),
              FUN = function(cur_cov){
                cov_values = colData(sce)[,cur_cov]
                pc_combinations = combn(1:pca_n_components, 2)
                if(is.numeric(cov_values)){
                  cov_order = findInterval(cov_values, sort(cov_values))
                  cell_cols = point_palette(length(cov_values))[cov_order]
                  tmp <- mapply(function(n1, n2){
                    plot(x = pcs[,n1], y = pcs[,n2], xlab = sprintf('QC PC%d', n1), ylab = sprintf('QC PC%d', n2), col = cell_cols, pch = 20, main = cur_cov)
                    legend("topright", col=point_palette(2), pch=20, legend=range(cov_values))
                    invisible(NULL)
                  }, pc_combinations[1,], pc_combinations[2,])
                } else {
                  fac_values = factor(cov_values)
                  if(length(levels(fac_values)) > 1 ){
                    tmp <- mapply(function(n1, n2){
                      plot(x = pcs[,n1], y = pcs[,n2], xlab = sprintf('QC PC%d', n1), ylab = sprintf('QC PC%d', n2), col = fac_values, pch = 20, main = cur_cov)
                      legend("topright", col=1:length(levels(fac_values)), pch=20, legend=levels(fac_values))
                      invisible(NULL)
                    }, pc_combinations[1,], pc_combinations[2,])
                  }
                }
                
                invisible(NULL)
              })

dev.off()

# clear memory
rm(pcs)
gc(reset = T)


verbose_print('plotting cells to drop.')
set.seed(100)
sce <- runPCA(sce, ncomponents = umap_n_reduced_dim, ntop = min(c(pca_n_features, dim(sce))), detect_outliers=F, exprs_values = 'counts', method = 'irlba')
sce <- runUMAP(sce, ncomponents = pca_n_components, use_dimred = 'PCA',  detect_outliers=F)
expr_cov_plt_fn = sprintf('%s_%d_expr_cov_raw_counts_pca_dropped.pdf', out_pfx, get_plot_count())
plot_expr_by_cov(sceobj = sce, ncomponents = pca_n_components, plt_fn = expr_cov_plt_fn, observed_covariates = c('final.drop'), plot_func = plotPCA)
expr_cov_plt_fn = sprintf('%s_%d_expr_cov_raw_counts_umap_dropped.pdf', out_pfx, get_plot_count())
plot_expr_by_cov(sceobj = sce, ncomponents = pca_n_components, plt_fn = expr_cov_plt_fn, observed_covariates = c('final.drop'), plot_func = plotUMAP)

keep <- !colData(sce)$final.drop
sce <- sce[,keep]  # filter low-quality cells

drop_count = c(low.libsize.drop=sum(low.libsize.drop), 
               high.libsize.drop=sum(high.libsize.drop),
               low.feature.drop=sum(low.feature.drop),
               high.feature.drop=sum(high.feature.drop),
               high.mito.drop=sum(high.mito.drop),
               pca.drop = sum(pca.drop),
               ind_pca.drop = sum(ind_pca.drop),
               Remaining=sum(keep))
verbose_print('Number of dropped cells')
print(drop_count)

verbose_print('plotting expression-covariate relationships in raw counts (just before normalization).')
set.seed(100)
sce <- runPCA(sce, ncomponents = umap_n_reduced_dim, ntop = min(c(pca_n_features, dim(sce))), detect_outliers=F, exprs_values = 'counts', method = 'irlba')
sce <- runUMAP(sce, ncomponents = pca_n_components, use_dimred = 'PCA',  detect_outliers=F)
expr_cov_plt_fn = sprintf('%s_%d_expr_cov_raw_counts_pca_before_normalization.pdf', out_pfx, get_plot_count())
plot_expr_by_cov(sceobj = sce, ncomponents = pca_n_components, plt_fn = expr_cov_plt_fn, observed_covariates = observed_covariates, plot_func = plotPCA)
expr_cov_plt_fn = sprintf('%s_%d_expr_cov_raw_counts_umap_before_normalization.pdf', out_pfx, get_plot_count())
plot_expr_by_cov(sceobj = sce, ncomponents = pca_n_components, plt_fn = expr_cov_plt_fn, observed_covariates = observed_covariates, plot_func = plotUMAP)


### normalize data
verbose_print('normalize data.')
future::plan(strategy = 'multicore', workers = n_cores)
options(future.globals.maxSize = max_mem)

set.seed(100)
vst_out <- sctransform::vst(umi = SingleCellExperiment::counts(sce), 
                            cell_attr = colData(sce), 
                            latent_var = vst_latent_var, 
                            batch_var = vst_batch_var,
                            latent_var_nonreg = vst_latent_var_nonreg,
                            n_genes = vst_n_genes,
                            n_cells = vst_n_cells,
                            return_gene_attr = FALSE, 
                            return_cell_attr = TRUE, 
                            return_corrected_umi = FALSE,
                            method = "poisson", 
                            do_regularize = TRUE,
                            residual_type = "pearson",
                            show_progress = TRUE)

# normalization_plot_fn = sprintf('%s_%d_normalization_plot.pdf', out_pfx, get_plot_count())
# pdf(normalization_plot_fn)
# sctransform::plot_model_pars(vst_out)
# ggplot(vst_out$gene_attr, aes(residual_mean)) + geom_histogram(binwidth=0.01)
# ggplot(vst_out$gene_attr, aes(residual_variance)) + geom_histogram(binwidth=0.1) + geom_vline(xintercept=1, color='red') + xlim(0, 10)
# ggplot(vst_out$gene_attr, aes(log10(gmean), residual_variance)) + geom_point(alpha=0.3, shape=16) + geom_density_2d(size = 0.3)
# 
# dev.off()

# put normalized data into sce
vst_out_fn = sprintf('%s_vst.rds', out_pfx)
saveRDS(vst_out, vst_out_fn)
corrected_genes = rownames(vst_out$y)
corrected_cells = colnames(vst_out$y)
if(length(corrected_genes) < nrow(sce)){
  verbose_print('removed genes: ')
  print(setdiff(rownames(sce), corrected_genes))
}
if(length(corrected_cells) < ncol(sce)){
  verbose_print('removed cells: ')
  print(setdiff(colnames(sce), corrected_cells))
}

sce = sce[corrected_genes, corrected_cells]
normcounts(sce) <- vst_out$y


### save corrected counts
verbose_print('temporarily saving normalized data.')
corrected_count_df = correct(vst_out, do_round = T, do_pos = T, show_progress = T)
assays(sce)[['corrected_count']] = corrected_count_df
corrected_count_fn = sprintf('%s_corrected_count.feather', out_pfx )
write_feather_df(as.data.frame(corrected_count_df), fn = corrected_count_fn, rownames.title = 'index')
intermediate_sce_fn = sprintf('%s_intermediate_sce_after_normalization.rds', out_pfx)
saveRDS(sce, file =  intermediate_sce_fn)
intermediate_image_fn = sprintf('%s_intermediate_image_after_normalization.RData', out_pfx)
save.image(file = intermediate_image_fn)

# clear memory
rm(vst_out)
gc(reset = T)


### plot after normalization
verbose_print('plotting expression-covariate relationships in normalized expression.')
set.seed(100)
sce <- runPCA(sce, exprs_values = 'normcounts', ncomponents = umap_n_reduced_dim, ntop = min(c(pca_n_features, dim(sce))), detect_outliers=F, method = 'irlba')
sce <- runUMAP(sce, ncomponents = pca_n_components, use_dimred = 'PCA',  detect_outliers=F)
expr_cov_plt_fn = sprintf('%s_%d_expr_cov_normalized_counts_pca.pdf', out_pfx, get_plot_count())
plot_expr_by_cov(sceobj = sce, ncomponents = pca_n_components, plt_fn = expr_cov_plt_fn, observed_covariates = observed_covariates, plot_func = plotPCA)
expr_cov_plt_fn = sprintf('%s_%d_expr_cov_normalized_counts_umap.pdf', out_pfx, get_plot_count())
plot_expr_by_cov(sceobj = sce, ncomponents = pca_n_components, plt_fn = expr_cov_plt_fn, observed_covariates = observed_covariates, plot_func = plotUMAP)

verbose_print('plotting expression-covariate relationships in corrected counts (after normalization).')
set.seed(100)
sce <- runPCA(sce, exprs_values = 'corrected_count', ncomponents = umap_n_reduced_dim, ntop = min(c(pca_n_features, dim(sce))), detect_outliers=F, method = 'irlba')
sce <- runUMAP(sce, ncomponents = pca_n_components, use_dimred = 'PCA',  detect_outliers=F)
expr_cov_plt_fn = sprintf('%s_%d_expr_cov_corrected_counts_pca.pdf', out_pfx, get_plot_count())
plot_expr_by_cov(sceobj = sce, ncomponents = pca_n_components, plt_fn = expr_cov_plt_fn, observed_covariates = observed_covariates, plot_func = plotPCA)
expr_cov_plt_fn = sprintf('%s_%d_expr_cov_corrected_counts_umap.pdf', out_pfx, get_plot_count())
plot_expr_by_cov(sceobj = sce, ncomponents = pca_n_components, plt_fn = expr_cov_plt_fn, observed_covariates = observed_covariates, plot_func = plotUMAP)

### observe NAs in the corrected count data
verbose_print('observing NAs in corrected count data.')
is_na_mat = is.na(corrected_count_df)
frac_na_per_gene = rowSums(is_na_mat) / ncol(is_na_mat)
frac_na_per_cell = colSums(is_na_mat) / nrow(is_na_mat)
is_zero = corrected_count_df == 0
frac_zero_per_gene = rowSums(is_zero, na.rm = T) / ncol(is_zero)
frac_zero_per_cell = colSums(is_zero, na.rm = T) / nrow(is_zero)

corrected_count_plt_fn = sprintf('%s_%d_corrected_count.pdf', out_pfx, get_plot_count())
pdf(corrected_count_plt_fn)
hist(frac_na_per_gene[frac_na_per_gene>0], breaks = 30, xlab = 'Fraction of NA per gene', main = sprintf('Total %d (out of %d) genes with NA', sum(frac_na_per_gene>0), length(frac_na_per_gene)))
hist(frac_na_per_cell[frac_na_per_cell>0], breaks = 30, xlab = 'Fraction of NA per cell', main = sprintf('Total %d (out of %d) cells with NA', sum(frac_na_per_cell>0), length(frac_na_per_cell)))
hist(frac_zero_per_gene[frac_zero_per_gene>0], breaks = 30, xlab = 'Fraction of zero per gene', main = sprintf('Total %d (out of %d) genes with zero', sum(frac_zero_per_gene>0), length(frac_zero_per_gene)))
hist(frac_zero_per_cell[frac_zero_per_cell>0], breaks = 30, xlab = 'Fraction of zero per cell', main = sprintf('Total %d (out of %d) cells with NA', sum(frac_zero_per_cell>0), length(frac_zero_per_cell)))
dev.off()

### filter genes/cells based on NA or Zero in corrected counts
verbose_print('filtering genes with either all 0\'s or some NAs in corrected count data.')
genes_wo_na = names(frac_na_per_gene[frac_na_per_gene==0])
genes_with_enough_cells = names(frac_zero_per_gene[frac_zero_per_gene <= (1-min_cell_per_gene/ncol(is_zero))])
cells_with_enough_genes = names(frac_zero_per_cell[frac_zero_per_cell <= (1-1/nrow(is_zero))])
sce = sce[intersect(genes_wo_na, genes_with_enough_cells), cells_with_enough_genes]

# clear memory
rm(corrected_count_df, is_na_mat, is_zero)
gc(reset = T)


### TODO: remove doublets (using doubletCells())
verbose_print('computing doublet scores.')
batches = unique(colData(sce)[,vst_batch_var])
doublet_scores = lapply(batches, function(batch){
  in_batch = colData(sce)[,vst_batch_var] == batch
  if(sum(in_batch) < 2)
    next
  
  batch_sce = sce[,in_batch]
  doublet_score <- doubletCells(batch_sce, assay.type="corrected_count", approximate = F) # irlba.args = list(nv = doublet_irlba_nv)
  names(doublet_score) = colnames(batch_sce)
  # clear memory
  rm(batch_sce)
  gc(reset = T)
  return(doublet_score)
})

doublet_scores = unlist(doublet_scores, recursive = T)
doublet_scores = doublet_scores[colnames(sce)]
colData(sce)$doublet_score = doublet_scores

verbose_print('plotting expression-doubleScore relationships in corrected counts (after normalization).')
set.seed(100)
sce <- runPCA(sce, exprs_values = 'normcounts', ncomponents = umap_n_reduced_dim, ntop = min(c(pca_n_features, dim(sce))), detect_outliers=F, method = 'irlba')
sce <- runUMAP(sce, ncomponents = pca_n_components, use_dimred = 'PCA',  detect_outliers=F)
expr_cov_plt_fn = sprintf('%s_%d_expr_cov_normalized_counts_doublet_pca.pdf', out_pfx, get_plot_count())
plot_expr_by_cov(sceobj = sce, ncomponents = pca_n_components, plt_fn = expr_cov_plt_fn, observed_covariates = c('doublet_score', 'total_counts', 'log10_total_counts'), plot_func = plotPCA)
expr_cov_plt_fn = sprintf('%s_%d_expr_cov_normalized_counts_doublet_umap.pdf', out_pfx, get_plot_count())
plot_expr_by_cov(sceobj = sce, ncomponents = pca_n_components, plt_fn = expr_cov_plt_fn, observed_covariates = c('doublet_score', 'total_counts', 'log10_total_counts'), plot_func = plotUMAP)

# TODO: remove doublets? too many cells, preferentially high library size cells, are filtered out.
# sce = sce[,colData(sce)<=max_doublet_score]  # filter doublets


### remove technical noise
verbose_print('analyzing biological and technical noise.')
# neds normalized logcounts
verbose_print('estimating technical noise')
tech_sce = SingleCellExperiment(list(counts=counts(sce), logcounts=logcounts(sce)))
tech_sce <- computeSumFactors(tech_sce)
tech_sce <- normalize(tech_sce)
#TODO: use normcounts?
tech_var_fit <- trendVar(tech_sce, parametric=TRUE, use.spikes=FALSE, loess.args=list(span=0.2), assay.type = 'logcounts')
tech_var_decomposed <- decomposeVar(tech_sce, tech_var_fit)
verbose_print('denoising technical noise')
sce = denoisePCA(x = sce, 
                 technical = tech_var_decomposed, 
                 value = 'lowrank', 
                 assay.type = 'normcounts',
                 approximate = T)

# clear memory
rm(tech_sce, tech_var_fit)
gc(reset = T)

### expr-cov plots after denoising
verbose_print('plotting expression-covariate relationships in denoised data.')
set.seed(100)
sce <- runPCA(sce, ncomponents = umap_n_reduced_dim, ntop = min(c(pca_n_features, dim(sce))), detect_outliers=F, exprs_values = 'lowrank', method = 'irlba')
sce <- runUMAP(sce, ncomponents = pca_n_components, use_dimred = 'PCA',  detect_outliers=F)
expr_cov_plt_fn = sprintf('%s_%d_expr_cov_denoised_pca.pdf', out_pfx, get_plot_count())
plot_expr_by_cov(sceobj = sce, ncomponents = pca_n_components, plt_fn = expr_cov_plt_fn, observed_covariates = observed_covariates, plot_func = plotPCA)
expr_cov_plt_fn = sprintf('%s_%d_expr_cov_denoised_umap.pdf', out_pfx, get_plot_count())
plot_expr_by_cov(sceobj = sce, ncomponents = pca_n_components, plt_fn = expr_cov_plt_fn, observed_covariates = observed_covariates, plot_func = plotUMAP)

### remove PCs and save in normalized_pc_corrected
verbose_print(sprintf('removing %d PCs from normalized data.', n_pcs_to_remove))
sce <- runPCA(sce, ncomponents = n_pcs_to_remove, ntop = nrow(sce), detect_outliers=F, exprs_values = 'normcounts', method = 'irlba')
pcs = reducedDim(sce, "PCA")
pc_removed_expr = removeBatchEffect(normcounts(sce), covariates = pcs)
assays(sce)[['pc_corrected']] = pc_removed_expr
for(ipc in 1:n_pcs_to_remove)
  colData(sce)[,sprintf('PC%d', ipc)] = pcs[,ipc]
# clear memory
rm(pcs, pc_removed_expr)
gc(reset = T)

### expr-cov plots after normalized_pc_corrected
verbose_print('plotting expression-covariate relationships in pc-removed data.')
set.seed(100)
sce <- runPCA(sce, ncomponents = umap_n_reduced_dim, ntop = min(c(pca_n_features, dim(sce))), detect_outliers=F, exprs_values = 'pc_corrected', method = 'irlba')
sce <- runUMAP(sce, ncomponents = pca_n_components, use_dimred = 'PCA',  detect_outliers=F)
expr_cov_plt_fn = sprintf('%s_%d_expr_cov_pc_corrected_pca.pdf', out_pfx, get_plot_count())
plot_expr_by_cov(sceobj = sce, ncomponents = pca_n_components, plt_fn = expr_cov_plt_fn, observed_covariates = observed_covariates, plot_func = plotPCA)
expr_cov_plt_fn = sprintf('%s_%d_expr_cov_pc_corrected_umap.pdf', out_pfx, get_plot_count())
plot_expr_by_cov(sceobj = sce, ncomponents = pca_n_components, plt_fn = expr_cov_plt_fn, observed_covariates = observed_covariates, plot_func = plotUMAP)


### function to write sce data into files
write_sce_data <- function(sceobj, assay.type, assay.fn, cell.meta.fn = sprintf('%s.cell.meta.txt', assay.fn), gene.meta.fn = sprintf('%s.gene.meta.txt', assay.fn)){
  # write assay (expr) data
  assay_df = as.data.frame(assays(sceobj)[[assay.type]])
  if(is.character(assay.fn)){
    if(endsWith(x = assay.fn, suffix = 'feather')){
      write_feather_df(assay_df, fn = assay.fn, rownames.title = 'index')
    } else {
      write_df(assay_df, file = assay.fn)
    }
  }
    
  # write cell metadata
  if(is.character(cell.meta.fn))
    write_df(colData(sceobj), file = cell.meta.fn)
  
  # write gene metadata
  if(is.character(gene.meta.fn))
    write_df(rowData(sceobj), file = gene.meta.fn)

  invisible(NULL)
}

### select genes with positive biological variances
verbose_print('saving data with (positive) biological variance.')
pos_bio_variance_decomposed = tech_var_decomposed[tech_var_decomposed$bio > 0, ,drop=F]
pos_bio_var_genes = rownames(pos_bio_variance_decomposed)
pos_bio_sce = sce[pos_bio_var_genes, ]
pos_bio_denoised_fn = sprintf('%s_biovar_normalized_denoised.feather', out_pfx)
pos_bio_counts_fn = sprintf('%s_biovar_counts.feather', out_pfx)
pos_bio_normalized_fn = sprintf('%s_biovar_normalized.feather', out_pfx)
pos_bio_pc_corrected_fn = sprintf('%s_biovar_normalized_pc_corrected.feather', out_pfx)
pos_bio_cell_meta_fn = sprintf('%s_biovar_cell_meta.txt', out_pfx)
pos_bio_gene_meta_fn = sprintf('%s_biovar_gene_meta.txt', out_pfx)
write_sce_data(pos_bio_sce, assay.type = 'counts', assay.fn = pos_bio_counts_fn, cell.meta.fn = NULL, gene.meta.fn = NULL)
write_sce_data(pos_bio_sce, assay.type = 'normcounts', assay.fn = pos_bio_normalized_fn, cell.meta.fn = pos_bio_cell_meta_fn, gene.meta.fn = pos_bio_gene_meta_fn)
write_sce_data(pos_bio_sce, assay.type = 'lowrank', assay.fn = pos_bio_denoised_fn, cell.meta.fn = NULL, gene.meta.fn = NULL)
write_sce_data(pos_bio_sce, assay.type = 'pc_corrected', assay.fn = pos_bio_pc_corrected_fn, cell.meta.fn = NULL, gene.meta.fn = NULL)

### select top variable genes
verbose_print(sprintf('saving data with %d highly variable genes.', n_top_variable_genes))

# get genes with min median expression
pos_bio_counts_mat = counts(pos_bio_sce)
min_expressed_pos_bio_counts_mat = filter_expr_by_tpm_read(expr.df = pos_bio_counts_mat, tpm.df = pos_bio_counts_mat, count.df = pos_bio_counts_mat, min.tpm = hvg_min_count, min.count = hvg_min_count, min.samples = hvg_min_sample_per_gene)
min_exprssed_genes = rownames(min_expressed_pos_bio_counts_mat)
candidate_hvg_variance_decomposed = pos_bio_variance_decomposed[
                                      rownames(pos_bio_variance_decomposed) %in% min_exprssed_genes & pos_bio_variance_decomposed$total >= quantile(pos_bio_variance_decomposed$total, probs = hvg_min_total_variance_quantile), , drop = F]
candidate_hvg_variance_decomposed = candidate_hvg_variance_decomposed[with(candidate_hvg_variance_decomposed, order(FDR, -bio)),,drop=F]
hvg = rownames(candidate_hvg_variance_decomposed)[1:min(n_top_variable_genes, nrow(candidate_hvg_variance_decomposed))]
hvg_sce = sce[hvg, ]
hvg_denoised_fn = sprintf('%s_hvg%d_normalized_denoised.feather', out_pfx, n_top_variable_genes)
hvg_counts_fn = sprintf('%s_hvg%d_counts.feather', out_pfx, n_top_variable_genes)
hvg_normalized_fn = sprintf('%s_hvg%d_normalized.feather', out_pfx, n_top_variable_genes)
hvg_pc_corrected_fn = sprintf('%s_hvg%d_normalized_pc_corrected.feather', out_pfx, n_top_variable_genes)
hvg_cell_meta_fn = sprintf('%s_hvg%d_cell_meta.txt', out_pfx, n_top_variable_genes)
hvg_gene_meta_fn = sprintf('%s_hvg%d_gene_meta.txt', out_pfx, n_top_variable_genes)
write_sce_data(hvg_sce, assay.type = 'lowrank', assay.fn = hvg_denoised_fn, cell.meta.fn = hvg_cell_meta_fn, gene.meta.fn = hvg_gene_meta_fn)
write_sce_data(hvg_sce, assay.type = 'counts', assay.fn = hvg_counts_fn, cell.meta.fn = NULL, gene.meta.fn = NULL)
write_sce_data(hvg_sce, assay.type = 'normcounts', assay.fn = hvg_normalized_fn, cell.meta.fn = NULL, gene.meta.fn = NULL)
write_sce_data(hvg_sce, assay.type = 'pc_corrected', assay.fn = hvg_pc_corrected_fn, cell.meta.fn = NULL, gene.meta.fn = NULL)


### expr-cov plots after denoising, using only top variable genes
verbose_print('plotting expression-covariate relationships in highly variable genes in denoised data.')
set.seed(100)
hvg_sce <- runPCA(hvg_sce, ncomponents = umap_n_reduced_dim, ntop = min(c(pca_n_features, dim(hvg_sce))), detect_outliers=F, exprs_values = 'lowrank', method = 'irlba')
hvg_sce <- runUMAP(hvg_sce, ncomponents = pca_n_components, use_dimred = 'PCA',  detect_outliers=F)
expr_cov_plt_fn = sprintf('%s_%d_expr_cov_denoised_hvg%d_pca.pdf', out_pfx, get_plot_count(), n_top_variable_genes)
plot_expr_by_cov(sceobj = hvg_sce, ncomponents = pca_n_components, plt_fn = expr_cov_plt_fn, observed_covariates = observed_covariates, plot_func = plotPCA)
expr_cov_plt_fn = sprintf('%s_%d_expr_cov_denoised_hvg%d_umap.pdf', out_pfx, get_plot_count(), n_top_variable_genes)
plot_expr_by_cov(sceobj = hvg_sce, ncomponents = pca_n_components, plt_fn = expr_cov_plt_fn, observed_covariates = observed_covariates, plot_func = plotUMAP)

### expr-cov plots after pc-correction, using only top variable genes
verbose_print('plotting expression-covariate relationships in highly variable genes in pc-removed data.')
set.seed(100)
hvg_sce <- runPCA(hvg_sce, ncomponents = umap_n_reduced_dim, ntop = min(c(pca_n_features, dim(hvg_sce))), detect_outliers=F, exprs_values = 'pc_corrected', method = 'irlba')
hvg_sce <- runUMAP(hvg_sce, ncomponents = pca_n_components, use_dimred = 'PCA',  detect_outliers=F)
expr_cov_plt_fn = sprintf('%s_%d_expr_cov_pc_corrected_hvg%d_pca.pdf', out_pfx, get_plot_count(), n_top_variable_genes)
plot_expr_by_cov(sceobj = hvg_sce, ncomponents = pca_n_components, plt_fn = expr_cov_plt_fn, observed_covariates = observed_covariates, plot_func = plotPCA)
expr_cov_plt_fn = sprintf('%s_%d_expr_cov_pc_corrected_hvg%d_umap.pdf', out_pfx, get_plot_count(), n_top_variable_genes)
plot_expr_by_cov(sceobj = hvg_sce, ncomponents = pca_n_components, plt_fn = expr_cov_plt_fn, observed_covariates = observed_covariates, plot_func = plotUMAP)

###
verbose_print('plotting expression patterns of highly variable genes.')
hvg_expr_plot_fn = sprintf('%s_%d_hvg%d_expr.pdf', out_pfx, get_plot_count(), n_top_variable_genes)
pdf(hvg_expr_plot_fn)
plotExpression(sce, features=hvg[1:10], exprs_values = 'counts') + fontsize
plotExpression(sce, features=hvg[1:10], exprs_values = 'logcounts') + fontsize
plotExpression(sce, features=hvg[1:10], exprs_values = 'normcounts') + fontsize
plotExpression(sce, features=hvg[1:10], exprs_values = 'lowrank') + fontsize
dev.off()


### save workspace after single-cell data processing
sc_image_fn = sprintf('%s_intermediate_image_after_single_cell_processing.RData', out_pfx)
save.image(file = sc_image_fn)

# clear memory
rm(pos_bio_sce, pos_bio_counts_mat, min_expressed_pos_bio_counts_mat)
gc(reset = T)

#################################################################################
### Make pseudo-bulk data & plot pseudo-bulk cells with covariates
verbose_print('computing pseudo-bulk expression.')
indiv_values = colData(sce)[,individual_cov]   # metadata_df[colnames(sce), individual_cov]
uniq_indiv = unique(indiv_values)
corrected_count_mat = assays(sce)[['corrected_count']]
pseudobulk_count_df = sapply(uniq_indiv, function(indiv_id){
  ind_count_df = corrected_count_mat[,indiv_values == indiv_id, drop = F]
  ind_count_sum_df = rowSums(ind_count_df)
  return(ind_count_sum_df)
})
colnames(pseudobulk_count_df) = uniq_indiv

pseudobulk_sce <- SingleCellExperiment(list(counts=as.matrix(pseudobulk_count_df),
                                            logcounts = log2(1+as.matrix(pseudobulk_count_df))))

# clear memory
rm(corrected_count_mat, pseudobulk_count_df)
gc(reset = T)

### make subject level covariates
verbose_print('generating individual-level covariates from cell-level covariates.')
sub_cov_df = aggregate(colData(sce)[,colnames(metadata_df)], by = list(colData(sce)[,individual_cov]), function(x){
  if(is.numeric(x))
    return(mean(x, na.rm = T))
  if(length(unique(x)) == 1)
    return(x[1])
  return(NA)
})
sub_cov_df = sub_cov_df[,-1,drop=F]
rownames(sub_cov_df) = sub_cov_df[,individual_cov]
sub_cov_df = sub_cov_df[colnames(pseudobulk_sce),]
for(metacol in colnames(sub_cov_df))
  colData(pseudobulk_sce)[,metacol] = sub_cov_df[,metacol]


rowData(pseudobulk_sce)$chr = rowData(sce)[rownames(pseudobulk_sce),'chr']
rowData(pseudobulk_sce)$gene_length = rowData(sce)[rownames(pseudobulk_sce),'gene_length']
mito <- which(rowData(pseudobulk_sce)$chr=="chrM")
pseudobulk_sce <- calculateQCMetrics(pseudobulk_sce, feature_controls=list(Mt=mito))

log_expr = logcounts(pseudobulk_sce)
mean_log_expr = rowMeans(log_expr)
variance_log_expr = apply(log_expr, 1, var)
rowData(pseudobulk_sce)$mean_log2_expr = mean_log_expr
rowData(pseudobulk_sce)$variance_log2_expr = variance_log_expr
# clear memory
rm(log_expr)
gc(reset = T)


### plot pseudobulk QC
verbose_print('plotting pseudobulk QC metrices.')
pseudo_qc_plot_fn = sprintf('%s_%d_pseudobulk_qc_plot.pdf', out_pfx, get_plot_count())
plot_qc(pseudobulk_sce, plt_fn = pseudo_qc_plot_fn, observed_covariates = observed_covariates)

# plot by covariates
verbose_print('plotting expression-covariate relationships in raw pseudobulk count data.')
set.seed(100)
pseudobulk_sce <- runPCA(pseudobulk_sce, ncomponents = umap_n_reduced_dim, ntop = min(c(pca_n_features, dim(pseudobulk_sce))), detect_outliers=F, exprs_values = 'counts', method = 'irlba')
pseudobulk_sce <- runUMAP(pseudobulk_sce, ncomponents = pca_n_components, use_dimred = 'PCA',  detect_outliers=F)
pseudo_expr_plot_fn = sprintf('%s_%d_pseudobulk_expr_cov_plot_before_normalization_pca.pdf', out_pfx, get_plot_count())
plot_expr_by_cov(pseudobulk_sce, plt_fn = pseudo_expr_plot_fn, observed_covariates = observed_covariates, plot_func = plotPCA)
pseudo_expr_plot_fn = sprintf('%s_%d_pseudobulk_expr_cov_plot_before_normalization_umap.pdf', out_pfx, get_plot_count())
plot_expr_by_cov(pseudobulk_sce, plt_fn = pseudo_expr_plot_fn, observed_covariates = observed_covariates, plot_func = plotUMAP)

# vst
verbose_print('normalizing pseudobulk data.')
set.seed(100)
pseudo_vst_out <- sctransform::vst(umi = SingleCellExperiment::counts(pseudobulk_sce), 
                            cell_attr = colData(pseudobulk_sce), 
                            latent_var = vst_latent_var, 
                            batch_var = NULL,
                            latent_var_nonreg = pseudobulk_vst_latent_var_nonreg,
                            n_genes = vst_n_genes,
                            n_cells = NULL,
                            return_gene_attr = FALSE, 
                            return_cell_attr = TRUE, 
                            return_corrected_umi = FALSE,
                            method = "poisson", 
                            do_regularize = TRUE,
                            residual_type = "pearson",
                            show_progress = TRUE)

pseudo_corrected_genes = rownames(pseudo_vst_out$y)
pseudo_corrected_cells = colnames(pseudo_vst_out$y)
if(length(pseudo_corrected_genes) < nrow(pseudobulk_sce)){
  verbose_print('removed genes: ')
  print(setdiff(rownames(pseudobulk_sce), pseudo_corrected_genes))
}
if(length(pseudo_corrected_cells) < ncol(pseudobulk_sce)){
  verbose_print('removed cells: ')
  print(setdiff(colnames(pseudobulk_sce), pseudo_corrected_cells))
}

pseudobulk_sce = pseudobulk_sce[pseudo_corrected_genes, pseudo_corrected_cells]
normcounts(pseudobulk_sce) <- pseudo_vst_out$y


# save pseudo_vst_out
pseudo_vst_out_fn = sprintf('%s_pseudobulk_vst.rds', out_pfx)
saveRDS(pseudo_vst_out, pseudo_vst_out_fn)

# clear memory
rm(pseudo_vst_out)
gc(reset = T)

# plot by covariates
verbose_print('plotting expression-covariate relationships in normalized pseudobulk data.')
set.seed(100)
pseudobulk_sce <- runPCA(pseudobulk_sce, ncomponents = umap_n_reduced_dim, ntop = min(c(pca_n_features, dim(pseudobulk_sce))), detect_outliers=F, exprs_values = 'normcounts', method = 'irlba')
pseudobulk_sce <- runUMAP(pseudobulk_sce, ncomponents = pca_n_components, use_dimred = 'PCA',  detect_outliers=F)
pseudo_expr_plot_fn = sprintf('%s_%d_pseudobulk_expr_cov_plot_after_normalization_pca.pdf', out_pfx, get_plot_count())
plot_expr_by_cov(pseudobulk_sce, plt_fn = pseudo_expr_plot_fn, observed_covariates = observed_covariates, plot_func = plotPCA)
pseudo_expr_plot_fn = sprintf('%s_%d_pseudobulk_expr_cov_plot_after_normalization_umap.pdf', out_pfx, get_plot_count())
plot_expr_by_cov(pseudobulk_sce, plt_fn = pseudo_expr_plot_fn, observed_covariates = observed_covariates, plot_func = plotUMAP)

### remove technical noise
verbose_print('analyzing contribution of biological and technical noise in pseudobulk data.')
# neds normalized logcounts
verbose_print('estimating technical noise')
pseudo_tech_sce = SingleCellExperiment(list(counts=counts(pseudobulk_sce),
                                            logcounts=logcounts(pseudobulk_sce)))
pseudo_tech_sce <- computeSumFactors(pseudo_tech_sce)
pseudo_tech_sce <- normalize(pseudo_tech_sce)
pseudo_tech_var_fit <- trendVar(pseudo_tech_sce, parametric=TRUE, use.spikes=FALSE, loess.args=list(span=0.2), assay.type = 'logcounts')
pseudo_tech_var_decomposed <- decomposeVar(pseudo_tech_sce, pseudo_tech_var_fit)
verbose_print('denoising technical noise')
pseudobulk_sce = denoisePCA(x = pseudobulk_sce, 
                 technical = pseudo_tech_var_decomposed, 
                 value = 'lowrank', 
                 assay.type = 'normcounts',
                 approximate = F)

### expr-cov plots after denoising
verbose_print('plotting expression-covariate relationships in denoised pseudobulk data.')
set.seed(100)
pseudobulk_sce <- runPCA(pseudobulk_sce, ncomponents = umap_n_reduced_dim, ntop = min(c(pca_n_features, dim(pseudobulk_sce))), detect_outliers=F, exprs_values = 'lowrank', method = 'irlba')
pseudobulk_sce <- runUMAP(pseudobulk_sce, ncomponents = pca_n_components, use_dimred = 'PCA',  detect_outliers=F)
expr_cov_plt_fn = sprintf('%s_%d_pseudobulk_expr_cov_denoised_pca.pdf', out_pfx, get_plot_count())
plot_expr_by_cov(sceobj = pseudobulk_sce, ncomponents = pca_n_components, plt_fn = expr_cov_plt_fn, observed_covariates = observed_covariates, plot_func = plotPCA)
expr_cov_plt_fn = sprintf('%s_%d_pseudobulk_expr_cov_denoised_umap.pdf', out_pfx, get_plot_count())
plot_expr_by_cov(sceobj = pseudobulk_sce, ncomponents = pca_n_components, plt_fn = expr_cov_plt_fn, observed_covariates = observed_covariates, plot_func = plotUMAP)

### remove PCs and save in normalized_pc_corrected
verbose_print(sprintf('removing %d PCs from normalized pseudobulk data.', n_pcs_to_remove))
pseudobulk_sce <- runPCA(pseudobulk_sce, ncomponents = n_pcs_to_remove, ntop = nrow(pseudobulk_sce), detect_outliers=F, exprs_values = 'normcounts', method = 'irlba')
pcs = reducedDim(pseudobulk_sce, "PCA")
pc_removed_expr = removeBatchEffect(normcounts(pseudobulk_sce), covariates = pcs)
assays(pseudobulk_sce)[['pc_corrected']] = pc_removed_expr
for(ipc in 1:n_pcs_to_remove)
  colData(pseudobulk_sce)[,sprintf('PC%d', ipc)] = pcs[,ipc]
# clear memory
rm(pseudo_tech_sce, pseudo_tech_var_fit, pcs, pc_removed_expr)
gc(reset = T)

### expr-cov plots after normalized_pc_corrected
verbose_print('plotting expression-covariate relationships in pc-removed pseudobulk data.')
set.seed(100)
pseudobulk_sce <- runPCA(pseudobulk_sce, ncomponents = umap_n_reduced_dim, ntop = min(c(pca_n_features, dim(pseudobulk_sce))), detect_outliers=F, exprs_values = 'pc_corrected', method = 'irlba')
pseudobulk_sce <- runUMAP(pseudobulk_sce, ncomponents = pca_n_components, use_dimred = 'PCA',  detect_outliers=F)
expr_cov_plt_fn = sprintf('%s_%d_pseudobulk_expr_cov_pc_corrected_pca.pdf', out_pfx, get_plot_count())
plot_expr_by_cov(sceobj = pseudobulk_sce, ncomponents = pca_n_components, plt_fn = expr_cov_plt_fn, observed_covariates = observed_covariates, plot_func = plotPCA)
expr_cov_plt_fn = sprintf('%s_%d_pseudobulk_expr_cov_pc_corrected_umap.pdf', out_pfx, get_plot_count())
plot_expr_by_cov(sceobj = pseudobulk_sce, ncomponents = pca_n_components, plt_fn = expr_cov_plt_fn, observed_covariates = observed_covariates, plot_func = plotUMAP)

### select genes with positive biological variances
verbose_print('saving pseudobulk data of genes with positive biologically variance.')
pseudo_pos_bio_variance_decomposed = pseudo_tech_var_decomposed[pseudo_tech_var_decomposed$bio > 0, ,drop=F]
pseudo_pos_bio_var_genes = rownames(pseudo_pos_bio_variance_decomposed)
pseudo_pos_bio_sce = pseudobulk_sce[pseudo_pos_bio_var_genes, ]
pseudo_pos_bio_denoised_fn = sprintf('%s_pseudobulk_biovar_normalized_denoised.feather', out_pfx)
pseudo_pos_bio_counts_fn = sprintf('%s_pseudobulk_biovar_counts.feather', out_pfx)
pseudo_pos_bio_normalized_fn = sprintf('%s_pseudobulk_biovar_normalized.feather', out_pfx)
pseudo_pos_bio_pc_corrected_fn = sprintf('%s_pseudobulk_biovar_normalized_pc_corrected.feather', out_pfx)
pseudo_pos_bio_cell_meta_fn = sprintf('%s_pseudobulk_biovar_cell_meta.txt', out_pfx)
pseudo_pos_bio_gene_meta_fn = sprintf('%s_pseudobulk_biovar_gene_meta.txt', out_pfx)
write_sce_data(pseudo_pos_bio_sce, assay.type = 'counts', assay.fn = pseudo_pos_bio_counts_fn, cell.meta.fn = NULL, gene.meta.fn = NULL)
write_sce_data(pseudo_pos_bio_sce, assay.type = 'normcounts', assay.fn = pseudo_pos_bio_normalized_fn, cell.meta.fn = pseudo_pos_bio_cell_meta_fn, gene.meta.fn = pseudo_pos_bio_gene_meta_fn)
write_sce_data(pseudo_pos_bio_sce, assay.type = 'lowrank', assay.fn = pseudo_pos_bio_denoised_fn, cell.meta.fn = NULL, gene.meta.fn = NULL)
write_sce_data(pseudo_pos_bio_sce, assay.type = 'pc_corrected', assay.fn = pseudo_pos_bio_pc_corrected_fn, cell.meta.fn = NULL, gene.meta.fn = NULL)

### select top variable genes within min-expressed genes
verbose_print(sprintf('saving pseudobulk data with %d highly variable genes.', n_top_variable_genes))
pseudo_pos_bio_counts_mat = counts(pseudo_pos_bio_sce)
min_expressed_pseudo_pos_bio_counts_mat = filter_expr_by_tpm_read(expr.df = pseudo_pos_bio_counts_mat, tpm.df = pseudo_pos_bio_counts_mat, count.df = pseudo_pos_bio_counts_mat, min.tpm = hvg_min_count, min.count = hvg_min_count, min.samples = hvg_min_sample_per_gene)
pseudo_min_exprssed_genes = rownames(min_expressed_pseudo_pos_bio_counts_mat)
pseudo_candidate_hvg_variance_decomposed = pseudo_pos_bio_variance_decomposed[
  rownames(pseudo_pos_bio_variance_decomposed) %in% pseudo_min_exprssed_genes & pseudo_pos_bio_variance_decomposed$total >= quantile(pseudo_pos_bio_variance_decomposed$total, probs = hvg_min_total_variance_quantile), , drop = F]
pseudo_candidate_hvg_variance_decomposed = pseudo_candidate_hvg_variance_decomposed[with(pseudo_candidate_hvg_variance_decomposed, order(FDR, -bio)),,drop=F]
pseudo_hvg = rownames(pseudo_candidate_hvg_variance_decomposed)[1:min(n_top_variable_genes, nrow(pseudo_candidate_hvg_variance_decomposed))]
pseudo_hvg_sce = pseudobulk_sce[pseudo_hvg, ]
pseudo_hvg_denoised_fn = sprintf('%s_pseudobulk_hvg%d_normalized_denoised.feather', out_pfx, n_top_variable_genes)
pseudo_hvg_counts_fn = sprintf('%s_pseudobulk_hvg%d_counts.feather', out_pfx, n_top_variable_genes)
pseudo_hvg_normalized_fn = sprintf('%s_pseudobulk_hvg%d_normalized.feather', out_pfx, n_top_variable_genes)
pseudo_hvg_pc_corrected_fn = sprintf('%s_pseudobulk_hvg%d_normalized_pc_corrected.feather', out_pfx, n_top_variable_genes)
pseudo_hvg_cell_meta_fn = sprintf('%s_pseudobulk_hvg%d_cell_meta.txt', out_pfx, n_top_variable_genes)
pseudo_hvg_gene_meta_fn = sprintf('%s_pseudobulk_hvg%d_gene_meta.txt', out_pfx, n_top_variable_genes)
write_sce_data(pseudo_hvg_sce, assay.type = 'lowrank', assay.fn = pseudo_hvg_denoised_fn, cell.meta.fn = pseudo_hvg_cell_meta_fn, gene.meta.fn = pseudo_hvg_gene_meta_fn)
write_sce_data(pseudo_hvg_sce, assay.type = 'counts', assay.fn = pseudo_hvg_counts_fn, cell.meta.fn = NULL, gene.meta.fn = NULL)
write_sce_data(pseudo_hvg_sce, assay.type = 'normcounts', assay.fn = pseudo_hvg_normalized_fn, cell.meta.fn = NULL, gene.meta.fn = NULL)
write_sce_data(pseudo_hvg_sce, assay.type = 'pc_corrected', assay.fn = pseudo_hvg_pc_corrected_fn, cell.meta.fn = NULL, gene.meta.fn = NULL)

### expr-cov plots after denoising, using only top variable genes
verbose_print('plotting expression-covariate relationships in denoised pseudobulk data with highly variable genes only.')
set.seed(100)
pseudo_hvg_sce <- runPCA(pseudo_hvg_sce, ncomponents = umap_n_reduced_dim, ntop = min(c(pca_n_features, dim(pseudo_hvg_sce))), detect_outliers=F, exprs_values = 'lowrank', method = 'irlba')
pseudo_hvg_sce <- runUMAP(pseudo_hvg_sce, ncomponents = pca_n_components, use_dimred = 'PCA',  detect_outliers=F)
expr_cov_plt_fn = sprintf('%s_%d_pseudobulk_expr_cov_denoised_hvg%d_pca.pdf', out_pfx, get_plot_count(), n_top_variable_genes)
plot_expr_by_cov(sceobj = pseudo_hvg_sce, ncomponents = pca_n_components, plt_fn = expr_cov_plt_fn, observed_covariates = observed_covariates, plot_func = plotPCA)
expr_cov_plt_fn = sprintf('%s_%d_pseudobulk_expr_cov_denoised_hvg%d_umap.pdf', out_pfx, get_plot_count(), n_top_variable_genes)
plot_expr_by_cov(sceobj = pseudo_hvg_sce, ncomponents = pca_n_components, plt_fn = expr_cov_plt_fn, observed_covariates = observed_covariates, plot_func = plotUMAP)

### expr-cov plots after pc-correction, using only top variable genes
verbose_print('plotting expression-covariate relationships in pc-removed pseudobulk data with highly variable genes only.')
set.seed(100)
pseudo_hvg_sce <- runPCA(pseudo_hvg_sce, ncomponents = umap_n_reduced_dim, ntop = min(c(pca_n_features, dim(pseudo_hvg_sce))), detect_outliers=F, exprs_values = 'pc_corrected', method = 'irlba')
pseudo_hvg_sce <- runUMAP(pseudo_hvg_sce, ncomponents = pca_n_components, use_dimred = 'PCA',  detect_outliers=F)
expr_cov_plt_fn = sprintf('%s_%d_pseudobulk_expr_cov_pc_corrected_hvg%d_pca.pdf', out_pfx, get_plot_count(), n_top_variable_genes)
plot_expr_by_cov(sceobj = pseudo_hvg_sce, ncomponents = pca_n_components, plt_fn = expr_cov_plt_fn, observed_covariates = observed_covariates, plot_func = plotPCA)
expr_cov_plt_fn = sprintf('%s_%d_pseudobulk_expr_cov_pc_corrected_hvg%d_umap.pdf', out_pfx, get_plot_count(), n_top_variable_genes)
plot_expr_by_cov(sceobj = pseudo_hvg_sce, ncomponents = pca_n_components, plt_fn = expr_cov_plt_fn, observed_covariates = observed_covariates, plot_func = plotUMAP)

###
verbose_print('plotting pseudobulk expression of highly variable genes.')
pseudo_hvg_expr_plot_fn = sprintf('%s_%d_pseudobulk_hvg%d_expr.pdf', out_pfx, get_plot_count(), n_top_variable_genes)
pdf(pseudo_hvg_expr_plot_fn)
print(plotExpression(pseudobulk_sce, features=pseudo_hvg[1:10], exprs_values = 'counts') + fontsize)
print(plotExpression(pseudobulk_sce, features=pseudo_hvg[1:10], exprs_values = 'logcounts') + fontsize)
print(plotExpression(pseudobulk_sce, features=pseudo_hvg[1:10], exprs_values = 'normcounts') + fontsize)
print(plotExpression(pseudobulk_sce, features=pseudo_hvg[1:10], exprs_values = 'lowrank') + fontsize)
print(plotExpression(pseudobulk_sce, features=pseudo_hvg[1:10], exprs_values = 'pc_corrected') + fontsize)
dev.off()

# clear memory
rm(pseudo_pos_bio_sce, pseudo_pos_bio_counts_mat, min_expressed_pseudo_pos_bio_counts_mat)
gc(reset = T)


### sanity check
sanity_check_plot <- function(raw_fn, normalized_fn, denoised_fn, pcremoved_fn, plt_fn, n.genes = 10){
  df0 = read_feather_df(fn = raw_fn, rownames.col = 1)
  df1 = read_feather_df(fn = normalized_fn, rownames.col = 1)
  df2 = read_feather_df(fn = denoised_fn, rownames.col = 1)
  df3 = read_feather_df(fn = pcremoved_fn, rownames.col = 1)
  
  all_genes = intersect(rownames(df0), rownames(df1))
  all_genes = intersect(all_genes, rownames(df2))
  all_genes = intersect(all_genes, rownames(df3))
  selected_genes = sample(all_genes, size=n.genes)
  
  all_cells = intersect(colnames(df0), colnames(df1))
  all_cells = intersect(all_cells, colnames(df2))
  all_cells = intersect(all_cells, colnames(df3))
  df0 = as.matrix(df0[,all_cells])
  df1 = as.matrix(df1[,all_cells])
  df2 = as.matrix(df2[,all_cells])
  df3 = as.matrix(df3[,all_cells])
  
  
  pdf(plt_fn)
  for(g in selected_genes){
    plot(x=df0[g,], y=df1[g,], xlab = sprintf('Raw:%s',g), ylab = sprintf('Normalized:%s',g), main = "", col = rgb(1,0,0,0.3))
    plot(x=df0[g,], y=df2[g,], xlab = sprintf('Raw:%s',g), ylab = sprintf('Denoised:%s',g), main = "", col = rgb(0,1,0,0.3))
    plot(x=df0[g,], y=df3[g,], xlab = sprintf('Raw:%s',g), ylab = sprintf('PC-removed:%s',g), main = "", col = rgb(0,0,1,0.3))
  }
  dev.off()
  
  # clear memory
  rm(df0, df1, df2, df3)
  gc(reset = T)
  
  invisible(NULL)
}

cell_sanity_plot_fn = sprintf('%s_%d_cell_sanity_gene_expression.pdf', out_pfx, get_plot_count())
sanity_check_plot(raw_fn = count_fn, normalized_fn = hvg_normalized_fn, denoised_fn = hvg_denoised_fn, pcremoved_fn = hvg_pc_corrected_fn, plt_fn = cell_sanity_plot_fn, n.genes = 10)

pseudo_sanity_plot_fn = sprintf('%s_%d_pseudobulk_sanity_gene_expression.pdf', out_pfx, get_plot_count())
sanity_check_plot(raw_fn = pseudo_hvg_counts_fn, normalized_fn = pseudo_hvg_normalized_fn, denoised_fn = pseudo_hvg_denoised_fn, pcremoved_fn = pseudo_hvg_pc_corrected_fn, plt_fn = pseudo_sanity_plot_fn, n.genes = 10)

### save final sce file
final_sce_fn = sprintf('%s_final_sce.rds', out_pfx)
saveRDS(sce, file =  final_sce_fn)

### save final workspace
pseudobulk_image_fn = sprintf('%s_final_image_after_pseudobulk_processing.RData', out_pfx)
save.image(file = pseudobulk_image_fn)
