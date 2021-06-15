fun_convert_seurat_to_adata = function(seurat_obj, slot="data", assay="RNA", add_reduction=TRUE){
  require(reticulate)
  require(Seurat)
  sc <- import("scanpy")
  X = t(as.matrix(GetAssayData(seurat_obj, slot=slot, assay=assay)))
  obs = seurat_obj[[]]
  var = GetAssay(seurat_obj)[[]]
  if(ncol(var)==0){
    var=as.data.frame(colnames(X))
    colnames(var)='genes'
    rownames(var)=genes$genes
    }
  adata <- sc$AnnData( X  = X, obs = obs, var = var)
  if(add_reduction){
    adata$obsm$update(umap = Embeddings(seurat_obj, "umap"))
    adata$obsm$update(pca = Embeddings(seurat_obj, "pca"))
    if("wnn.umap" %in% names(seurat_obj@reductions)){
      adata$obsm$update(wnn.umap = Embeddings(seurat_obj, "wnn.umap"))
    }
    if("phate" %in% names(seurat_obj@reductions)){
      adata$obsm$update(phate = Embeddings(seurat_obj, "phate"))
    }
  }
  adata
}
