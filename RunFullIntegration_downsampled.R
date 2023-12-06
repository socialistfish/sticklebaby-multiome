library(Seurat)
library(Signac)
library(dyplr)

###########################################################################
# Script to run integration between samples and modalities WITH downsampling
# Aim: Generate reference data set to project full dataset onto
# 
###########################################################################

RunFullIntegration <- function (samples = all, data.list, seed=42, samplesize=1000){
    # Part 1: Create or set list of SeuratObjects named data.list
    if (!(missing(data.list))){
        data.list <- data.list
    }
    else {
        if (samples == all){
            samplelist <- c("sample01", "sample06", "sample09", "sample11","sample18", "sample20","sample03", "sample05", "sample13", "sample16", "sample21", "sample23","sample04", "sample07", "sample12", "sample15", "sample17", "sample24","sample02", "sample08", "sample10", "sample14", "sample19", "sample22")
        }
        else {samplelist = samples}
        samplelist <- samples
        data.list <- CreateDataList(samplelist=samplelist)
    }

    # Part 2: Downsample each SeuratObject in list to samplesize (default 1k) (note: might change this to # of cells in smallest object later)
    # The two loops are a workaround since using lapply or a single loop returned non-random sample - debug!
        y <- {}
        x <- {}
        subset <- {}
        for(sample in names(data.list)){
                x <- data.list[[sample]]
                y <- sample(Cells(x), size=samplesize, replace=F)
                subset[[sample]] <- y
        }
        rm(x,y)

        print(subset)
        for(sample in names(subset)){
                x <- data.list[[sample]]
                data.list[[sample]] <- x[,subset[[sample]]]
        }

        rm(x,y)

    message("Downsampling finished")
    print(data.list)
    # Part 3: Processing GEX
    message("Start processing for GEX")
    # Normalisation
    for(sample in names(data.list)){DefaultAssay(data.list[[sample]]) <- "RNA"}

    data.list <- lapply(data.list, SCTransform, return.only.var.genes=F) #Seurat 5
    message("GEX normalisation finished.")
    data.list <- lapply(data.list, RunPCA, npcs=100, return.only.var.genes = F, assay="SCT")
    d <- 100 # number of PCs to use downstream
    message("GEX PCA finished.")

    # Integration across samples for GEX
    for(sample in names(data.list)){DefaultAssay(data.list[[sample]]) <- "SCT"}

    gex.features <- SelectIntegrationFeatures(object.list = data.list, nfeatures = 3000)
    data.list <- PrepSCTIntegration(
        object.list = data.list,
        anchor.features = gex.features,
        verbose = FALSE
    )
gex.anchors <- FindIntegrationAnchors(
        object.list = data.list,
        normalization.method = "SCT",
        anchor.features = gex.features,
        reduction = "rpca"
    )
    data.integrated <- IntegrateData(
        anchorset = gex.anchors, dims = 1:d,
        features.to.integrate = gex.features,    #Reduce(intersect,lapply(data.list,rownames)),
        normalization.method = "SCT"
    )

    message("GEX integration finished")

    DefaultAssay(data.integrated) <- "integrated"

    data.integrated <- RunPCA(data.integrated, npcs=d)
    data.integrated <- FindNeighbors(data.integrated, reduction = "pca", dims = 1:d)

    data.integrated <- FindClusters(data.integrated, resolution = 1, algorithm = 1, reduction = "pca")
    #data.integrated <- BuildClusterTree(data.integrated, reorder = T, reorder.numeric = T, dims = 1:d)
    data.integrated <- RunUMAP(data.integrated, reduction = "pca", seed.use = 404, dims = 1:d, reduction.key = "UMAPgex_", reduction.name="UMAPgex")

    message("GEX processing finished")

    # Part 4: Processing ATAC
    message("Start processing for ATAC")
    # I will skip re-calling of my peaks since I already called my peaks on pseudobulk data using macs2
    for(sample in names(data.list)){DefaultAssay(data.list[[sample]]) <- "ATAC"}

    data.list <- lapply(data.list, FindTopFeatures, assay = "ATAC")
    data.list <- lapply(data.list, RunTFIDF, assay = "ATAC")
    data.list <- lapply(data.list, RunSVD)


    # Integrate across samples in ATAC 
    message("Dropping GEX assays from data.list!")
    for(sample in names(data.list)){data.list[['RNA']] <- NULL}
    for(sample in names(data.list)){data.list[['SCT']] <- NULL}

    atac.merged <- merge(data.list[[1]],y=data.list[2:length(data.list)])

    DefaultAssay(atac.merged) <- "ATAC"
    atac.merged <- FindTopFeatures(atac.merged)
    atac.merged <- RunTFIDF(atac.merged)
    atac.merged <- RunSVD(atac.merged)

    atac.anchors <- FindIntegrationAnchors(
        object.list = data.list,
        anchor.features=Reduce(intersect,lapply(data.list,rownames)),
        reduction = "rlsi",
        dims = 2:50
    )

    atac.integrated <- IntegrateEmbeddings(
        anchorset = atac.anchors,
        reductions = atac.merged[["lsi"]],
        new.reduction.name = "integrated_lsi",
        dims.to.integrate = 1:50
    )

    data.integrated@reductions$integrated_lsi <- atac.integrated@reductions$integrated_lsi
    message("ATAC integration finished")

    rm(atac.integrated, atac.merged)
    DefaultAssay(data.integrated) <- "ATAC"
    data.integrated <- RunUMAP(data.integrated, reduction = "integrated_lsi", seed.use=404, dims = 2:50, reduction.key = "UMAPatac_", reduction.name="UMAPatac")
    data.integrated <- FindNeighbors(data.integrated, reduction = "integrated_lsi", dims = 1:50)
    data.integrated <- FindClusters(data.integrated, resolution = 1, algorithm=1, reduction='integrated_lsi')
    #data.integrated <- BuildClusterTree(object = data.integrated, reorder = TRUE, reduction='integrated_lsi',reorder.numeric = TRUE, dims = 1:50)
    data.integrated$atac_clusters <- data.integrated$seurat_clusters
    #save(data.integrated, file="dataintegrated.eachassay_v4.RObj")
    message("ATAC processing finished")
    rm(atac.integrated, atac.merged)


    # Part 5: Integration across modalities
    message("Start multi-modal integration")
    data.integrated <- FindMultiModalNeighbors(
        object = data.integrated,
        reduction.list = list("pca", "integrated_lsi"),
        dims.list = list(1:100, 2:40),
        #modality.weight.name = "RNA.weight",
        verbose = TRUE
    )
    data.integrated <- RunUMAP(
    object = data.integrated,
    nn.name = "weighted.nn",
    assay = "integrated",
    seed.use = 404,
    dims = 1:30,
    verbose = TRUE,
    reduction.key="umapwnn",
        reduction.name="UMAPwnn"
    )

    data.integrated <- FindClusters(
    object=data.integrated,
    resolution = 1.0,
    random.seed= 404,
    graph.name = "wsnn"
    )

   # data.integrated <- BuildClusterTree(
    #object = data.integrated, 
    #reorder = TRUE, 
    #reorder.numeric = TRUE, 
    #reduction="wsnn",
    #dims = 1:d
    #)
    data.integrated$wnn_clusters <- data.integrated$seurat_clusters

    #save(data.integrated, file="dataintegrated.all_v4.RObj")

    return(data.integrated)
}
