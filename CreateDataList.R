CreateDataList <- function (samplelist, path="/fml/obi/projects/PJ064_scMultiome_Environment/cellranger-runs/transdecoder_", GEX_lower=500, GEX_upper=10000, ATAC_lower=1000, ATAC_upper=100000. nucleosome_signal=2, TSS.enrichment=2){

    gffRangedData <- import.gff("/fml/chones/genome/gbdb/gasAcu1/Gasterosteus_aculeatus.BROADS1.104.gasAcu1.lifted.nosort.gtf")
    annotation<-as(gffRangedData, "GRanges")
    genome(annotation) <- "gasAcu1"

    data.list <- list()

    for (sample in samplelist){
        print(paste0("Adding ", sample," to data.list"))

        rawdata <- Read10X_h5(paste0("/fml/obi/projects/PJ064_scMultiome_Environment/cellranger-runs/transdecoder_",sample,"/outs/filtered_feature_bc_matrix.h5"))
        rna_counts <- rawdata[["Gene Expression"]]
        atac_counts <- rawdata[["Peaks"]]
        x <- CreateSeuratObject(counts = rna_counts)
        frag.file <- paste0("../",sample,"/outs/atac_fragments.tsv.gz")
        chrom_assay <- CreateChromatinAssay(counts=atac_counts, sep = c(":", "-"), fragments=frag.file, annotation=annotation, genome="gasAcu1")
        x[["ATAC"]] <- chrom_assay

        DefaultAssay(x) <- "ATAC"
        x <- NucleosomeSignal(x)
        x <- TSSEnrichment(x)

        x <- subset(
        x = x,
        subset = nCount_RNA > GEX_lower &
        nCount_RNA < GEX_upper &
        nCount_ATAC < ATAC_upper &
        nCount_ATAC > ATAC_lower &
        nucleosome_signal < nucleosome_signal &
        TSS.enrichment > TSS.enrichment
        )

        x@meta.data$orig.ident <- sample
        x <- RenameCells(x, add.cell.id = sample)

        data.list[[paste0(sample)]] <- x
        rm(x, chrom_assay, frag.file, rawdata, rna_counts, atac_counts)
        print(paste0("Finished adding ", sample," to data.list"))
    }

    return(data.list)
}
