library("dplyr"); library("Seurat")
statistics_folder="dimredclean_statistics"

files<-list.files(statistics_folder, full.names=TRUE)
# ignore statistics to do with tracks or signal intensity inside a surface.
ignore<-grepl("*_Track*|*Intensity*|*Center*|*Angle*|*Image*", files) 
surface_files<-files[!ignore]

readData<-function(file) {
    x <- read.csv(file, skip=2, header=TRUE)
    if ("Collection" %in% colnames(x)) { # remove irrelevant columns in the statistics file. e.g. trackID is not useful.
        x<-x %>% dplyr::select(., -c(Collection) )
    }
    if ("Time" %in% colnames(x)) {
        x<-x %>% dplyr::select(., -c(Time) )
    }
    if ("Image" %in% colnames(x)) {
        x<-x %>% dplyr::select(., -c(Image) )
    }
    if ("Channel" %in% colnames(x)) {
        x<-x %>% dplyr::select(., -c(Channel) )
    }
    if ("TrackID" %in% colnames(x)) {
        x<-x %>% dplyr::select(., -c(TrackID) )
    }
    if ("Category" %in% colnames(x)) {
        x<-x %>% dplyr::select(., -c(Category) )
    }
    if ("Unit" %in% colnames(x)) {
        x<-x %>% dplyr::select(., -c(Unit) )
    }
    if ("X" %in% colnames(x)) {
        x<-x %>% dplyr::select(., -c(X) )
    }
    x <- x
}

data<-lapply(surface_files, readData)

# get most common numrows in the data. This is a shortcut to get rid of stats that aren't for individual cells. 
nrows<-lapply(data, nrow) %>% unlist() %>% table() %>% which.max () %>% names() %>% as.numeric()

is_cell_data<-lapply(data, function(x) {
    nrow(x)==nrows
})

cell_data<-data[unlist(is_cell_data)]

# check order is the same for each object
id_order<-cell_data[[1]]$ID
order_df<-lapply(cell_data, function(x) {
    x$ID == id_order
})

# order is same, so remove "ID" column from each as it is not neeeded.
cell_data<-lapply(cell_data, function(x) {
    x<-x %>% dplyr::select(., -ID)
})

cell_matrix<-bind_cols(cell_data) %>% t()
colnames(cell_matrix) <- id_order
# convert to matrix, with last column (ID) as cell name, so colnames, and rownames as the colnames of each df,

# Normalization function.
norm <- function(x){
    x<-(x-min(x))/(max(x)-min(x))
}

# Final cleaning of data.
cell_matrix <- t(apply(cell_matrix, 1, FUN=norm))
cell_matrix <- na.omit(cell_matrix)

# Import into Seurat and run standard dimension reduction workflow.
seuobj<-CreateSeuratObject(cell_matrix, project = "SeuratProject", assay = "RNA",
                   min.cells = 0, min.features = 0, names.field = 1,
                meta.data = NULL)

seuobj<-NormalizeData(seuobj)
seuobj<-ScaleData(seuobj)

seuobj<-FindVariableFeatures(seuobj)
seuobj<-RunPCA(seuobj, pcs=5)
seuobj<-RunUMAP(seuobj, dims = 1:5)
seuobj<-FindNeighbors(seuobj, dims=1:5)
seuobj<-FindClusters(seuobj)

marker_genes<-FindAllMarkers(seuobj)
top_markers <- marker_genes$gene %>% head(12)

# Feature Plots
DimPlot(seuobj)
FeaturePlot(seuobj, features=c(seuobj@assays$RNA@var.features), slot="counts")
FeaturePlot(seuobj, features=c("Volume", "Area", "Ellipticity..oblate.", "Ellipticity..prolate.","Number.of.Voxels"), slot="scale.data")
DimHeatmap(seuobj, dims = 1:5, cells = 243, balanced = TRUE)

# Extract tiles from cluster X
cell_tiles <- row.names(seuobj@meta.data[seuobj@meta.data$seurat_clusters == 3, ])




