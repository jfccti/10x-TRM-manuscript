jianing_dat=readRDS("/Users/aleksandar/Dropbox/jianing_R01_documents/integrated8.rpca0919.rds")

##run VIPER
##load aracne networks
nets=readRDS("~/Documents/Documents/MED SCHOOL/summer 2018 (drake:califano lab)/Kidney:Bladder:Prostate:GBM/itherapy/kidney_prostate_bladder_gbm/itherapy_pruned.rds")
dat=jianing_dat@assays$integrated@scale.data
vp_meta_1<-viper(dat[,1:5000], nets, method = 'none')
dat=dat[,5001:ncol(dat)]
vp_meta_2<-viper(dat[,1:5000], nets, method = 'none')
dat=dat[,5001:ncol(dat)]
vp_meta_3<-viper(dat[,1:5000], nets, method = 'none')
dat=dat[,5001:ncol(dat)]
vp_meta_4<-viper(dat[,1:5000], nets, method = 'none')
dat=dat[,5001:ncol(dat)]
vp_meta_5<-viper(dat[,1:5000], nets, method = 'none')
dat=dat[,5001:ncol(dat)]
vp_meta_6<-viper(dat[,1:5000], nets, method = 'none')
dat=dat[,5001:ncol(dat)]
vp_meta_7<-viper(dat[,1:5000], nets, method = 'none')
dat=dat[,5001:ncol(dat)]
vp_meta_8<-viper(dat[,1:5000], nets, method = 'none')
dat=dat[,5001:ncol(dat)]
vp_meta_9<-viper(dat[,1:5000], nets, method = 'none')
dat=dat[,5001:ncol(dat)]
vp_meta_10<-viper(dat[,1:5000], nets, method = 'none')
dat=dat[,5001:ncol(dat)]
vp_meta_11<-viper(dat, nets, method = 'none')
vp_meta=cbind(vp_meta_1,vp_meta_2,vp_meta_3,vp_meta_4,vp_meta_5,vp_meta_6,vp_meta_7,vp_meta_8,vp_meta_9,vp_meta_10,vp_meta_11)
rm(vp_meta_1,vp_meta_2,vp_meta_3,vp_meta_4,vp_meta_5,vp_meta_6,vp_meta_7,vp_meta_8,vp_meta_9,vp_meta_10,vp_meta_11)
##VIPER clustering
cbcMRs <- CBCMRs(vp_meta)
jianing.integrated.meta.vp <- CreateSeuratObject(counts = vp_meta[cbcMRs,])
jianing.integrated.meta.vp@assays$RNA@scale.data=as.matrix(jianing.integrated.meta.vp@assays$RNA@data)
jianing.integrated.meta.vp <- RunPCA(jianing.integrated.meta.vp,features=rownames(jianing.integrated.meta.vp))
jianing.integrated.meta.vp <- RunUMAP(jianing.integrated.meta.vp, dims = 1:50, verbose = FALSE,metric="correlation")
jianing.integrated.meta.vp <- FindNeighbors(jianing.integrated.meta.vp, dims = 1:50, verbose = FALSE)
jianing.integrated.meta.vp <- FindClusters(jianing.integrated.meta.vp, resolution=seq(0.01,1,by=0.01), verbose = FALSE,algorithm=1) 
clust=jianing.integrated.meta.vp@meta.data[,which(grepl("RNA_snn_res.",colnames(jianing.integrated.meta.vp@meta.data)))]
mat=jianing.integrated.meta.vp@assays$RNA@scale.data
out=sil_subsample(mat,clust)
means=out[[1]]
sd=out[[2]]
x=seq(0.01,1,by=0.01)
errbar(x,means,means+sd,means-sd,ylab="mean silhouette score",xlab="resolution parameter")
lines(x,means)
best=tail(x[which(means[5:length(means)]==max(means[5:length(means)]))+4],n=1)
legend("topright",paste("Best",best,sep = " = "))
jianing.integrated.meta.vp$seurat_clusters=jianing.integrated.meta.vp@meta.data[,which(colnames(jianing.integrated.meta.vp@meta.data)==paste("RNA_snn_res.",best,sep=""))]
Idents(jianing.integrated.meta.vp) <- "seurat_clusters"
plot(DimPlot(jianing.integrated.meta.vp, reduction = "umap",label = TRUE,label.size=7,repel=T) + NoLegend())
markers.vp <- FindAllMarkers(jianing.integrated.meta.vp, only.pos = TRUE, min.pct = 0, logfc.threshold = 0,test.use = "t")
top10 <- markers.vp %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
geneHeatmap_plot(jianing.integrated.meta.vp@assays$RNA@counts,jianing.integrated.meta.vp$seurat_clusters,top10$gene,n_top_genes_per_cluster = 5,scaled = T)
l=jianing.integrated.meta.vp$blueprint_labels
l[which(jianing.integrated.meta.vp$blueprint_pvals>0.1)]=NA
l[which(l %in% names(which(table(l)<150)))]=NA
jianing.integrated.meta.vp$l=l
Idents(jianing.integrated.meta.vp) <- "l"
plot(DimPlot(jianing.integrated.meta.vp, reduction = "umap",label=TRUE,repel=T,label.size=5)+NoLegend())
saveRDS(jianing.integrated.meta.vp, file = "/Users/aleksandar/Dropbox/jianing_R01_documents/jianing.integrated.meta.vp.rds")


jianing.integrated.meta.vp$categories=jianing_dat$categories[colnames(jianing.integrated.meta.vp)]
plot(DimPlot(jianing.integrated.meta.vp, reduction = "umap",label=TRUE,repel=T,label.size=5,group.by="categories")+NoLegend())
plot(DimPlot(jianing.integrated.meta.vp[,jianing_dat$categories != "Others"], reduction = "umap",label=TRUE,repel=T,label.size=5,group.by="categories"))

table(jianing_dat$repos)

jianing.integrated.meta.vp$gene_exp_clusters=jianing_dat$integrated_snn_res.0.5[colnames(jianing.integrated.meta.vp)]


#non transplant control D251
#queiscent : Pt15_POD1194, Pt13_POD1032_IEL,Pt13_POD1032_LPL, Pt16_POD1004_IEL, Pt16_POD1004_LPL, Pt21_POD626  vs rejection: all others 

table(jianing_dat$categories)


table(jianing_dat$integrated_snn_res.0.5)


library(Seurat)
devtools::install_github('satijalab/seurat-data')


library(SeuratData)
library(ggplot2)
library(patchwork)

DimPlot(jianing.integrated.meta.vp, reduction = "umap")
DoHeatmap(jianing.integrated.meta.vp, features= c("BACH2", "TRAK2", "MXI1", "NFKBIZ", "ZC3H12D", "SETD7", "MORC3","KLF9", "NSD3", "STAT5B", "EOMES", "TBX21", "RUNX3","IFNG","MYBL1", "ZEB2", "ZSCAN22", "FOXP4", "NFATC2", "SMAD5", "SMAD7","LMO4", "NR4A2","NR4A3", "ETV3", "MAFF",  "ZFP36L1", "ID2", "MSC","ELF4", "NR3C1", "STAT4", "MXD1", "EGR1", "EGR3","NR4A1", "BTG2", "PPP1R10","MED29","CSRNP1", "LBH", "NSMCE3", "POU3F1", "FOSB","JUND","PHF1","ID3", "EAF1",  "IRF7", "POLR3G", "ZC3H7B", "FOSL1", "ZNF607", "MAMLD1"), size=4)


