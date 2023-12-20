# ------------------------------------------------------------------------------
# NeuroConnectivity - Profiling
#
# Author: Marlies Verschuuren
# Modified by: Marlies Verschuuren
# Creation date: 2023-12-01
# Last Modified: 2023-12-20
# ------------------------------------------------------------------------------

#--1. User Settings-------------------------------------------------------------
#----1.1. Select Directories----------------------------------------------------
dir.morph="/Users/marliesverschuuren/Library/CloudStorage/OneDrive-UniversiteitAntwerpen/Projects/DeVosLab/NeuroConnectivity/Results_PLA/Neuro_Morph_MergedData"
dir.func="/Users/marliesverschuuren/Library/CloudStorage/OneDrive-UniversiteitAntwerpen/Projects/DeVosLab/NeuroConnectivity/Results_PLA/Neuro_Func_MergedData"

dir.output="/Users/marliesverschuuren/Library/CloudStorage/OneDrive-UniversiteitAntwerpen/Projects/DeVosLab/NeuroConnectivity/Results_PLA"

#--2. Packages And Settings-----------------------------------------------------
#----2.1. Packages--------------------------------------------------------------
if (!require("tidyverse")) {install.packages("tidyverse"); require("tidyverse")}
if (!require("data.table")) {install.packages("data.table"); require("data.table")}
if(any(grepl("package:plyr", search()))) detach("package:plyr") else message("plyr not loaded") #Detach Plyr packages to avoid problems with dplyr package
if (!require("gridExtra")) {install.packages("gridExtra"); require("gridExtra")}                #Function: marrangegrob
if (!require("pheatmap")) {install.packages("pheatmap"); require("pheatmap")}       
if (!require("factoextra")) {install.packages("factoextra"); require("factoextra")}      

#----2.2. Create Result Directories---------------------------------------------
dir.plot=file.path(dir.output, "Neuro_Profile_Plots")
dir.create(dir.plot, showWarnings = FALSE)

dir.mergedData=file.path(dir.output, "Neuro_Profile_MergedData")
dir.create(dir.mergedData, showWarnings = FALSE)

#--3. Read Norm Data ----------------------------------------------------------
data.func.well.norm=read.table(file=file.path(dir.func,"Data_Func_Well_Norm.txt"),sep="\t",dec=".", header=TRUE)
data.morph.well.norm=read.table(file=file.path(dir.morph,"Data_Morph_Well_Norm.txt"),sep="\t",dec=".", header=TRUE)
data.connectivity=left_join(data.func.well.norm,data.morph.well.norm)

#--4. Hierarchical clustering --------------------------------------------------
save_pheatmap_png <- function(x, filename, width=30, height=30, res = 300) {
  png(filename, width = width, height = height, res = res, units = "cm")
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

#----4.1. Morph ----------------------------------------------------------------
data.hc=data.morph.well.norm
var.features=names(data.hc)[-c(1:7)]
var.group=names(data.hc)[c(1:7)]

#Prepare data matrix
data.matrix=as.matrix(data.hc[,c(var.features)])
data.matrix = scale(data.matrix)
rownames(data.matrix)=seq(1:nrow(data.matrix))
data.hc$ID=str_replace(data.hc$ID,":","_")
data.hc$ID=str_replace(data.hc$ID,"-","_")
rowannotation=data.frame(REP=data.hc[,c("Rep")],ID=data.hc[,c("ID")])
colannotation=data.frame(ASSAY=ifelse(grepl(var.features, pattern="_SC"),"Morph","Func"))
rownames(colannotation)=var.features

#Adapt IDs, Reps and colors
unique(data.hc$ID)
col = list(
  ID = c(B27_NA_NA = "grey30", B27_AO_NA_NA = "#3ca17e", B27_1_10_AO_NA_NA = "#6fa8dc"),
  REP = c(Rep1 = "#be8989", Rep2 = "#8b5656", Rep3="#573636"),
  ASSAY = c(Morph = "grey10", Func = "grey90")
)

#Color heatmap
paletteLength=100
colHeatmap=colorRampPalette(c("#0b2f37","#134f5c","#89a7ad","white","#f0bd87","#e69138","#e66e38"))(paletteLength)
breaksHeatmap=c(seq(min(data.matrix, na.rm=TRUE), 0, length.out=ceiling(paletteLength/2) + 1),
                seq(max(data.matrix, na.rm=TRUE)/paletteLength, max(data.matrix, na.rm=TRUE), length.out=floor(paletteLength/2)))

#Heatmap
p1=pheatmap(data.matrix, annotation_row = rowannotation, annotation_col = colannotation,
            clustering_distance_rows = "euclidean" , #"euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
            clustering_distance_cols = "euclidean" , #"euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
            clustering_method = 'ward.D', #"ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median" or "centroid".
            annotation_colors = col, color=colHeatmap, breaks=breaksHeatmap, border_color="white", na_col = "grey50",
            show_colnames = TRUE, show_rownames = FALSE,
            cutree_rows = length(unique(data.hc$ID)))

save_pheatmap_png(p1, file.path(dir.plot,"HC_Morph.jpg"))

#----4.2. Func ----------------------------------------------------------------
data.hc=data.func.well.norm
var.features=names(data.hc)[-c(1:7)]
var.group=names(data.hc)[c(1:7)]

#Prepare data matrix
data.matrix=as.matrix(data.hc[,c(var.features)])
data.matrix = scale(data.matrix)
rownames(data.matrix)=seq(1:nrow(data.matrix))
data.hc$ID=str_replace(data.hc$ID,":","_")
data.hc$ID=str_replace(data.hc$ID,"-","_")
rowannotation=data.frame(REP=data.hc[,c("Rep")],ID=data.hc[,c("ID")])
colannotation=data.frame(ASSAY=ifelse(grepl(var.features, pattern="_SC"),"Morph","Func"))
rownames(colannotation)=var.features


#Adapt IDs, Reps and colors
unique(data.hc$ID)
col = list(
  ID = c(B27_NA_NA = "grey30", B27_AO_NA_NA = "#3ca17e", B27_1_10_AO_NA_NA = "#6fa8dc"),
  REP = c(Rep1 = "#be8989", Rep2 = "#8b5656", Rep3="#573636"),
  ASSAY = c(Morph = "grey10", Func = "grey90")
)

#Color heatmap
paletteLength=100
colHeatmap=colorRampPalette(c("#0b2f37","#134f5c","#89a7ad","white","#f0bd87","#e69138","#e66e38"))(paletteLength)
breaksHeatmap=c(seq(min(data.matrix, na.rm=TRUE), 0, length.out=ceiling(paletteLength/2) + 1),
                seq(max(data.matrix, na.rm=TRUE)/paletteLength, max(data.matrix, na.rm=TRUE), length.out=floor(paletteLength/2)))

#Heatmap
p2=pheatmap(data.matrix, annotation_row = rowannotation, annotation_col = colannotation, 
            clustering_distance_rows = "euclidean" , #"euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
            clustering_distance_cols = "euclidean" , #"euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
            clustering_method = 'ward.D', #"ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median" or "centroid".
            annotation_colors = col, color=colHeatmap, breaks=breaksHeatmap, border_color="white", na_col = "grey50",
            show_colnames = TRUE, show_rownames = FALSE,
            cutree_rows = length(unique(data.hc$ID)))
save_pheatmap_png(p2, file.path(dir.plot,"HC_Func.jpg"))


#----4.3. Combined ----------------------------------------------------------------
data.hc=data.connectivity
var.features=names(data.hc)[-c(1:7)]
var.group=names(data.hc)[c(1:7)]

#Prepare data matrix
data.matrix=as.matrix(data.hc[,c(var.features)])
data.matrix = scale(data.matrix)
rownames(data.matrix)=seq(1:nrow(data.matrix))
data.hc$ID=str_replace(data.hc$ID,":","_")
data.hc$ID=str_replace(data.hc$ID,"-","_")
rowannotation=data.frame(REP=data.hc[,c("Rep")],ID=data.hc[,c("ID")])
colannotation=data.frame(ASSAY=ifelse(grepl(var.features, pattern="_SC"),"Morph","Func"))
rownames(colannotation)=var.features

#Adapt IDs, Reps and colors
unique(data.hc$ID)
col = list(
  ID = c(B27_NA_NA = "grey30", B27_AO_NA_NA = "#3ca17e", B27_1_10_AO_NA_NA = "#6fa8dc"),
  REP = c(Rep1 = "#be8989", Rep2 = "#8b5656", Rep3="#573636"),
  ASSAY = c(Morph = "grey10", Func = "grey90")
)

#Color heatmap
paletteLength=100
colHeatmap=colorRampPalette(c("#0b2f37","#134f5c","#89a7ad","white","#f0bd87","#e69138","#e66e38"))(paletteLength)
breaksHeatmap=c(seq(min(data.matrix, na.rm=TRUE), 0, length.out=ceiling(paletteLength/2) + 1),
                seq(max(data.matrix, na.rm=TRUE)/paletteLength, max(data.matrix, na.rm=TRUE), length.out=floor(paletteLength/2)))

#Heatmap
p3=pheatmap(data.matrix, annotation_row = rowannotation, annotation_col = colannotation,
            clustering_distance_rows = "euclidean" , #"euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
            clustering_distance_cols = "euclidean" , #"euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
            clustering_method = 'ward.D', #"ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median" or "centroid".
            annotation_colors = col, color=colHeatmap, breaks=breaksHeatmap, border_color="white", na_col = "grey50",
            show_colnames = TRUE, show_rownames = FALSE,
            cutree_rows = length(unique(data.hc$ID)))

save_pheatmap_png(p3, file.path(dir.plot,"HC_Combined.jpg"))

#--5. PCA ----------------------------------------------------------------------
#Adapt colors
col=c("#6fa8dc","grey10","#3ca17e")

theme_pca = theme_minimal(base_size=8)+
  theme(panel.grid = element_blank(),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        legend.position="bottom",
        legend.key.size = unit(0.2, 'cm'), 
        legend.key.height = unit(0.2, 'cm'), 
        legend.key.width = unit(0.2, 'cm'),
        legend.title = element_text(size=4),
        legend.text = element_text(size=(4)))

#----5.1. Morph ----------------------------------------------------------------
data.pca=data.morph.well.norm
var.features=names(data.pca)[-c(1:7)]
var.group=names(data.pca)[c(1:7)]

#Remove samples with NA
data.pca=data.pca[complete.cases(data.pca[,var.features]),]

pca = data.pca[,c(var.features)] %>%
  prcomp(retx=T, scale=T)

data.rot=as.data.frame(pca$x) %>%
  cbind(data.pca[,var.group])

pctPC1=round((pca$sdev^2/sum(pca$sdev^2)*100)[1],2)
pctPC2=round((pca$sdev^2/sum(pca$sdev^2)*100)[2],2)

data.plot=data.rot
p1=ggplot(data=data.plot,aes(x=PC1,y=PC2, fill=ID, color=ID)) +
  stat_ellipse(level=0.68, geom="polygon", alpha=0.7, color=NA)+
  stat_ellipse(level=0.68, geom="polygon", fill=NA)+
  geom_point(aes(shape=Rep),size=1, alpha=0.7, color="black")+
  scale_color_manual(values=col,guide=guide_legend(nrow=2)) +
  scale_fill_manual(values=col,guide=guide_legend(nrow=2)) +
  scale_shape_manual(values=seq(21,21+length(unique(data.plot$Rep)),1))+
  xlab(paste("PC1 ( ", paste(pctPC1,"% explained var.)",sep=""),sep=""))+
  ylab(paste("PC2 ( ", paste(pctPC2,"% explained var.)",sep=""),sep=""))+
  ggtitle("Morph")+
  theme_pca

p1

ggsave(file=file.path(dir.plot,"PCA_Morph.jpg"),width = 7, height=7, units = "cm")
write.table(x=data.rot,file=file.path(dir.mergedData,"Data_Profile_PCA_Morph_Scores.txt"), sep='\t', dec=".", col.names = TRUE, row.names = FALSE)
write.table(x=pca$rotation,file=file.path(dir.mergedData,"Data_Profile_PCA_Morph_Rotation.txt"), sep='\t', dec=".", col.names = TRUE, row.names = TRUE)

#----5.2. Func ----------------------------------------------------------------
data.pca=data.func.well.norm
var.features=names(data.pca)[-c(1:7)]
var.group=names(data.pca)[c(1:7)]

#Remove samples with NA
data.pca=data.pca[complete.cases(data.pca[,var.features]),]

pca = data.pca[,c(var.features)] %>%
  prcomp(retx=T, scale=T)

data.rot=as.data.frame(pca$x) %>%
  cbind(data.pca[,var.group])

pctPC1=round((pca$sdev^2/sum(pca$sdev^2)*100)[1],2)
pctPC2=round((pca$sdev^2/sum(pca$sdev^2)*100)[2],2)

data.plot=data.rot
p2=ggplot(data=data.plot,aes(x=PC1,y=PC2, fill=ID, color=ID)) +
  stat_ellipse(level=0.68, geom="polygon", alpha=0.7, color=NA)+
  stat_ellipse(level=0.68, geom="polygon", fill=NA)+
  geom_point(aes(shape=Rep),size=1, alpha=0.7, color="black")+
  scale_color_manual(values=col,guide=guide_legend(nrow=2)) +
  scale_fill_manual(values=col,guide=guide_legend(nrow=2)) +
  scale_shape_manual(values=seq(21,21+length(unique(data.plot$Rep)),1))+
  xlab(paste("PC1 ( ", paste(pctPC1,"% explained var.)",sep=""),sep=""))+
  ylab(paste("PC2 ( ", paste(pctPC2,"% explained var.)",sep=""),sep=""))+
  ggtitle("Func")+
  theme_pca

p2
ggsave(file=file.path(dir.plot,"PCA_Func.jpg"),width = 7, height=7, units = "cm")
write.table(x=data.rot,file=file.path(dir.mergedData,"Data_Profile_PCA_Func_Scores.txt"), sep='\t', dec=".", col.names = TRUE, row.names = FALSE)
write.table(x=pca$rotation,file=file.path(dir.mergedData,"Data_Profile_PCA_Func_Rotation.txt"), sep='\t', dec=".", col.names = TRUE, row.names = TRUE)


#----5.3. Combined----------------------------------------------------------------
data.pca=data.connectivity
var.features=names(data.pca)[-c(1:7)]
var.group=names(data.pca)[c(1:7)]

#Remove samples with NA
data.pca=data.pca[complete.cases(data.pca[,var.features]),]

pca = data.pca[,c(var.features)] %>%
  prcomp(retx=T, scale=T)

data.rot=as.data.frame(pca$x) %>%
  cbind(data.pca[,var.group])

pctPC1=round((pca$sdev^2/sum(pca$sdev^2)*100)[1],2)
pctPC2=round((pca$sdev^2/sum(pca$sdev^2)*100)[2],2)

data.plot=data.rot

p3=ggplot(data=data.plot,aes(x=PC1,y=PC2, fill=ID, color=ID)) +
  stat_ellipse(level=0.68, geom="polygon", alpha=0.7, color=NA)+
  stat_ellipse(level=0.68, geom="polygon", fill=NA)+
  geom_point(aes(shape=Rep),size=1, alpha=0.7, color="black")+
  scale_color_manual(values=col,guide=guide_legend(nrow=2)) +
  scale_fill_manual(values=col,guide=guide_legend(nrow=2)) +
  scale_shape_manual(values=seq(21,21+length(unique(data.plot$Rep)),1))+
  xlab(paste("PC1 ( ", paste(pctPC1,"% explained var.)",sep=""),sep=""))+
  ylab(paste("PC2 ( ", paste(pctPC2,"% explained var.)",sep=""),sep=""))+
  ggtitle("Combined")+
  theme_pca
p3
ggsave(file=file.path(dir.plot,"PCA_Combined.jpg"),width = 7, height=7, units = "cm")

write.table(x=data.rot,file=file.path(dir.mergedData,"Data_Profile_PCA_Combined_Scores.txt"), sep='\t', dec=".", col.names = TRUE, row.names = FALSE)
write.table(x=pca$rotation,file=file.path(dir.mergedData,"Data_Profile_PCA_Combined_Rotation.txt"), sep='\t', dec=".", col.names = TRUE, row.names = TRUE)

grid.arrange(grobs = list(p1,p2,p3), ncol=3)
ggsave(file=file.path(dir.plot,"PCA.jpg"),grid.arrange(grobs = list(p1,p2,p3), ncol=3),width = 20, height=7, units = "cm")


