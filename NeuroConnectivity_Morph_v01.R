# ------------------------------------------------------------------------------
# NeuroConnectivity - Morphological analysis
#
# Author: Marlies Verschuuren
# Modified by: Marlies Verschuuren
# Creation date: 2022-01-13
# Last Modified: 2023-12-20
# ------------------------------------------------------------------------------

#--1. User Settings-------------------------------------------------------------
#----1.1. Select Directories----------------------------------------------------
# Input: "Data" folder with structure: Rep > Plate > Morph > Output
#                                      Rep > Plate > PlateLayout.txt
dir.input="/Users/marliesverschuuren/Documents/UA_DataSets/NeuroConnectivity/PLA/Data"
dir.output="/Users/marliesverschuuren/Library/CloudStorage/OneDrive-UniversiteitAntwerpen/Projects/DeVosLab/NeuroConnectivity/Results_PLA"

#----1.2. Control condition---------------------------------------------------
ctrl.condition="B27_NA_NA" #Condition_Treatment_Concentration

#--2. Packages And Settings-----------------------------------------------------
#----2.1. Packages--------------------------------------------------------------
if (!require("tidyverse")) {install.packages("tidyverse"); require("tidyverse")}
if (!require("data.table")) {install.packages("data.table"); require("data.table")}
if (!require("RColorBrewer")) {install.packages("RColorBrewer"); require("RColorBrewer")}       #Function: brewer.pal
if(any(grepl("package:plyr", search()))) detach("package:plyr") else message("plyr not loaded") #Detach Plyr packages to avoid problems with dplyr package
if (!require("gridExtra")) {install.packages("gridExtra"); require("gridExtra")}                #Function: marrangegrob

#----2.2. Create Result Directories---------------------------------------------
dir.plot=file.path(dir.output, "Neuro_Morph_Plots")
dir.create(dir.plot, showWarnings = FALSE)

dir.mergedData=file.path(dir.output, "Neuro_Morph_MergedData")
dir.create(dir.mergedData, showWarnings = FALSE)

#--3. Read Morph Data ----------------------------------------------------------
#----3.1. Read data from summary files -----------------------------------------
data = data.frame()
folders.rep=list.files(path=file.path(dir.input))
data.list=list()
i=1
#Loop over all rep folders 
for (rep in folders.rep){ 
  folders.plate=list.files(path=file.path(dir.input,rep))
  #Loop over all Plate folders in Rep folder
  for (plate in folders.plate){
    data.plate = data.frame()
    files = list.files(path = file.path(dir.input,rep,plate,"Morph","Output"), pattern = "summary", full.names = TRUE)
    
    expr="well [A-Z][0-9]{1,2}" #### ADAPT IF NEEDED
    #Check regexpr:
    #File=files[1]
    #str_extract(File,expr)
    
    data.plate = tibble(File = files) %>%
      mutate(data = lapply(File, fread)) %>%
      unnest(data)%>%
      mutate(Rep=rep,
             Plate=plate,
             File=File,
             Image=substring(str_extract(File,"Output/.+"),8,nchar(str_extract(File,"Output/.+"))-12),
             Well=substring(str_extract(File,expr),6,9)) #### ADAPT IF NEEDED
    
    #Add 0 if well is defined in format B2
    data.plate$Well=ifelse(nchar(data.plate$Well)==2,paste(substring(data.plate$Well,1,1),"0",substring(data.plate$Well,2,2),sep=""),data.plate$Well)

    #Merge with plate layout
    layout=read.table(file.path(dir.input,rep,plate,"PlateLayout.txt"),header=TRUE, sep='\t', fill=TRUE)
    data.plate=data.plate%>%
      left_join(layout)
    
    data.list[[i]]=data.plate
    i=i+1
  }
}

#Bind all replicates and plates
data.morph.image=rbindlist(data.list)%>%
  as.data.frame()

#Check number of images
check=data.morph.image%>%
  group_by(Rep,Plate,Condition,Treatment,Concentration,Well)%>%
  summarise(ImageCount=n_distinct(File))

#----3.2. Define variables  ----------------------------------------------------
data.morph.image$ID=paste(data.morph.image$Condition,data.morph.image$Treatment, data.morph.image$Concentration, sep="_")

#Spot density
if(any(grepl(pattern = "Neurites_",names(data.morph.image)))){
  var.spot=names(data.morph.image)[grepl(pattern = "Spot_SC",names(data.morph.image))]
  spot_channel=substring(unique(str_extract(var.spot, pattern = "Spot_SC[0-9]")),8,9)
  for(c in spot_channel){
    var.spot=paste("Spot_SC",c,"_Nr",sep="")
    var.new=paste("Spot_SC",c,"_SpotDensity",sep="")
    var.neurite=names(data.morph.image)[grepl(x = names(data.morph.image), pattern="Neurites_SC[0-9]_Area_MC[0-9]")]
    
    data.morph.image[,var.new] = data.morph.image[,var.spot]/data.morph.image[,var.neurite[1]]
  }
} else if (any(grepl(pattern = "Search_",names(data.morph.image)))){
  var.spot=names(data.morph.image)[grepl(pattern = "Spot_SC",names(data.morph.image))]
  spot_channel=substring(unique(str_extract(var.spot, pattern = "Spot_SC[0-9]")),8,9)
  for(c in spot_channel){
    var.spot=paste("Spot_SC",c,"_Nr",sep="")
    var.new=paste("Spot_SC",c,"_SpotDensity",sep="")
    var.neurite=names(data.morph.image)[grepl(x = names(data.morph.image), pattern="Search_SC[0-9]_Area_MC[0-9]")]
    
    data.morph.image[,var.new] = data.morph.image[,var.spot]/data.morph.image[,var.neurite[1]]
  }
}

var.group=c("Rep","Plate","Well","Condition","Treatment","Concentration","ID")
var.features=names(data.morph.image)[!names(data.morph.image) %in% c(var.group,"Image","File","V1")]

#Exclude nuclei measurements in other channels as well as summed features
if(any(grepl(pattern = "Nuclei",var.features))){
  var.nucl=var.features[grepl(pattern = "Nuclei",var.features)]
  nucl_channel=substring(unique(str_extract(var.nucl, pattern = "Nuclei_SC[0-9]")),10,11)
  expr=paste("Nuclei_SC",nucl_channel,".+MC1",sep="")
  var.excl=c()
  for(i in 1:4){
    if(i!=as.numeric(nucl_channel)){
      expr=paste("Nuclei_SC",nucl_channel,".+MC",i,sep="")
      var.excl=c(var.excl,var.features[grepl(x = var.features, pattern=expr)])
    }
  }
  expr=paste("Nuclei_SC",nucl_channel,".+Sum",sep="")
  var.excl=c(var.excl,var.features[grepl(x = var.features, pattern=expr)],
             var.features[grepl(x = var.features, pattern="_X_")],
             var.features[grepl(x = var.features, pattern="_Y_")],
             var.features[grepl(x = var.features, pattern="_XM_")],
             var.features[grepl(x = var.features, pattern="_YM_")],
             var.features[grepl(x = var.features, pattern="FeretAngle")],
             var.features[grepl(x = var.features, pattern="FeretX")],
             var.features[grepl(x = var.features, pattern="FeretY")],
             var.features[grepl(x = var.features, pattern=paste("Neurites_.+MC",nucl_channel,sep=""))],
             var.features[grepl(x = var.features, pattern=paste("Search_.+MC",nucl_channel,sep=""))])
  var.features=var.features[!(var.features %in% var.excl)]
}

#Exclude spot measurements in other channels as well as summed features
if(any(grepl(pattern = "Spot",var.features))){
  var.spot=var.features[grepl(pattern = "Spot_SC",var.features)]
  spot_channel=substring(unique(str_extract(var.spot, pattern = "Spot_SC[0-9]")),8,9)
  var.excl=c()
  for(c in spot_channel){
    for(i in 1:4){
      if(i!=as.numeric(c )){
        expr=paste("Spot_SC",c ,".+MC",i,sep="")
        var.excl=c(var.excl,var.features[grepl(pattern = expr,var.features)])
      }
    }
    expr=paste("Spot_SC",c,".+Sum",sep="")
    var.excl=c(var.excl,
               var.features[grepl(x = var.features, pattern=expr)],
               var.features[grepl(x = var.features, pattern=paste("Search_.+Area_MC",c,sep=""))],
               var.features[grepl(x = var.features, pattern=paste("Neurites_.+Area_MC",c,sep=""))])
  }
  var.features=var.features[!(var.features %in% var.excl)]
}

#Exclude search area measurements
if(any(grepl(pattern = "Search",var.features))){
  var.excl=var.features[grepl(pattern = "Search",var.features)]
  var.features=var.features[!(var.features %in% var.excl)]
}

#--4. Plot Raw Image Data-----------------------------------------------------
data.temp=data.morph.image
col=brewer.pal(length(unique(data.temp$ID)),"Set3")
i=1
plot_list = list()
var=var.features[1]
for(var in var.features){
  data.plot=data.temp[,c(var.group,var)]
  names(data.plot)[ncol(data.plot)]="variable"
  p=ggplot(data.plot,aes(x = Plate,y = variable, fill = ID)) +
    geom_violin(alpha=0.8)+
    scale_fill_manual(values=col, guide=guide_legend(nrow=2)) +
    facet_grid(Rep~.)+
    ylab(var)+
    theme_minimal(base_size = 10)+
    theme(legend.position = "bottom")
  p
  plot_list[[i]] = p
  i=i+1
}
ggsave(file=file.path(dir.plot,"Morph_Image_Violin.pdf"),marrangeGrob(plot_list, nrow=1, ncol=1),width = 30, height=21, units = "cm")

#--5. Well Data-----------------------------------------------------------------
#----5.1. Calculate Well Averages-----------------------------------------------
data.morph.well = data.morph.image %>%
  group_by_at(var.group)%>%
  summarise_at(var.features, mean, na.rm = TRUE)

#----5.2. Plate Overview  --------------------------------------------------------
data.temp=data.morph.well
data.temp$Column=as.numeric(substring(data.temp$Well, 2, nchar(as.character(data.temp$Well))))
data.temp$Column=factor(data.temp$Column, levels =(c(1,2,3,4,5,6,7,8,9,10,11,12)))
data.temp$Row=substring(data.temp$Well, 1,1)
data.temp$Row=factor(data.temp$Row, levels =(c("H","G","F","E","D","C","B","A")))

i=1
plot_list = list()
var=var.features
for(var in var.features){
  data.plot=data.temp[,c(var.group,"Column","Row",var)]
  names(data.plot)[ncol(data.plot)]="variable"
  p=ggplot(data.plot,aes(x = Column,y = Row, fill=variable)) +
    facet_grid(Rep~Plate)+
    geom_tile()+
    scale_x_discrete(drop=FALSE)+
    scale_y_discrete(drop=FALSE)+
    scale_fill_viridis_c()+
    coord_equal()+
    ggtitle(var)+
    xlab("Col")+
    ylab("Row")+
    theme_minimal(base_size = 10)
  p
  plot_list[[i]] = p
  i=i+1
}
ggsave(file=file.path(dir.plot,"Morph_Well_PlateOverview.pdf"),marrangeGrob(plot_list, nrow=1, ncol=1),width = 30, height=21, units = "cm")

data.plot=data.temp
col=brewer.pal(length(unique(data.plot$ID)),"Set3")
ggplot(data.plot,aes(x = Column,y = as.factor(Row), fill=ID)) +
  geom_tile()+
  scale_fill_manual(values=col)+
  scale_x_discrete(drop=FALSE)+
  scale_y_discrete(drop=FALSE)+
  coord_equal()+
  facet_grid(Rep~Plate)+
  xlab("Col")+
  ylab("Row")+
  theme_minimal(base_size = 10)+
  theme(legend.position = "bottom")
ggsave(file=file.path(dir.plot,"Morph_Well_PlateOverview_Layout.pdf"),width = 30, height=21, units = "cm")

#----5.3. Bar Plots-------------------------------------------------------------------
data.temp=data.morph.well
col=brewer.pal(length(unique(data.temp$ID)),"Set3")
i=1
plot_list = list()
for(var in var.features){
  data.plot=data.temp[,c(var.group,var)]
  names(data.plot)[ncol(data.plot)]="variable"
  p=ggplot(data.plot,aes(x = Plate,y = variable, fill=ID, group=ID)) +
    stat_summary(fun.data = mean_sdl, geom="errorbar", fun.args = list(mult=1), position=position_dodge(width = 0.9, preserve = "single"), width=.5) +
    stat_summary(fun = mean, geom = "bar", position=position_dodge(width = 0.9, preserve = "single"), width=0.85) +
    geom_point(position = position_dodge(width = .9), shape=21, fill="grey50", color="white")+
    scale_fill_manual(values=col, guide=guide_legend(nrow=2))+
    facet_grid(Rep~.)+
    ylab(var)+
    theme_minimal(base_size = 10)+
    theme(legend.position = "bottom")
  p
  plot_list[[i]] = p
  i=i+1
}
ggsave(file=file.path(dir.plot,"Morph_Well_BarPlot_Rep_Plate.pdf"),marrangeGrob(plot_list, nrow=1, ncol=1),width = 30, height=21, units = "cm")

i=1
plot_list = list()
for(var in var.features){
  data.plot=data.temp[,c(var.group,var)]
  names(data.plot)[ncol(data.plot)]="variable"
  p=ggplot(data.plot,aes(x = NA,y = variable, fill=ID, group=ID)) +
    stat_summary(fun.data = mean_sdl, geom="errorbar", fun.args = list(mult=1), position=position_dodge(width = 0.9, preserve = "single"), width=.5) +
    stat_summary(fun = mean, geom = "bar", position=position_dodge(width = 0.9, preserve = "single"), width=0.85) +
    geom_point(position = position_dodge(width = .9), shape=21, fill="grey50", color="white")+
    scale_fill_manual(values=col, guide=guide_legend(nrow=2))+
    ylab(var)+
    theme_minimal(base_size = 10)+
    theme(legend.position = "bottom")
  p
  plot_list[[i]] = p
  i=i+1
}
ggsave(file=file.path(dir.plot,"Morph_Well_BarPlot.pdf"),marrangeGrob(plot_list, nrow=2, ncol=2),width = 30, height=21, units = "cm")

#--6. Normalise Data -----------------------------------------------------------
#----6.1. Calculate Normalised Well Averages------------------------------------
unique(data.morph.image$ID)
data.morph.image.norm = data.morph.image%>%
  group_by(Rep,Plate)%>%
  mutate_at(c(var.features),funs(ZScore = (. - mean(.[ID==ctrl.condition], na.rm=TRUE))/sd(.[ID==ctrl.condition],na.rm=TRUE)))%>%
  select(Rep,Plate,Well,Image,Condition,Treatment,Concentration,ID,contains("ZScore"))

var.features.zscore=paste(var.features,"_ZScore",sep="")
data.morph.well.norm=data.morph.image.norm%>%
  group_by_at(c(var.group))%>%
  summarise_at(var.features.zscore, mean, na.rm=TRUE)%>%
  ungroup()

#----6.2. Plate Overview  --------------------------------------------------------
data.temp=data.morph.well.norm
data.temp$Column=as.numeric(substring(data.temp$Well, 2, nchar(as.character(data.temp$Well))))
data.temp$Column=factor(data.temp$Column, levels =(c(1,2,3,4,5,6,7,8,9,10,11,12)))
data.temp$Row=substring(data.temp$Well, 1,1)
data.temp$Row=factor(data.temp$Row, levels =(c("H","G","F","E","D","C","B","A")))

i=1
plot_list = list()
for(var in var.features.zscore){
  data.plot=data.temp[,c(var.group,"Column","Row",var)]
  names(data.plot)[ncol(data.plot)]="variable"
  p=ggplot(data.plot,aes(x = Column,y = Row, fill=variable)) +
    facet_grid(Rep~Plate)+
    geom_tile()+
    scale_x_discrete(drop=FALSE)+
    scale_y_discrete(drop=FALSE)+
    scale_fill_viridis_c()+
    coord_equal()+
    ggtitle(var)+
    xlab("Col")+
    ylab("Row")+
    theme_minimal(base_size = 10)
  p
  plot_list[[i]] = p
  i=i+1
}
ggsave(file=file.path(dir.plot,"Morph_Well_Norm_PlateOverview.pdf"),marrangeGrob(plot_list, nrow=1, ncol=1),width = 30, height=21, units = "cm")

#----6.3. Bar Plots-------------------------------------------------------------------
data.temp=data.morph.well.norm
col=brewer.pal(length(unique(data.temp$ID)),"Set3")
i=1
plot_list = list()
for(var in var.features.zscore){
  data.plot=data.temp[,c(var.group,var)]
  names(data.plot)[ncol(data.plot)]="variable"
  p=ggplot(data.plot,aes(x = Plate,y = variable, fill=ID, group=ID)) +
    stat_summary(fun.data = mean_sdl, geom="errorbar", fun.args = list(mult=1), position=position_dodge(width = 0.9, preserve = "single"), width=.5) +
    stat_summary(fun = mean, geom = "bar", position=position_dodge(width = 0.9, preserve = "single"), width=0.85) +
    geom_point(position = position_dodge(width = .9), shape=21, fill="grey50", color="white")+
    scale_fill_manual(values=col, guide=guide_legend(nrow=2))+
    facet_grid(Rep~.)+
    ylab(var)+
    theme_minimal(base_size = 10)+
    theme(legend.position = "bottom")
  p
  plot_list[[i]] = p
  i=i+1
}
ggsave(file=file.path(dir.plot,"Morph_Well_Norm_BarPlot_Rep_Plate.pdf"),marrangeGrob(plot_list, nrow=1, ncol=1),width = 30, height=21, units = "cm")

i=1
plot_list = list()
for(var in var.features.zscore){
  data.plot=data.temp[,c(var.group,var)]
  names(data.plot)[ncol(data.plot)]="variable"
  p=ggplot(data.plot,aes(x = NA,y = variable, fill=ID, group=ID)) +
    stat_summary(fun.data = mean_sdl, geom="errorbar", fun.args = list(mult=1), position=position_dodge(width = 0.9, preserve = "single"), width=.5) +
    stat_summary(fun = mean, geom = "bar", position=position_dodge(width = 0.9, preserve = "single"), width=0.85) +
    geom_point(position = position_dodge(width = .9), shape=21, fill="grey50", color="white")+
    scale_fill_manual(values=col, guide=guide_legend(nrow=2))+
    ylab(var)+
    theme_minimal(base_size = 10)+
    theme(legend.position = "bottom")
  p
  plot_list[[i]] = p
  i=i+1
}
ggsave(file=file.path(dir.plot,"Morph_Well_Norm_BarPlot.pdf"),marrangeGrob(plot_list, nrow=2, ncol=2),width = 30, height=21, units = "cm")

#--7. Export Data Frames -------------------------------------------------------
write.table(x=data.morph.image[,c(var.group,"Image",var.features)],file=file.path(dir.mergedData,"Data_Morph_Image.txt"), sep='\t', dec=".", col.names = TRUE, row.names = FALSE)
write.table(x=data.morph.image.norm[,c(var.group,"Image",var.features.zscore)],file=file.path(dir.mergedData,"Data_Morph_Image_Norm.txt"), sep='\t', dec=".",col.names = TRUE, row.names = FALSE)
write.table(x=data.morph.well[,c(var.group,var.features)],file=file.path(dir.mergedData,"Data_Morph_Well.txt"), sep='\t', dec=".",col.names = TRUE, row.names = FALSE)
write.table(x=data.morph.well.norm[,c(var.group,var.features.zscore)],file=file.path(dir.mergedData,"Data_Morph_Well_Norm.txt"), sep='\t', dec=".",col.names = TRUE, row.names = FALSE)

