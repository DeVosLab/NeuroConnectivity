# ------------------------------------------------------------------------------
# NeuroConnectivity - Functional calcium analysis
#
# Author: Winnok H. De Vos
# Modified by: Marlies Verschuuren
# Creation date: 2019-12-13
# Last Modified: 2023-12-20
# ------------------------------------------------------------------------------

#--1. User settings-------------------------------------------------------------
#----1.1. Select directories----------------------------------------------------
# Input: Folder with structure: Rep > Plate > Func > Output
#                                      Rep > Plate > PlateLayout.txt
dir.input="/Users/marliesverschuuren/Documents/UA_DataSets/NeuroConnectivity/PLA/Data"
dir.output="/Users/marliesverschuuren/Library/CloudStorage/OneDrive-UniversiteitAntwerpen/Projects/DeVosLab/NeuroConnectivity/Results_PLA"

#----1.2. Settings analysis------------------------------------------------------
interval = 0.5 # (500 ms)
ctrl.condition="B27_NA_NA" #Condition_Treatment_Concentration
peakheight=1.05 #Height normalised peak (1.05 = 5% increase from median intensity)
peakdistance=5 #Number of frames peak
activePeakNr=5 #Number of peaks to be considered active 

#--2. Packages and Settings-----------------------------------------------------
#----2.1. Packages--------------------------------------------------------------
if (!require("tidyverse")) {install.packages("tidyverse"); require("tidyverse")}
if (!require("data.table")) {install.packages("data.table"); require("data.table")}
if (!require("RColorBrewer")) {install.packages("RColorBrewer"); require("RColorBrewer")}   #Function: brewer.pal
if (!require("plyr")) {install.packages("plyr"); require("plyr")}                           #Function: ddply >> Not compatible with dplyr >> Specify dplyr functions
if (!require("gridExtra")) {install.packages("gridExtra"); require("gridExtra")}            #Function: marrangegrob
if (!require("pracma")) {install.packages("pracma"); require("pracma")}                     #Find peaks
if (!require("doParallel")) {install.packages("doParallel"); require("doParallel")}         #Parallel computing on multiple cores

#----2.2. Create result directories---------------------------------------------
dir.plot=file.path(dir.output, "Neuro_Func_Plots")
dir.create(dir.plot, showWarnings = FALSE)

dir.mergedData=file.path(dir.output, "Neuro_Func_MergedData")
dir.create(dir.mergedData, showWarnings = FALSE)

#----2.3. Functions-------------------------------------------------------------
# bleach correction (requires index x, intensity y)
bleachcorrect = function(x, y){
  params = c()
  it = 0
  model = NULL
  while(TRUE & it<=1)
  { 
    # try exp fit
    try(model<-nls(y ~ SSasymp(x, yf, y0, log_alpha)), silent = T)
    it = it + 1
  }
  if(is.null(model)) 
  {
    # if exp fit doesn't work, lin fit 
    model = lm(y ~ x)
    p = summary(model)$coefficients[2,4]
  }else p=0
  if(p<0.05){ # fit needs to have a significant slope coeff
    y_bleach = predict(model)
    y_corr = y/y_bleach*mean(y,na.rm=T)  
  }else y_corr=y # no correction applied
}

# normalization functions 
norm.fun = function(x) x/median(x, na.rm=T) # divide by the median of the time trace

#--3. Read Func Data ----------------------------------------------------------------
data.roi = data.frame()
folders.rep=list.files(path=file.path(dir.input))
data.list=list()
i=1
#Loop over all rep folders 
for (rep in folders.rep){ 
  folders.plate=list.files(path=file.path(dir.input,rep))
  #Loop over all Plate folders in Rep folder
  for (plate in folders.plate){
    data.plate = data.frame()
    files = list.files(path = file.path(dir.input,rep,plate,"Func","Output"), pattern = "results", full.names = TRUE)
    
    expr="_[A-Z][0-9]{2}_" #### ADAPT IF NEEDED
    #Check regexpr:
    #File=files[1]
    #str_extract(File,expr)
    
    #Read plate
    data.plate = tibble(File = files) %>%
      dplyr::mutate(data = lapply(File, fread)) %>%
      unnest(data)%>%
      dplyr::mutate(Rep=rep,
               Plate=plate,
               File=File,
               Image=substring(str_extract(File,"Output/.+"),8,nchar(str_extract(File,"Output/.+"))-12),
               Well=substring(str_extract(File,expr),2,4)) %>% #### ADAPT IF NEEDED
      gather(ROI_Image,Mean,-Rep,-Plate,-File,-Image,-Well,-V1)
    
    #Omit missing data
    data.plate = na.omit(data.plate)
    
    #Change variable names
    data.plate$ROI_Image=as.factor(substring( data.plate$ROI_Image,5,nchar(data.plate$ROI_Image)))
    names(data.plate)[which(names(data.plate)=="V1")]="Frame" 
    
    #Add image index and roi index within well
    data.plate=data.plate%>%
      dplyr::group_by(Rep,Plate,Well)%>%
      dplyr::mutate(ImageNr = cumsum(!duplicated(Image)),
                    ROI_Well = cumsum(!duplicated(ROI_Image)))
    
    #Add 0 if well is defined in format B2
    data.plate$Well=ifelse(nchar(data.plate$Well)==2,paste(substring(data.plate$Well,1,1),"0",substring(data.plate$Well,2,2),sep=""),data.plate$Well)
    
    #Merge with plate layout
    layout=read.table(file.path(dir.input,rep,plate,"PlateLayout.txt"),header=TRUE, sep='\t', fill=TRUE)
    data.plate=data.plate%>%
      dplyr::left_join(layout)
    
    data.list[[i]]=data.plate
    i=i+1
  }
}

#Bind all replicates and plates, and assign roi index for whole data set
data.roi=rbindlist(data.list)%>%
  as.data.frame()%>%
  dplyr::group_by(Rep,Plate,Well,ImageNr,ROI_Image)%>%
  dplyr::mutate(ROI_Data = cur_group_id())%>%
  dplyr::ungroup()

#Specify time in seconds
data.roi$Time = data.roi$Frame*interval

#Change type variables
data.roi$Rep = as.factor(data.roi$Rep)
data.roi$Plate = as.factor(data.roi$Plate)
data.roi$Well = as.factor(data.roi$Well)
data.roi$ROI_Image = as.factor(data.roi$ROI_Image)
data.roi$ROI_Image=factor(data.roi$ROI_Image, levels = sort(unique(as.numeric(data.roi$ROI_Image))))
data.roi$ROI_Well = as.factor(data.roi$ROI_Well)
data.roi$ROI_Well=factor(data.roi$ROI_Well, levels = sort(unique(as.numeric(data.roi$ROI_Well))))
data.roi$ROI_Data = as.factor(data.roi$ROI_Data)
data.roi$ROI_Data=factor(data.roi$ROI_Data, levels = sort(unique(as.numeric(data.roi$ROI_Data))))

#--4. Plot Raw Data -------------------------------------------------------------
data.temp=data.roi
i=1
plot_list = list()
for(rep in unique(data.temp$Rep)){
  data.plot.rep=data.temp[data.temp$Rep==rep,]
  for(plate in unique(data.plot.rep$Plate)){
    data.plot.plate=data.plot.rep[data.plot.rep$Plate==plate,]
    for(well in unique(data.plot.plate$Well)){
      data.plot.well=data.plot.plate[data.plot.plate$Well==well,]
      for(image in unique(data.plot.well$Image)){
        data.plot=data.plot.well[data.plot.well$Image==image,]
        p=ggplot(data.plot,aes(Time, Mean, group = as.numeric(ROI_Image), colour = as.numeric(ROI_Image))) +
          geom_path() +
          scale_colour_viridis_c()+
          facet_wrap(ROI_Image~.) +
          theme_minimal()+
          ggtitle(paste(rep,plate,well,image,sep="_"))+
          theme(legend.position = "bottom")
        p
        plot_list[[i]] = p
        i=i+1
      }
    }
  }
}
ggsave(file=file.path(dir.plot,"Func_ROI_Raw.pdf"),marrangeGrob(plot_list, nrow=1, ncol=1),width = 30, height=21, units = "cm")

data.temp=data.roi
i=1
plot_list = list()
for(rep in unique(data.temp$Rep)){
  data.plot.rep=data.temp[data.temp$Rep==rep,]
  for(plate in unique(data.plot.rep$Plate)){
    data.plot.plate=data.plot.rep[data.plot.rep$Plate==plate,]
    for(well in unique(data.plot.plate$Well)){
      data.plot=data.plot.plate[data.plot.plate$Well==well,]
      p=ggplot(data.plot,aes(Time, Mean, group = as.numeric(ROI_Image), colour = as.numeric(ROI_Image))) +
        geom_path() +
        facet_wrap(Image~.) +
        scale_colour_viridis_c()+
        theme_minimal()+
        ggtitle(paste(rep,plate,well,sep="_"))+
        theme(legend.position = "bottom")
      p
      plot_list[[i]] = p
      i=i+1
    }
  }
}
ggsave(file=file.path(dir.plot,"Func_Image_Raw.pdf"),marrangeGrob(plot_list, nrow=1, ncol=1),width = 30, height=21, units = "cm")

#--5. Bleach correction, normalisation and smoothing per ROI--------------------
data.roi = ddply(data.roi, "ROI_Data", transform, corr.mean = bleachcorrect(Frame, Mean))
data.roi = ddply(data.roi, "ROI_Data", transform, norm.mean = norm.fun(corr.mean))
data.roi = ddply(data.roi, "ROI_Data", transform, smooth.mean = loess(norm.mean ~ Frame, span = 0.02)$fitted) # Very fine smoothing to blunt noise (enhances peak detection) - currently uses 2% points for smoothing

#--6. Plot Normalised data -----------------------------------------------------
data.temp=data.roi
i=1
plot_list = list()
col=c("MeanInt"="black","CorrInt"="dodgerblue","NormInt"="firebrick","SmoothInt"="goldenrod")
for(rep in unique(data.temp$Rep)){
  data.plot.rep=data.temp[data.temp$Rep==rep,]
  for(plate in unique(data.plot.rep$Plate)){
    data.plot.plate=data.plot.rep[data.plot.rep$Plate==plate,]
    for(well in unique(data.plot.plate$Well)){
      data.plot.well=data.plot.plate[data.plot.plate$Well==well,]
      for(image in unique(data.plot.well$Image)){
        data.plot=data.plot.well[data.plot.well$Image==image,]
        factorAxis=mean(data.plot$Mean)
        p=ggplot(data.plot,aes(group = as.numeric(ROI_Image))) +
          geom_path(aes(Time, Mean, color = "MeanInt"), linewidth=0.2) +
          geom_path(aes(Time, corr.mean, color="CorrInt"),linewidth=0.2) +
          geom_path(aes(Time, norm.mean * factorAxis, color="NormInt"),linewidth=0.2) +
          geom_path(aes(Time, smooth.mean * factorAxis, color="SmoothInt"),linewidth=0.2) +
          scale_color_manual(values = col)+
          scale_y_continuous(name = "MeanInt or CorrInt", sec.axis = sec_axis(~./factorAxis, name="NormInt or SmoothInt"))+
          facet_wrap(ROI_Image~.) +
          labs(color="Legend")+
          theme_minimal()+
          ggtitle(paste(rep,plate,well,image,sep="_"))+
          theme(legend.position = "bottom")
        p
        plot_list[[i]] = p
        i=i+1
      }
    }
  }
}
ggsave(file=file.path(dir.plot,"Func_ROI_Norm.pdf"),marrangeGrob(plot_list, nrow=1, ncol=1),width = 30, height=21, units = "cm")

#--7. Peak detection -----------------------------------------------------------
#----7.1. Peak Detection ROI----------------------------------------------------
data.roi.peak = data.frame()
total.roi.nr = length(unique(data.roi$ROI_Data))

#Define cores for parallel computing
cores=detectCores()
cl = makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

#Peak detection
data.roi.peak = foreach(i=1:total.roi.nr, .combine=rbind) %dopar% {
  temp = data.roi[which(data.roi$ROI_Data == levels(data.roi$ROI_Data)[i]),]
  temp$duration = NA
  temp$interval = NA
  #findpeaks returns a matrix where each row represents one peak found. The first column gives the height, the second the position/index where the maximum is reached, the third and forth the indices of where the peak begins and ends
  peaks = as.data.frame(pracma::findpeaks(temp$smooth.mean, nups = 2, ndowns = 2, minpeakheight = peakheight, minpeakdistance=peakdistance))
  if(dim(peaks)[1]!=0){
    peaks=peaks[order(peaks$V2),]
  }
  peaks$duration = (peaks$V4 - peaks$V3) * interval # exact duration of detected wave in sec
  temp$peak = ifelse(temp$Frame%in%(peaks$V2),temp$smooth.mean,NA) # labels peak with smooth mean value
  temp$peak.start = ifelse(temp$Frame%in%(peaks$V3),temp$smooth.mean,NA) # labels peak with smooth mean value  
  temp$peak.stop = ifelse(temp$Frame%in%(peaks$V4),temp$smooth.mean,NA) # labels peak with smooth mean value  
  for (f in 1: dim(peaks)[1])
  {
    temp$duration[peaks$V4[f]] = peaks$duration[f]
    temp$interval[peaks$V2[f]] = ifelse(f==1, peaks$V2[f]* interval, (peaks$V2[f]-peaks$V2[f-1]) * interval )
  }
  temp #Equivalent to data.roi.peak = rbind(data.roi.peak, temp)
}

#stop cluster
stopCluster(cl)

#----7.2. Peak Detection Average Trace Image------------------------------------
data.image=data.roi%>%
  dplyr::group_by(File,Frame,Rep,Plate,Well,Image,ImageNr,Condition,Treatment, Concentration,Time)%>%
  dplyr::summarise(norm.mean=mean(norm.mean,na.rm = TRUE),
                   smooth.mean=mean(smooth.mean,na.rm = TRUE),
                   Mean=mean(Mean,na.rm = TRUE))%>%
  dplyr::ungroup()%>%
  dplyr::mutate(Image_ID=paste(Rep,Plate,Well,ImageNr,sep="_")
)
data.image$Image_ID=as.factor(data.image$Image_ID)

data.image.peak = data.frame()
total.roi.nr = length(unique(data.image$Image_ID))

#Define cores for parallel computing
cores=detectCores()
cl = makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

#Peak detection
data.image.peak = foreach(i=1:total.roi.nr, .combine=rbind) %dopar% {
  temp = data.image[which(data.image$Image_ID == levels(data.image$Image_ID)[i]),]
  temp$duration = NA
  temp$interval = NA
  #findpeaks returns a matrix where each row represents one peak found. The first column gives the height, the second the position/index where the maximum is reached, the third and forth the indices of where the peak begins and ends
  peaks = as.data.frame(pracma::findpeaks(temp$smooth.mean, nups = 2, ndowns = 2, minpeakheight = peakheight, minpeakdistance=peakdistance))
  if(dim(peaks)[1]!=0){
    peaks=peaks[order(peaks$V2),]
  }
  peaks$duration = (peaks$V4 - peaks$V3) * interval # exact duration of detected wave in sec
  temp$peak = ifelse(temp$Frame%in%(peaks$V2),temp$smooth.mean,NA) # labels peak with smooth mean value
  temp$peak.start = ifelse(temp$Frame%in%(peaks$V3),temp$smooth.mean,NA) # labels peak with smooth mean value  
  temp$peak.stop = ifelse(temp$Frame%in%(peaks$V4),temp$smooth.mean,NA) # labels peak with smooth mean value  
  for (f in 1: dim(peaks)[1])
  {
    temp$duration[peaks$V4[f]] = peaks$duration[f]
    temp$interval[peaks$V2[f]] = ifelse(f==1, peaks$V2[f]* interval, (peaks$V2[f]-peaks$V2[f-1]) * interval )
  }
  temp #Equivalent to data.roi.peak = rbind(data.roi.peak, temp)
}

#stop cluster
stopCluster(cl)

#----7.3. Plot peaks -----------------------------------------------------------
#Plot peaks per trace
data.temp=data.roi.peak
col=colorRampPalette(c("black","violetred3"))(20)
i=1
for(rep in unique(data.temp$Rep)){
  data.plot.rep=data.temp[data.temp$Rep==rep,]
  for(plate in unique(data.plot.rep$Plate)){
    data.plot.plate=data.plot.rep[data.plot.rep$Plate==plate,]
    for(well in unique(data.plot.plate$Well)){
      data.plot.well=data.plot.plate[data.plot.plate$Well==well,]
      for(image in unique(data.plot.well$Image)){
        data.plot=data.plot.well[data.plot.well$Image==image,]
        p=ggplot(data.plot) +
          geom_path(aes(Time,smooth.mean),color = "black", linewidth=0.5) +
          geom_point(aes(Time,peak), color="red", size=0.5)+
          geom_point(aes(Time,peak.start), color="orange", size=0.5)+
          facet_wrap(~ROI_Image) +
          ggtitle(paste(rep,plate,well,image,sep="_"))+
          theme_minimal() +
          theme(legend.position = "none")
        p
        plot_list[[i]] = p
        i=i+1
      }
    }
  }
}
ggsave(file=file.path(dir.plot,"Func_ROI_Peak.pdf"),marrangeGrob(plot_list, nrow=1, ncol=1),width = 30, height=21, units = "cm")

#Plot peaks on average trace per image
data.temp=data.image.peak
col=colorRampPalette(c("black","violetred3"))(20)
i=1
for(rep in unique(data.temp$Rep)){
  data.plot.rep=data.temp[data.temp$Rep==rep,]
  for(plate in unique(data.plot.rep$Plate)){
    data.plot.plate=data.plot.rep[data.plot.rep$Plate==plate,]
    for(well in unique(data.plot.plate$Well)){
      data.plot=data.plot.plate[data.plot.plate$Well==well,]
      p=ggplot(data.plot) +
        geom_path(aes(Time,smooth.mean),color = "black", linewidth=0.5) +
        geom_point(aes(Time,peak), color="red", size=0.5)+
        geom_point(aes(Time,peak.start), color="orange", size=0.5)+
        facet_wrap(.~Image,scales = "free_y") +
        ggtitle(paste(rep,plate,well,sep="_"))+
        theme_minimal() +
        theme(legend.position = "none")
      p
      plot_list[[i]] = p
      i=i+1
    }
  }
}

ggsave(file=file.path(dir.plot,"Func_Image_AvTrace_Peak.pdf"),marrangeGrob(plot_list, nrow=1, ncol=1),width = 30, height=21, units = "cm")

# Plot peaks per well in a tile plot 
data.temp=data.roi.peak
i=1
plot_list = list()
for(rep in unique(data.temp$Rep)){
  data.plot.rep=data.temp[data.temp$Rep==rep,]
  for(plate in unique(data.plot.rep$Plate)){
    data.plot=data.plot.rep[data.plot.rep$Plate==plate,]
    
    p=ggplot(data.plot, aes(Time, as.factor(ROI_Well))) +
      geom_tile(data = subset(data.plot, peak>0), fill="black") +
      facet_wrap(Well~Image,scales = "free_y") +
      ggtitle(paste(rep,plate,sep="_"))+
      theme_minimal(base_size = 10) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.y=element_blank(),
            legend.position="none") 
    p
    plot_list[[i]] = p
    i=i+1
  }
}
ggsave(file=file.path(dir.plot,"Func_ROI_Peak_WellTile.pdf"),marrangeGrob(plot_list, nrow=1, ncol=1),width = 30, height=21, units = "cm")

#--8. Activity stats -----------------------------------------------------------
#----8.1. Image Data: Synchronicity between ROIs within an image  --------------
cor.data = data.frame()
data.roi.peak$Image_ID=paste(data.roi.peak$Rep,data.roi.peak$Plate,data.roi.peak$Well,data.roi.peak$ImageNr,sep="_")
data.roi.peak$Image_ID=as.factor(data.roi.peak$Image_ID)
image.nr = length(unique(data.roi.peak$Image_ID))
for (i in 1:image.nr){
  temp = data.roi.peak[which(data.roi.peak$Image_ID == levels(data.roi.peak$Image_ID)[i]),]
  df = as.data.frame(cbind(temp$norm.mean,temp$Time,temp$ROI_Data))
  names(df) = c("norm.mean","Time","ROI_Data")
  df.d = reshape2::dcast(data = df,formula = Time~ROI_Data, value.var = "norm.mean") # convert dataframe to matrix of ROIs (only retain the norm values per ROI as a function of time)
  df.matrix = as.matrix(df.d[, -1]) # remove time
  correlations = cor(df.matrix, use="pairwise.complete.obs") 
  correlations = correlations[upper.tri(correlations, diag=F)] # only retain unique correlations
  synchronous.fraction = sum(correlations>0.7)/length(correlations)*100 # consider correlated neurons the ones with more than 70% correlation
  av.correlation = mean(correlations,na.rm=T)
  cor.data = rbind(cor.data,c(levels(data.roi.peak$Image_ID)[i],synchronous.fraction,av.correlation))
}
names(cor.data) = c("Image_ID","synchronous.fraction","av.correlation")
cor.data$synchronous.fraction=as.numeric(cor.data$synchronous.fraction)
cor.data$av.correlation=as.numeric(cor.data$av.correlation)

#----8.2. Image Data: Global parameters per ROI  -------------------------------
features.roi = ddply(data.roi.peak, c("Rep","Plate","Well","Image","ImageNr","Image_ID","Condition","Treatment","Concentration","ROI_Data"), summarise,
  baseline.intensity = median(Mean),
  peak.nr = length(peak[!is.na(peak)]),
  peak.frequency = peak.nr/max(Time),
  peak.duration = mean(duration, na.rm=T),
  peak.interval = mean(interval, na.rm=T),
  peak.duration.variability = sd(duration, na.rm=T)/mean(duration, na.rm=T),
  peak.interval.variability = sd(interval, na.rm=T)/mean(interval, na.rm=T),
  dynamic.range = max(Mean)- min(Mean), 
  dynamic.range.norm = max(norm.mean) - min(norm.mean))

features.image.a = ddply(features.roi, c("Rep","Plate","Well","Image","ImageNr","Image_ID","Condition","Treatment","Concentration"), summarise,
                     baseline.intensity = mean(baseline.intensity),
                     peak.frequency = mean(peak.frequency, na.rm=T),
                     peak.duration = mean(peak.duration, na.rm=T),
                     peak.interval = mean(peak.interval, na.rm=T),
                     peak.duration.variability = mean(peak.duration.variability, na.rm=T),
                     peak.interval.variability = mean(peak.interval.variability, na.rm=T),
                     dynamic.range = mean(dynamic.range),
                     dynamic.range.norm = mean(dynamic.range.norm))

# summarize roi data per image - with a cutoff for inactive ROIs showing less than 10 peaks across the time window
features.image.b = ddply(features.roi, c("Rep","Plate","Well","Image","ImageNr","Image_ID","Condition","Treatment","Concentration"), summarise,
        active.fraction = sum(peak.nr>activePeakNr)/length(peak.nr)*100)

features.image.c = ddply(features.roi[features.roi$peak.nr>activePeakNr,], c("Rep","Plate","Well","Image","ImageNr","Image_ID","Condition","Treatment","Concentration"), summarise,
                     act.baseline.intensity = mean(baseline.intensity),
                     act.peak.frequency = mean(peak.frequency),
                     act.peak.duration = mean(peak.duration),
                     act.peak.interval = mean(peak.interval),
                     act.peak.duration.variability = mean(peak.duration.variability),
                     act.peak.interval.variability = mean(peak.interval.variability),
                     act.dynamic.range = mean(dynamic.range),
                     act.dynamic.range.norm = mean(dynamic.range.norm))

features.image.bc=left_join(features.image.b,features.image.c)

#----8.3. Image Data: Global parameters of average trace  -------------------------------
features.image.d=ddply(data.image.peak, c("Rep","Plate","Well","Image","ImageNr","Image_ID","Condition","Treatment","Concentration"), summarise,
  avTrace.baseline.intensity = median(Mean),
  avTrace.peak.nr = length(peak[!is.na(peak)]),
  avTrace.peak.frequency = avTrace.peak.nr/max(Time), # peaks per sec
  avTrace.dynamic.range = max(Mean) - min(Mean), 
  avTrace.dynamic.range.norm = max(norm.mean) - min(norm.mean))

features.image.d=subset(features.image.d, select=-c(avTrace.peak.nr))

# combine all image data
data.func.image = left_join(features.image.a,features.image.bc)%>%
  dplyr::left_join(features.image.d)%>%
  dplyr::left_join(.,cor.data)

data.func.image$ID=paste(data.func.image$Condition,data.func.image$Treatment, data.func.image$Concentration, sep="_")

#----8.4. Well Data: -----------------------------------------------------------
var.features=names(data.func.image)[!names(data.func.image)%in% c("Rep","Plate","Well","Image","ImageNr","Image_ID","Condition","Treatment","Concentration","ID")]
data.func.well=data.func.image%>%
  dplyr::group_by(Rep,Plate,Well,Condition,Treatment,Concentration,ID)%>%
  dplyr::summarise(across(var.features, mean, na.rm = TRUE))

#----8.5. Normalise to control condition----------------------------------------
unique(data.func.image$ID)

data.func.image.norm = data.func.image%>%
  dplyr::group_by(Rep,Plate)%>%
  dplyr:: mutate_at(c(var.features),funs(ZScore = (. - mean(.[ID==ctrl.condition], na.rm=TRUE))/sd(.[ID==ctrl.condition],na.rm=TRUE)))%>%
  dplyr::select(Rep,Plate,Well,ImageNr,Image_ID,Condition,Treatment,Concentration,ID,contains("ZScore"))

var.features.zscore=paste(var.features,"_ZScore",sep="")
data.func.well.norm=data.func.image.norm%>%
  dplyr::group_by(Rep,Plate,Well,Condition,Treatment,Concentration,ID)%>%
  dplyr::summarise(across(var.features.zscore, mean, na.rm = TRUE))

#--9. Result plot --------------------------------------------------------------
#----9.1. Well Data ------------------------------------------------------------
data.temp=data.func.well
var.group=c("Rep","Plate","Well","Condition","Treatment","Concentration","ID")
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

ggsave(file=file.path(dir.plot,"Func_Well_Barplot_Rep_Plate.pdf"),marrangeGrob(plot_list, nrow=1, ncol=1),width = 30, height=21, units = "cm")

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

ggsave(file=file.path(dir.plot,"Func_Well_Barplot.pdf"),marrangeGrob(plot_list, nrow=2, ncol=2),width = 30, height=21, units = "cm")

#----9.2. Well Norm Data -------------------------------------------------------
data.temp=data.func.well.norm
var.group=c("Rep","Plate","Well","Condition","Treatment","Concentration","ID")
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

ggsave(file=file.path(dir.plot,"Func_Well_Norm_Barplot_Rep_Plate.pdf"),marrangeGrob(plot_list, nrow=1, ncol=1),width = 30, height=21, units = "cm")

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

ggsave(file=file.path(dir.plot,"Func_Well_Norm_Barplot.pdf"),marrangeGrob(plot_list, nrow=2, ncol=2),width = 30, height=21, units = "cm")

#--10. Export Data Frames -------------------------------------------------------
write.table(x=data.roi.peak,file=file.path(dir.mergedData,"Data_Func_ROI.txt"), sep='\t', dec=".", col.names = TRUE, row.names = FALSE)
write.table(x=data.func.image,file=file.path(dir.mergedData,"Data_Func_Image.txt"), sep='\t', dec=".", col.names = TRUE, row.names = FALSE)
write.table(x=data.func.well,file=file.path(dir.mergedData,"Data_Func_Well.txt"), sep='\t', dec=".",col.names = TRUE, row.names = FALSE)
write.table(x=data.func.image.norm,file=file.path(dir.mergedData,"Data_Func_Image_Norm.txt"), sep='\t', dec=".", col.names = TRUE, row.names = FALSE)
write.table(x=data.func.well.norm,file=file.path(dir.mergedData,"Data_Func_Well_Norm.txt"), sep='\t', dec=".",col.names = TRUE, row.names = FALSE)

