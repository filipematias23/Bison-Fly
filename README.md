# [Bison-Fly](https://github.com/filipematias23/Bison-Fly): UAV pipeline at NDSU Spring Wheat Breeding Program

The **Bison-Fly** tutorial is the *NDSU Spring Wheat UAV Pipeline* developed in partnership with *Drone2Phenome* ([D2P](https://www.ag2pi.org/)). In this pipeline we are presenting step by step how we have been applying UAV data on our breeding program. This is an open source R code from where anyone can use and adept for different crops. Hope this pipeline help you with your research and all suggestions to improve it are welcome. Help us to build this platform in partnership with all D2P community. 

<p align="center">
  <img src="https://raw.githubusercontent.com/filipematias23/images/master/readme/BisonFly.jpg" width="40%" height="40%">
</p>

<div id="menu" />

## Resources
  
   * [Introduction](#Intro)
   * [Image analysis in R](#R1)
   * [1. Selecting the targeted field](#P1)
   * [3. Rotating the image](#P2)

<br />

<div id="Intro" />

## Introduction 

> Welcome to the **NDSU Spring Wheat Breeding Program** lead by [Dr. Andrew Green](https://www.ndsu.edu/agriculture/ag-home/directory/andrew-green). In this tutorial, we will be presenting how we have been implementing UAV data on our breeding program. This pipeline is a compilation of [R]( https://www.r-project.org/) codes and functions from different packages to implement high-throughput phenotyping methods and understand spring wheat traits.

<p align="center">
  <img src="https://raw.githubusercontent.com/filipematias23/images/master/readme/BF_1a.jpg">
</p>

<br />

> During 2021, our team used the drone [Inspire 2 from DJI]( https://www.dji.com/inspire-2) combined with the multispectral sensor [Sentera 6X]( https://sentera.com/data-capture/6x-multispectral/). For instance, our phenotyping kit is composed by three set of batteries that gives us an autonomy around one hour of flying time. The kit also has one [reflectance panel]( https://support.sentera.com/portal/en/kb/articles/capturing-your-refelectance-panel) for radiometric calibration. This step is important to transform digital numbers to radiance and then reflectance. In addition, this step allows comparing data from different days or locations, which is very common on plant breeding programs. 

<p align="center">
  <img src="https://raw.githubusercontent.com/filipematias23/images/master/readme/BF_2.jpg" width="70%" height="70%">
</p>

<br />

> The flight mission was planed using the software [FieldAgent from Sentera]( https://sentera.com/fieldagent/). The idea is drawing the trial using the desktop version before going to field. Once in the flying area or trial was established, the flights parameters were set for 75% of overlapping, 15 mph, and flying height of 200 feet. There are other great open source options for DJI drones as [Pix4Dcapture]( https://support.pix4d.com/hc/en-us/articles/202557269-Pix4Dcapture-Getting-Started) and [DroneDeploy]( https://www.dronedeploy.com/product/mobile/). 

<br />

> For georeferencing, we have been using white plastic bucket lids (3.6 gallons) as geographic control points (GCPs). These lids are fixed in the ground with metal sticks in each corners of the field trial at the beginning of planting season. Observe that we normally paint the central part of the lid to facilitate finding where to click to georeferencing this point during the orthomosaic step. We also make some holes to drain the water from rain or irrigation. The GCPs stay in the field during the entire season and must be imagery alongside with the plants for each flight or data collection. This is important because allows a better overlayer of the experimental plots and facilitate drawing the grid plot polygons only once. Please for more information read this link: [https://github.com/OpenDroneMap/FIELDimageR#P5](https://github.com/OpenDroneMap/FIELDimageR#P5). 

<p align="center">
  <img src="https://raw.githubusercontent.com/filipematias23/images/master/readme/BF_3.jpg">
</p>

<br />

> Another simple idea to facilitate finding the corners on the field trial and draw the experimental plot polygons is to use four of the same plastic lids as extra points between borders and plots. This is an easy way to find where to click when using the function fieldShape() from [FIELDimageR]( https://github.com/OpenDroneMap/FIELDimageR#P5). Also, this method helps to identify where are the boundaries of multiple experimental trials plated in sequence as highlighted in the example below from our breeding program.

<p align="center">
  <img src="https://raw.githubusercontent.com/filipematias23/images/master/readme/BF_4.jpg">
</p>

<br />

> In this tutorial we evaluate a population of 106 white wheat lines on second year of yield trial (breeding pipeline). At the total 15 flights were collected in 2021 at NDSU Casselton Seed Farm in North Dakota (46.88051067579893, -97.2370382539691). All genotypes had at least two replicates, some check lines had three replicates (e.g., ALPINE, BARLOW, ELGIN-ND, FALLER, GLENN, SYINGMAR, NDSW0932, NDSW14098). The population was evaluated for maturity day (MAT_DAY), plant height (HT), days to heading (DH), lodging (LODG), and yield (YLD). All flights were performed around noon on days without clouds. The flights dates were (1) 05/11, (2) 05/23, (3) 05/26, (4) 06/01, (5) 06/07, (6) 06/10, (7) 06/17, (8) 06/22, (9) 06/25, (10) 06/28, (11) 07/01, (12) 07/13, (13) 07/19, (14) 07/20, and (15) 07/26.  

<p align="center">
  <img src="https://raw.githubusercontent.com/filipematias23/images/master/readme/BF_5.jpg">
</p>

<br />

> As an open source software to do the orthomosaic step we suggest the [OpenDroneMap]( https://www.opendronemap.org/) with the click and point interface [WebODM]( https://www.opendronemap.org/webodm/). A workshop tutorial video is available at [Phenome-Force YouTube Channel]( https://www.youtube.com/watch?v=U-bsA7QjzYE&t=4842s). More information is avaliable [HERE](https://github.com/OpenDroneMap/FIELDimageR/blob/master/README.md#P18).

<p align="center">
  <img src="https://raw.githubusercontent.com/filipematias23/images/master/readme/ODM_FIELDimageR_New.jpg" width="70%" height="70%">
</p>

<div id="R1" />

[Menu](#menu)

## Image analysis in R

> The first step is to [DOWNLOAD]() the data used in this tutorial by cliking [HERE]().

```r
##########################################
### Bison-Fly: Plant Breeding Pipeline ###
##########################################

### Necessary packages ###
library(FIELDimageR)
library(raster)
library(rgdal)
library(ggplot2)
library(DescTools)
library(lme4)
library(emmeans)
library(reshape2)
library(car)

#####################################
### NDSU - Spring Wheat UAV Data  ###
#####################################

### List of orthomosaics ###
Field <- list.files("./5band/") # 15 5band orthomosaics
Field_RGB <- list.files("./RGB/") # 15 RGB orthomosaics
Field_DSM <- list.files("./DSM/") # 15 DSM orthomosaics

######################################################
### Using one orthomosaic as base to draw the grid ###
######################################################

### Suggested strategy of using lids to highlight where to click ###
Test <- stack(paste("./5band/",Field[4],sep = ""))
plotRGB(FIELDimageR:::RGB.rescale(Test,3))

##########################
### Plot Polygons Grid ###
##########################

Data <- read.csv("EX_DATA.csv",header = T,fileEncoding="UTF-8-BOM")
Map <- read.csv("EX_MAP.csv",header = F,fileEncoding="UTF-8-BOM")
Shapefile <- fieldShape(mosaic = Test,
                        ncols = 11,
                        nrows = 20,
                        fieldData = Data,
                        ID = "PLOT",
                        fieldMap = Map)
Shapefile<-Shapefile$fieldShape
plotRGB(FIELDimageR:::RGB.rescale(Test,3))
plot(Shapefile, border="green",add=T)

### Saving Shapefile ###
library(rgdal)
writeOGR(Shapefile, ".", "Shapefile", driver="ESRI Shapefile")

```

<p align="center">
  <img src="https://raw.githubusercontent.com/filipematias23/images/master/readme/BF_6.jpg" width="60%" height="60%">
</p>

<br />

```r
############
### Mask ###
############

### Removing soil if necessary ###
Mask <- fieldMask(Test,
                  Blue = 1, Green = 2, Red = 3, RedEdge = 4, NIR = 5, # Layers position according to each sensor
                  index = "NDVI",
                  cropValue = 0.7,
                  cropAbove = FALSE)
```

<p align="center">
  <img src="https://raw.githubusercontent.com/filipematias23/images/master/readme/BF_7.jpg" width="80%" height="80%">
</p>

<br />

```r
################################################
### Extracting data for 14 flights in a loop ###
################################################

DataTotal<-NULL
for(i in 2:length(Field)){
  print(paste("===",Field[i],"==="))
  EX <- stack(paste("./5band/",Field[i],sep = ""))
  EX.1 <- fieldMask(EX,plot = F)
  EX.I <- fieldIndex(mosaic = EX,Red = 3,Green = 2,Blue = 1,RedEdge = 4,NIR = 5,
                     index = c("NGRDI","BGI","GLI","NDVI","NDRE","CIG","CIRE"),plot = F)
  crs(Shapefile)<-crs(EX.I)
  EX.I<- fieldInfo(mosaic = EX.I,
                   fieldShape = Shapefile,
                   buffer = -0.05,
                   n.core = 3)
  # Canopy Cover
  EX.I<-fieldArea(mosaic = EX.1$mask, 
                  fieldShape = EX.I$fieldShape,plot = F)
  # EPH
  DSM0 <- stack(paste("./DSM/",Field_DSM[1],sep = ""))
  DSM1 <- stack(paste("./DSM/",Field_DSM[i],sep = ""))
  
  # Canopy Height Model (CHM):
  DSM0 <- resample(DSM0, DSM1)
  CHM <- DSM1-DSM0
  CHM <- fieldMask(CHM, mask = Mask$mask, plot=F)
  CHM <- CHM$newMosaic
  
  # Extracting the estimate plant height average (EPH):
  EPH <- fieldInfo(CHM, fieldShape = EX.I$fieldShape, fun = "quantile",n.core = 3,plot = F)
  DataTotal1<-data.frame(Date=Field[i],EPH$fieldShape@data)
  EPH.1090.A<-extract(x = CHM, y = EPH$fieldShape)
  EPH.1090<-do.call(rbind,lapply(EPH.1090.A, quantile, probs = c(0.1,0.9), na.rm=TRUE))
  DataTotal1$'Height_10'<-EPH.1090[,1]
  DataTotal1$'Height_90'<-EPH.1090[,2]
  
  # Data Table:
  DataTotal<-rbind(DataTotal,DataTotal1)
  
  # Making plots:
  fieldPlot(EX.I$fieldShape,
            mosaic = EX,
            color = c("red","green"),
            #min.lim = 0.45,max.lim = 0.75,
            fieldAttribute = "NDVI")
}

### Correcting Names ###
Data.names<-gsub("layer",'Height_0',colnames(DataTotal))
Data.names<-gsub("NA..1",'Height_50',Data.names)
Data.names<-gsub("NA..2",'Height_75',Data.names)
Data.names<-gsub("NA..3",'Height_100',Data.names)
Data.names<-gsub("\\NA.",'Height_25',Data.names)
Data.names<-gsub("objArea",'Canopy',Data.names)
colnames(DataTotal)<-Data.names
DataTotal<-subset(DataTotal, select = -ID.1)
DataTotal$DAP<-as.numeric(do.call(rbind,strsplit(DataTotal$Date,split = "_"))[,4])
head(DataTotal)

### Saving extracted data in a .CSV ###
write.csv(DataTotal,"DataTotal.csv",row.names = F,col.names = T)
# DataTotal<-read.csv("DataTotal.csv",header = T)
```

<p align="center">
  <img src="https://raw.githubusercontent.com/filipematias23/images/master/readme/BF_8.jpg">
</p>

<br />

```r
####################
### Heritability ###
####################

### Agronomical Traits ###

Data <- read.csv("EX_DATA.csv",header = T,fileEncoding="UTF-8-BOM")
Data$RANGE<-as.factor(Data$RANGE)
Data$ROW<-as.factor(Data$ROW)
Data$NAME<-as.factor(Data$NAME)

Trait<-c("MAT_DAY","HT","DH","LODG","YLD")
H2.AG<-NULL
for(t in 1:length(Trait)){
  Data1<-droplevels(Data[!is.na(Data[,colnames(Data)==Trait[t]]),])
  mod<-lmer(eval(parse(text = paste(Trait[t],"~RANGE+ROW+(1|NAME)",sep=""))),data = Data1)
  Var1<-as.data.frame(VarCorr(mod))$vcov
  names(Var1)<-as.data.frame(VarCorr(mod))$grp
  H2.AG<-rbind(H2.AG,data.frame(Trait=Trait[t],
                          H2=round(c(Var1[1]/sum(Var1[1],
                                           Var1[2]/2 #2 replicates
                          )),3)
                          ))
  mod<-lm(eval(parse(text = paste(Trait[t],"~RANGE+ROW+NAME",sep=""))),data = Data1)
  Adj.Mean<-emmeans(mod, ~ NAME)
  if(t==1){
    Pheno.AG<-as.data.frame(Adj.Mean)[,c(1,2)]
  }
  if(t!=1){
    Pheno.AG<-merge(Pheno.AG,as.data.frame(Adj.Mean)[,c(1,2)],by="NAME")
  }
}
colnames(Pheno.AG)<-c("NAME",Trait)
head(Pheno.AG)

ggplot(data = H2.AG, 
       aes(x = Trait,
           y = H2*100,
           fill=as.factor(Trait))) +
  geom_bar(stat="identity", position = "dodge") +
  scale_fill_grey(start=0.2, end=0.8)+
  ylim(c(0,100))+
  labs(y="H2 (%)",
       x="", 
       fill="Agronomical Traits") +
  geom_text(aes(label=paste(H2*100,"%")),size=5, position=position_dodge(width=0.9), vjust=-0.25)+
  theme_bw()+
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.text = element_text(color="black",size=18),
        legend.title = element_text(color="black",size=18),
        axis.text.y = element_text(color="black",size=18),
        axis.title = element_text(color="black",size=18),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        strip.text = element_text(color="black",size=18),
        strip.background = element_rect(fill="white")) 
```




```r
### UAV Traits ###

DataTotal<-read.csv("DataTotal.csv",header = T)
DataTotal$RANGE<-as.factor(DataTotal$RANGE)
DataTotal$ROW<-as.factor(DataTotal$ROW)
DataTotal$NAME<-as.factor(DataTotal$NAME)

DAP<-unique(DataTotal$DAP)
Trait.UAV<-c("Blue","Green","Red","RedEdge","NIR", # Single Bands
             "NGRDI","BGI","GLI","NDVI","NDRE","CIG","CIRE", # Vegetation indices
             "Canopy", # Canopy cover
             "Height_0","Height_10","Height_25","Height_50","Height_75","Height_90","Height_100" #Estimated Plant Height
)

H2.UAV<-NULL
Pheno.UAV<-list()
for(t1 in 1:length(DAP)){
  Data<-droplevels(DataTotal[DataTotal$DAP==DAP[t1],])
  for(t2 in 1:length(Trait.UAV)){
    mod<-lmer(eval(parse(text = paste(Trait.UAV[t2]," ~ RANGE+ROW+(1|NAME)",sep=""))),data = Data)
    H2.UAV<-rbind(H2.UAV,cbind(DAP=DAP[t1],
                               Trait=Trait.UAV[t2],
                               H2=round(as.data.frame(VarCorr(mod))$vcov[1]/sum(as.data.frame(VarCorr(mod))$vcov[1],
                                                                                as.data.frame(VarCorr(mod))$vcov[2]/2),3)))
    mod<-lm(eval(parse(text = paste(Trait.UAV[t2],"~RANGE+ROW+NAME",sep=""))),data = Data)
    Adj.Mean<-emmeans(mod, ~ NAME)
    if(t2==1){
      Pheno.UAV.1<-as.data.frame(Adj.Mean)[,c(1,2)]
      colnames(Pheno.UAV.1)<-c("NAME",Trait.UAV[1:t2])
    }
    if(t2!=1){
      Pheno.UAV.1<-merge(Pheno.UAV.1,as.data.frame(Adj.Mean)[,c(1,2)],by="NAME")
      colnames(Pheno.UAV.1)<-c("NAME",Trait.UAV[1:t2])
      }
  }
  Pheno.UAV[[t1]]<-Pheno.UAV.1
}
names(Pheno.UAV)<-DAP

H2.UAV<-as.data.frame(H2.UAV)
H2.UAV$H2<-as.numeric(as.character(H2.UAV$H2))
H2.UAV$DAP<-as.numeric(as.character(H2.UAV$DAP))
H2.UAV$Trait<-factor(H2.UAV$Trait,levels = Trait.UAV)

ggplot(data = H2.UAV, 
       aes(x = DAP,
           y = H2*100)) +
  facet_wrap(~Trait)+
  geom_line(size=1) +
  geom_point(size=2)+
  ylim(c(0,100))+
  labs(y="H2 (%)",
       x="Days After Planting (DAP)") +
  theme_bw()+
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.text = element_text(color="black",size=18),
        legend.title = element_text(color="black",size=18),
        axis.text.y = element_text(color="black",size=10),
        axis.title = element_text(color="black",size=18),
        axis.text.x = element_text(color="black",size=10),
        strip.text = element_text(color="black",size=18),
        strip.background = element_rect(fill="white")) 
```




```r 
##################################
### Area Under the Curve (AUC) ###
##################################

### Choosing some DAP to investigate: ###

unique(DataTotal$DAP)[order(unique(DataTotal$DAP))]
#DAP<-c(16,19,25,31,34,41,46,49,52,55,67,73,74,80) 
DAP<-c(19,34,46,55,73,80)
Data<-DataTotal[DataTotal$DAP%in%DAP,]

### Data visualization - NDVI ###

ggplot(Data, 
       aes(x = NDVI,
           fill=as.factor(DAP))) +
  geom_density(alpha=.5,position = 'identity') +
  facet_wrap(~DAP, ncol = 1)+
  scale_fill_grey(start=1, end=0)+
  labs(y="#genotypes",
       x="NDVI", 
       fill="DAP") +
  theme_bw()+
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.text = element_text(color="black",size=18),
        legend.title = element_text(color="black",size=18),
        axis.text.y = element_text(color="black",size=14),
        axis.title = element_text(color="black",size=18),
        axis.text.x = element_text(color="black",size=14),
        strip.text = element_text(color="black",size=18),
        strip.background = element_rect(fill="white")) 
```



```r
### Visualization AUC ###

library(reshape2)
library(plyr)
library(DescTools)

unique(Data$NAME) # "FALLER","SYINGMAR","NDVITPRO"

Data1<-Data[as.character(Data$NAME)%in%c("FALLER","SYINGMAR","NDVITPRO"),]

Data2<-ddply(Data1,NAME~DAP,summarise,NDVI=mean(NDVI))

ggplot(data=Data2, aes(x=as.numeric(DAP), y= NDVI, col= NAME, group=NAME)) +
  geom_point(size=6)+
  geom_line(size=1.2) +
  scale_color_grey(start=0.8, end=0.2)+
  labs(x="Days After Planting (DAP)", fill="", col="")+
  theme_linedraw()+
  theme_bw()+
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.text = element_text(color="black",size=18),
        legend.title = element_text(color="black",size=18),
        axis.text.y = element_text(color="black",size=18),
        axis.title = element_text(color="black",size=18),
        axis.text.x = element_text(color="black",size=18),
        strip.text = element_text(color="black",size=18),
        strip.background = element_rect(fill="white")) 
```



```r
### AUC - NDVI ###

Trait<-c("NGRDI","NDRE","CIRE","Canopy","Height_50","Height_90") 

DataAUC<-fieldAUC(data = Data,
                  trait = Trait,
                  keep.columns = colnames(Data)[2:27],
                  frame = "long")

DataAUC$AUC<-as.numeric(DataAUC$AUC)
DataAUC$TRAIT<-factor(DataAUC$TRAIT,levels = Trait)
DataAUC$NAME<-as.factor(DataAUC$NAME)
DataAUC$RANGE<-as.factor(DataAUC$RANGE)
DataAUC$ROW<-as.factor(DataAUC$ROW)
# write.csv(DataAUC,"DataAUC.csv",row.names = F,col.names = T)

ggplot(DataAUC, aes(x = AUC)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.4,position = 'identity', fill="gray") +
  facet_wrap(~TRAIT,scales = "free")+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) 
```



```r
### Basic mixed model: ###

DataAUC$RANGE<-as.factor(DataAUC$RANGE)
DataAUC$ROW<-as.factor(DataAUC$ROW)
DataAUC$NAME<-as.factor(DataAUC$NAME)

H2.AUC<-NULL
for(t in 1:length(Trait)){
  Data1<-droplevels(DataAUC[as.character(DataAUC$TRAIT)==Trait[t],])
  mod<-lmer(AUC~RANGE+ROW+(1|NAME),data = Data1)
  Var1<-as.data.frame(VarCorr(mod))$vcov
  names(Var1)<-as.data.frame(VarCorr(mod))$grp
  H2.AUC<-rbind(H2.AUC,data.frame(Trait=Trait[t],
                                  H2=round(c(Var1[1]/sum(Var1[1],
                                                         Var1[2]/2 #2 replicates
                                  )),3)
  ))
  mod<-lm(AUC~RANGE+ROW+NAME,data = Data1)
  Adj.Mean<-emmeans(mod, ~ NAME)
  if(t==1){
    Pheno.AUC<-as.data.frame(Adj.Mean)[,c(1,2)]
    colnames(Pheno.AUC)<-c("NAME",Trait[1:t])
  }
  if(t!=1){
    Pheno.AUC<-merge(Pheno.AUC,as.data.frame(Adj.Mean)[,c(1,2)],by="NAME")
    colnames(Pheno.AUC)<-c("NAME",Trait[1:t])
    }
}
head(Pheno.AUC)

H2.AUC<-as.data.frame(H2.AUC)
H2.AUC$Trait<-factor(H2.AUC$Trait,Trait)
H2.AUC$H2<-as.numeric(H2.AUC$H2)

ggplot(data = H2.AUC, 
       aes(x = Trait,
           y = H2*100,
           fill=as.factor(Trait))) +
  geom_bar(stat="identity", position = "dodge") +
  scale_fill_grey(start=0.2, end=0.8)+
  ylim(c(0,100))+
  labs(y="H2 (%)",
       x="Trait", 
       fill="UAV_AUC") +
  geom_text(aes(label=paste(H2*100,"%")), position=position_dodge(width=0.9), vjust=-0.25)+
  theme_bw() 
```


```r
###########
### PCA ###
###########

library(factoextra)
Pheno.PCA<-merge(Pheno.AG,Pheno.AUC,by="NAME")
Pheno.PCA.1<-Pheno.PCA[,c("MAT_DAY","HT","DH","LODG","YLD",
                          "NDRE","Canopy","Height_90")]
rownames(Pheno.PCA.1)<-Pheno.PCA$NAME
Pheno.PCA.2 <- prcomp(Pheno.PCA.1, center = TRUE, scale = TRUE)

checks<-c("NDSW0932","NDSW14098","NDVITPRO","SYINGMAR","ALPINE","BARLOW","ELGIN-ND","FALLER","GLENN","MAX")
groups <- as.character(Pheno.PCA$NAME)
groups[groups%in%checks]<-"Checks"
groups[groups!="Checks"]<-"Lines"
groups.text <- as.character(Pheno.PCA$NAME)
groups.text[!groups.text%in%checks]<-""

require("ggrepel")

fviz_pca_biplot(Pheno.PCA.2,
             col.ind = groups, # color by groups
             legend.title = "",
             palette = c("black", "gray55"),
             repel = TRUE,
             geom.ind = "point",mean.point=F,
             pointshape = 21,
             pointsize = 3,
             fill.ind = groups
             )+
  geom_text_repel(aes(label = groups.text),size = 3.5)+
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.text = element_text(color="black",size=18),
        legend.title = element_text(color="black",size=18),
        axis.text.y = element_text(color="black",size=14),
        axis.title = element_text(color="black",size=18),
        axis.text.x = element_text(color="black",size=14),
        strip.text = element_text(color="black",size=18),
        strip.background = element_rect(fill="white"))
```




```r
###################
### Correlation ###
###################

library(agricolae)
library(corrplot)

### Specific flight ###

# 31 (MAT_DAY)
# 46 (YLD)
# 49 (DH)
# 74 (LODG)

Pheno.UAV.2<-Pheno.UAV$`74`[,c("NAME","NGRDI", "NDRE", "CIRE","Canopy","Height_50","Height_90")]
Pheno.COR<-merge(Pheno.AG,Pheno.UAV.2,by="NAME")
Pheno.COR.1<-scale(Pheno.COR[,-1],scale = T)
rownames(Pheno.COR.1)<-Pheno.COR[,1]
r<-correlation(Pheno.COR.1)
r$correlation
round(r$pvalue,2)
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(r$correlation, 
         p.mat = r$pvalue,
         sig.level = 0.05, # Level of significance
         method="color", col=col(200),  
         type="upper", order="hclust",addCoef.col = "black", 
         tl.col="black", tl.srt=45, 
         insig = "blank", 
         diag=FALSE)

### AUC ###
Pheno.COR<-merge(Pheno.AG,Pheno.AUC,by="NAME")
Pheno.COR.1<-scale(Pheno.COR[,-1],scale = T)
rownames(Pheno.COR.1)<-Pheno.COR[,1]
r<-correlation(Pheno.COR.1)
r$correlation
round(r$pvalue,2)
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(r$correlation, 
         p.mat = r$pvalue,
         sig.level = 0.1,
         method="color", col=col(200),  
         type="upper", order="hclust",addCoef.col = "black", 
         tl.col="black", tl.srt=45, 
         insig = "blank", 
         diag=FALSE)
```





```r
###################
### Heading Day ###
###################

### Simple Regression ###
DAP<-c(41,49,52)
Data<-droplevels(DataTotal[as.numeric(DataTotal$DAP)%in%DAP,])
Data[is.na(as.character(Data$DAP)),]
Data<-Data[,c("DAP","DH","Height_90")] #Other options: c("CIG","CIRE","Canopy","NGRDI","BGI","GLI","NDVI","NDRE","Height_50","Height_75","Height_90")

Data$DAP<-as.factor(Data$DAP)

Data.1<-melt(Data,
             value.name = "Index",
             measure.vars = c("Height_90"))

Data.2<-melt(Data.1,
             value.name = "Days to Heading",
             measure.vars = c("DH"))

colnames(Data.2)<-c("DAP","Index","Index.var","Trait","Trait.var")

Data.2$DAP<-as.factor(Data.2$DAP)
Data.2$Index<-as.factor(Data.2$Index)
Data.2$Trait<-as.factor(Data.2$Trait)
Data.2$Index.var<-as.numeric(as.character(Data.2$Index.var))
Data.2$Trait.var<-as.numeric(as.character(Data.2$Trait.var))

ggplot(data = Data.2,
       aes(x = Index.var,
           y = Trait.var,
           colour=Index)) +
  facet_grid(DAP~Trait, scales = "free_y")+
  geom_smooth(method=lm) +
  geom_point(size = 2) +
labs(y="Days to Heading (day of the year)",
     x="Index") +
  theme_bw()+
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.text = element_text(color="black",size=18),
        legend.title = element_blank(),
        axis.text.y = element_text(color="black",size=10),
        axis.title = element_text(color="black",size=18),
        axis.text.x = element_text(color="black",size=10),
        strip.text = element_text(color="black",size=18),
        strip.background = element_rect(fill="white"))
```





```r
################
### Maturity ###
################

### Simple Regression ###
DAP<-c(49,52,67,74,80)
Data<-DataTotal[as.numeric(DataTotal$DAP)%in%DAP,]
Data[is.na(as.character(Data$DAP)),]
Data<-Data[,c("DAP","MAT_DAY","CIG","CIRE")] 

Data$DAP<-as.factor(Data$DAP)

Data.1<-melt(Data,
             value.name = "Index",
             measure.vars = c("CIG","CIRE"))

Data.2<-melt(Data.1,
             value.name = "Trait",
             measure.vars = c("MAT_DAY"))

colnames(Data.2)<-c("DAP","Index","Index.var","Trait","Trait.var")

Data.2$DAP<-as.factor(Data.2$DAP)
Data.2$Index<-as.factor(Data.2$Index)
Data.2$Trait<-as.factor(Data.2$Trait)
Data.2$Index.var<-as.numeric(as.character(Data.2$Index.var))
Data.2$Trait.var<-as.numeric(as.character(Data.2$Trait.var))

ggplot(data = Data.2, 
       aes(x = Index.var,
           y = Trait.var,
           colour=Index)) + 
  facet_grid(DAP~Trait, scales = "free")+
  geom_smooth(method=lm) + 
  geom_point(size = 2) +
  labs(y="Maturity (day of the year)",
       x="Index") +
  theme_bw()+
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.text = element_text(color="black",size=18),
        legend.title = element_blank(),
        axis.text.y = element_text(color="black",size=10),
        axis.title = element_text(color="black",size=18),
        axis.text.x = element_text(color="black",size=10),
        strip.text = element_text(color="black",size=18),
        strip.background = element_rect(fill="white")) 
```



```r
###############
### Lodging ###
###############

### (FB) Flight before: 07/19/2021 (73 DAP) ###
### (FA1) Flight after_1: 07/20/2021 (74 DAP) ###
### (FA2) Flight after_2: 07/26/2021 (80 DAP) ###

### Simple Regression ###
DAP<-c(73,74,80)
Data<-DataTotal[as.numeric(DataTotal$DAP)%in%DAP,]
Data[is.na(as.character(Data$DAP)),]
Data<-Data[,c("DAP","LODG","Height_90")] #Other options: 

Data$DAP<-as.factor(Data$DAP)

Data.1<-melt(Data,
             value.name = "Index",
             measure.vars = c("Height_90"))

Data.2<-melt(Data.1,
             value.name = "Trait",
             measure.vars = c("LODG"))

colnames(Data.2)<-c("DAP","Index","Index.var","Trait","Trait.var")

Data.2$DAP<-as.factor(Data.2$DAP)
Data.2$Index<-as.factor(Data.2$Index)
Data.2$Trait<-as.factor(Data.2$Trait)
Data.2$Index.var<-as.numeric(as.character(Data.2$Index.var))
Data.2$Trait.var<-as.factor(as.character(Data.2$Trait.var))

ggplot(data = Data.2, 
       aes(y = Index.var,
           x = Trait.var,
           fill=Index)) + 
  facet_grid(DAP~Trait, scales = "free")+
  geom_boxplot() +
  labs(y="Estimate Plant Height",
       x="Lodging (Score)") +
  scale_fill_grey(start=0.8, end=0.2)+
  theme_bw()+
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.text = element_text(color="black",size=18),
        legend.title = element_blank(),
        axis.text.y = element_text(color="black",size=10),
        axis.title = element_text(color="black",size=18),
        axis.text.x = element_text(color="black",size=10),
        strip.text = element_text(color="black",size=18),
        strip.background = element_rect(fill="white")) 
```




```r
### Lodging Visualization ###

library(RStoolbox)
library(gridExtra)

###########
### RGB ###
###########

# 49 DAP #
lodg.52DAP <- stack("./Lodging/RGB/10_DAP_52_2021_Casselton_YT_06-28_rgb.tif")
p1<-ggRGB(lodg.52DAP, r=1, g=2, b=3, stretch = 'lin') +
  theme(axis.text.y = element_text(angle = 45),
        axis.title = element_blank())

# 74 DAP #
lodg.74DAP <- stack("./Lodging/RGB/15_DAP_74_2021_Casselton_YT_07-20_rgb.tif")
p2<-ggRGB(lodg.74DAP, r=1, g=2, b=3, stretch = 'lin') +
  theme(axis.text.y = element_text(angle = 45),
        axis.title = element_blank())

grid.arrange(p1, p2, ncol=2)
```




```r
###########
### DSM ###
###########

dev.off()
par(mfrow=c(1,2))

# 49 DAP #
lodg.52DAP.dsm <- stack("./Lodging/DSM/10_DAP_52_2021_Casselton_YT_06-28_dem.tif")
plot(lodg.52DAP.dsm)

# 74 DAP #
lodg.74DAP.dsm <- stack("./Lodging/DSM/15_DAP_74_2021_Casselton_YT_07-20_dem.tif")
plot(lodg.74DAP.dsm)
```



```r
#####################
### Drawing Lines ###
#####################

# Observing EPH profile for 2 rows:
Draw <- fieldDraw(mosaic = lodg.74DAP.dsm,
                     ndraw = 1)

# Profile plot:
dev.off()
par(mfrow=c(1,2))

# Profile plot:
plot(x = Draw$drawData$x, 
     y = Draw$drawData[,4], 
     type="l", col="red",lwd=1,
     #ylim=c(261,262),
     xlab="Distance (m)", 
     ylab="EPH (m)")

#RGB plot:
plotRGB(lodg.74DAP)
lines(Draw$drawData$x,Draw$drawData$y, type="l", col="red",lwd=2)
par(mfrow=c(1,1))
```




```r
################
### Modeling ###
################

library(car)
dev.off()

### Covariates in the model ###

# No covariate: 
Data <- read.csv("EX_DATA.csv",header = T,fileEncoding="UTF-8-BOM")
Data$RANGE<-as.factor(Data$RANGE)
Data$ROW<-as.factor(Data$ROW)
Data$NAME<-as.factor(Data$NAME)
mod<-lmer(YLD~RANGE+ROW+(1|NAME),data = Data)
(Yield.AIC<-AIC(mod))
qqPlot(residuals(mod))

# Single flight co-variate: 
Trait<-c("NGRDI","NDRE","CIRE","Canopy","Height_50","Height_90")
Data.AIC<-NULL
for(i in 1:length(Trait)){
DataTotal<-read.csv("DataTotal.csv",header = T)
DAP<-c(55)
Data<-DataTotal[DataTotal$DAP%in%DAP,]
Data$RANGE<-as.factor(Data$RANGE)
Data$ROW<-as.factor(Data$ROW)
Data$NAME<-as.factor(Data$NAME)
mod<-lmer(eval(parse(text = paste("YLD~",Trait[i],"+RANGE+ROW+(1|NAME)",sep=""))),data = Data)
Data.AIC<-rbind(Data.AIC,cbind(Trait=Trait[i],AIC=AIC(mod), Model="55DAP"))
qqPlot(residuals(mod))

# AUC co-variate: 
DataAUC<-read.csv("DataAUC.csv",header = T)
Data<-DataAUC[as.character(DataAUC$TRAIT)%in%Trait[i],]
Data$RANGE<-as.factor(Data$RANGE)
Data$ROW<-as.factor(Data$ROW)
Data$NAME<-as.factor(Data$NAME)
mod<-lmer(YLD~AUC+RANGE+ROW+(1|NAME),data = Data)
Data.AIC<-rbind(Data.AIC,cbind(Trait=Trait[i],AIC=AIC(mod), Model="AUC"))
qqPlot(residuals(mod))
}

Data.AIC<-as.data.frame(Data.AIC)
Data.AIC$AIC<-as.numeric(Data.AIC$AIC)

ggplot(data = Data.AIC, 
       aes(x = Trait,
           y = AIC,
           fill=as.factor(Trait))) +
  geom_bar(stat="identity", position = "dodge") +
  facet_wrap(~Model,scales = "free")+
  geom_hline(yintercept = Yield.AIC, col="red",linetype = "dashed", size=1)+
  scale_fill_grey(start=0.2, end=0.8)+
  labs(y="AIC",
       x="", 
       fill="UAV Traits") +
  theme_bw()+
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.text = element_text(color="black",size=18),
        legend.title = element_text(color="black",size=18),
        axis.text.y = element_text(color="black",size=18),
        axis.title = element_text(color="black",size=18),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        strip.text = element_text(color="black",size=18),
        strip.background = element_rect(fill="white")) 

```




```r
### Main trait selection (indirect selection) ###

# Yield selection based on BLUP (observed data): 
Data <- read.csv("EX_DATA.csv",header = T,fileEncoding="UTF-8-BOM")
Data$RANGE<-as.factor(Data$RANGE)
Data$ROW<-as.factor(Data$ROW)
Data$NAME<-as.factor(Data$NAME)
mod<-lmer(YLD~RANGE+ROW+(1|NAME),data = Data)
BLUP.Pheno<-as.matrix(ranef(mod)$NAME)
Sel.Pheno<-rownames(BLUP.Pheno)[order(BLUP.Pheno,decreasing = T)]

# Yield selection based on BLUP (UAV data): 
DataAUC<-read.csv("DataAUC.csv",header = T)
Trait<-c("NGRDI","NDRE","CIRE","Canopy","Height_50","Height_90")
Data.SC<-NULL
for(i in 1:length(Trait)){
  Data<-droplevels(DataAUC[as.character(DataAUC$TRAIT)%in%Trait[i],])
  Data$RANGE<-as.factor(Data$RANGE)
  Data$ROW<-as.factor(Data$ROW)
  Data$NAME<-as.factor(Data$NAME)
  mod<-lmer(AUC~RANGE+ROW+(1|NAME),data = Data)
  BLUP.AUC<-as.matrix(ranef(mod)$NAME)
  Sel.AUC<-rownames(BLUP.AUC)[order(BLUP.AUC,decreasing = T)]
  
  # Selection coincidence:
  n.sel<-round(length(Sel.Pheno)*0.25,0)
  Data.SC<-rbind(Data.SC,cbind(Trait=Trait[i],SC=sum(Sel.Pheno[1:n.sel]%in%Sel.AUC[1:n.sel])/n.sel))
}

Data.SC<-as.data.frame(Data.SC)
Data.SC$SC<-as.numeric(Data.SC$SC)

ggplot(data = Data.SC, 
       aes(x = Trait,
           y = SC*100,
           fill=as.factor(Trait))) +
  geom_bar(stat="identity", position = "dodge") +
  scale_fill_grey(start=0.2, end=0.8)+
  ylim(c(0,100))+
  labs(y="Indirect selection coincidence (%)",
       x="", 
       fill="UAV Traits") +
  geom_text(aes(label=paste(round(SC*100,2),"%")),size=5, position=position_dodge(width=0.9), vjust=-0.25)+
  theme_bw()+
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.text = element_text(color="black",size=18),
        legend.title = element_text(color="black",size=18),
        axis.text.y = element_text(color="black",size=18),
        axis.title = element_text(color="black",size=18),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        strip.text = element_text(color="black",size=18),
        strip.background = element_rect(fill="white")) 

```

