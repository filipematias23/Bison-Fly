# [Bison-Fly](https://github.com/filipematias23/Bison-Fly): UAV pipeline at NDSU Spring Wheat Breeding Program

The **Bison-Fly** tutorial is the *NDSU Spring Wheat UAV Pipeline* developed in partnership with *Drone2Phenome* ([D2P](https://www.ag2pi.org/)). In this pipeline we present step-by-step how we leverage UAV data in our breeding program. This is open source R code that anyone may use and adopt for different crops and situations. We hope this pipeline helps you with your research, and all suggestions to improve it are welcome. Help us to build this platform in partnership with the D2P community. 

<p align="center">
  <img src="https://raw.githubusercontent.com/filipematias23/images/master/readme/BisonFly.jpg" width="40%" height="40%">
</p>

<div id="menu" />

## Resources
  
   * [Introduction](#intro)
   * [Image analysis in R](#r1)
   * [1. Agronomic Traits](#p1)
   * [2. UAV Traits](#p2)
   * [3. Area Under the Curve (AUC)](#p3)
   * [4. Principal component analysis (PCA)](#p4)
   * [5. Correlation (r)](#p5)
   * [6. Heading Day](#p6)
   * [7. Maturity](#p7)
   * [8. Lodging](#p8)
   * [9. Statistical Applications](#p9)
   * [Contact](#pC)

<br />

<div id="intro" />

## Introduction 

> Welcome to the **NDSU Spring Wheat Breeding Program** led by [Dr. Andrew Green](https://www.ndsu.edu/agriculture/ag-home/directory/andrew-green). In this tutorial, we present how we leverage UAV data in our breeding program. This pipeline is a compilation of [R]( https://www.r-project.org/) code and functions from different packages to implement high-throughput phenotyping methods to understand spring wheat traits.

<p align="center">
  <img src="https://raw.githubusercontent.com/filipematias23/images/master/readme/BF_1a.jpg">
</p>

<br />

> Our pipeline starts by having at least one authorized drone pilot certified by the Federal Aviation Administration ([FAA]( https://www.faa.gov/)). This person needs to pass the Remote Pilot of Small Unmanned Aircraft Systems exam ([Part 107]( https://www.faa.gov/uas/commercial_operators/become_a_drone_pilot/)). The next step is getting the institution approval to collect data using UAVs. In the case of NDSU, any UAV flight related with research requires approval from the Research Operations Activities Office by filling out the [Flight Plan Form]( https://www.ndsu.edu/research/for_researchers/unmanned_aircraft_systems/). Receiving approval normally takes some time and needs to be planned ahead of the season for sufficient time before planting starts. 

<br />

> During 2021, our team used the drone [Inspire 2 from DJI]( https://www.dji.com/inspire-2) combined with the multispectral sensor [Sentera 6X]( https://sentera.com/data-capture/6x-multispectral/). Our phenotyping kit inclues three set of batteries that gives us autonomy for around one hour of flying time. The kit also has one [reflectance panel]( https://support.sentera.com/portal/en/kb/articles/capturing-your-refelectance-panel) for radiometric calibration. This step is important to transform digital numbers to radiance and then reflectance. In addition, this step allows comparing data from different days or locations, which is very common in plant breeding programs. 

<p align="center">
  <img src="https://raw.githubusercontent.com/filipematias23/images/master/readme/BF_2.jpg" width="70%" height="70%">
</p>

<br />

> The flight mission was planned using the software [FieldAgent from Sentera]( https://sentera.com/fieldagent/). This software is used to draw the trial using the desktop application before going to field. Once the flying area was established, the flight parameters were set for 75% of overlapping, 15 mph, and flying height of 200 feet. There are other great open source options for DJI drones as [Pix4Dcapture]( https://support.pix4d.com/hc/en-us/articles/202557269-Pix4Dcapture-Getting-Started) and [DroneDeploy]( https://www.dronedeploy.com/product/mobile/). 

<br />

> For georeferencing, we use white plastic bucket lids (5 gallons) as ground control points (GCPs). These lids are fixed in the ground with metal landscaping stakes in the corners of the field trial at the beginning of the planting season. We opted to paint the central part of the lid to facilitate finding where to georeference during the orthomosaic step. We also made holes to drain water from the lids. The GCPs remain in the field during the entire growing season and must be photographed alongside the plants for each flight or other image collection. GCPs enable a better overlap of the experimental plots and facilitate drawing grid plot polygons only once. Please find more information at this link: [https://github.com/OpenDroneMap/FIELDimageR#P5](https://github.com/OpenDroneMap/FIELDimageR#P5). 

<p align="center">
  <img src="https://raw.githubusercontent.com/filipematias23/images/master/readme/BF_3.jpg">
</p>

<br />

> Another simple way to find the corners on a field trial and draw the experimental plot polygons is to use four of the same plastic lids as extra points between borders and plots. This is an easy way to find where to click when using the function fieldShape() from [FIELDimageR]( https://github.com/OpenDroneMap/FIELDimageR#P5). Also, this method helps to identify where the boundaries of multiple experimental trials planted in sequence as highlighted in the example below from our breeding program.

<p align="center">
  <img src="https://raw.githubusercontent.com/filipematias23/images/master/readme/BF_4.jpg">
</p>

<br />

> In this tutorial we evaluate a population of 106 white wheat lines in a second year yield trial (breeding pipeline). A the total 15 flights were flown in 2021 at NDSU Casselton Seed Farm in North Dakota (46.88051067579893, -97.2370382539691). All genotypes had at least two replicates, some check lines had three replicates (e.g., ALPINE, BARLOW, ELGIN-ND, FALLER, GLENN, SYINGMAR, NDSW0932, NDSW14098). The population was evaluated for maturity day (MAT_DAY), plant height (HT), days to heading (DH), lodging (LODG), and yield (YLD). All flights were performed around noon on days without clouds. The flights dates were (1) 05/11, (2) 05/23, (3) 05/26, (4) 06/01, (5) 06/07, (6) 06/10, (7) 06/17, (8) 06/22, (9) 06/25, (10) 06/28, (11) 07/01, (12) 07/13, (13) 07/19, (14) 07/20, and (15) 07/26.  

<p align="center">
  <img src="https://raw.githubusercontent.com/filipematias23/images/master/readme/BF_5.jpg">
</p>

<br />

> Open source software we suggest for creating an orthomosaic is [OpenDroneMap]( https://www.opendronemap.org/) with the click and point interface [WebODM]( https://www.opendronemap.org/webodm/). A workshop tutorial video is available at [Phenome-Force YouTube Channel]( https://www.youtube.com/watch?v=U-bsA7QjzYE&t=4842s). More information is avaliable [HERE](https://github.com/OpenDroneMap/FIELDimageR/blob/master/README.md#P18).

<p align="center">
  <img src="https://raw.githubusercontent.com/filipematias23/images/master/readme/ODM_FIELDimageR_New.jpg" width="70%" height="70%">
</p>

<div id="r1" />

[Menu](#menu)

## Image analysis in R

> The first step is to **DOWNLOAD** the data used in this tutorial by clicking [**HERE**](https://drive.google.com/file/d/1_Uj3oaiSv31WpbGyyUCfKiKhAqQ_jPUN/view?usp=sharing). The folder has 15 Multispectral orthomosaics, 15 RGB orthomosaics, 15 DSM orthomosaics, 1 CSV file with agronomical field data, and 1 CSV file with the plots' ID map. 

<p align="center">
  <img src="https://raw.githubusercontent.com/filipematias23/images/master/readme/BF_0.jpg" width="50%" height="50%">
</p>

<br />

* Orthomosaics metadata is avaliable [**HERE**](https://docs.google.com/spreadsheets/d/1E2jmO3Su38PYQ-5FHZkYjH4VgSlFvvNa/edit?usp=sharing&ouid=108197715650347930652&rtpof=true&sd=true) 

<p align="center">
  <img src="https://raw.githubusercontent.com/filipematias23/images/master/readme/BF_meta.jpg">
</p>

<br />

* Hourly weather and soil data from 7-May-2021 to 7-August-2021 is available [**HERE**](https://drive.google.com/file/d/1gedPmoVb7LVdZBxmGrN1vllHQqZaalpc/view?usp=sharing). This data was collected in Prosper/ND (NDSU Station) at the closest station from Casselton/ND (experimental fields in this tutorial). 

<p align="center">
  <img src="https://raw.githubusercontent.com/filipematias23/images/master/readme/BF_weather.jpg">
</p>

<br />

> The following code is an example how to prepare the images to extract the UAV data with biological meaning for posterior breeding applications.

**Steps:**
* Necessary packages
* Uploading orthomosaics and unzip (save all files in the same directory)
* Drawing polygons grid ("Shapefile")
* Removing soil (Evaluate canopy coverage)
* Extracting data for 14 flights in a loop (UAV traits: calculating vegetation indices and plant height)

<br />

> In order to run this tutorial, you need to install [R](https://www.r-project.org/) and [RStudio](https://www.rstudio.com/). For Windows users who have an R version higher than 4.0, you need to install [RTools](https://cran.r-project.org/bin/windows/Rtools/rtools40.html).

```r
##########################################
### Bison-Fly: Plant Breeding Pipeline ###
##########################################

### Install Packages ###
requiredPackages = c("devtools","terra","sf","ggplot2","DescTools","lme4","emmeans","reshape2","car","plyr","factoextra","ggrepel","agricolae","corrplot","RStoolbox","gridExtra",'mapview','stars','caret','mapedit',"dplyr","fields","leafem","leafsync","lwgeom","BiocManager","git2r","exactextractr")
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p,character.only = TRUE)
}

devtools::install_github("OpenDroneMap/FIELDimageR")
devtools::install_github("filipematias23/FIELDimageR.Extra")

### Necessary packages ###
library(FIELDimageR)
library(FIELDimageR.Extra)
library(terra)
library(sf)
library(ggplot2)
library(DescTools)
library(lme4)
library(emmeans)
library(reshape2)
library(car)
library(plyr)
library(factoextra)
library(ggrepel)
library(agricolae)
library(corrplot)
library(RStoolbox)
library(gridExtra)

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
Test <- rast(paste("./5band/",Field[4],sep = ""))
fieldView(Test)

```

> The plot shape file can be drawn by selecting at least four points at the corners of the experiment. The number of columns and rows must be informed. At this point the experimental borders can be eliminated as shwing in the example bellow. You can download the Grid-Polygons-Plot in a shapefile format for comparations with your method on [HERE](https://drive.google.com/file/d/13kDE2vc2tZe1oFkEM_qdGet70kWIE8zI/view?usp=sharing)

```r
##########################
### Plot Polygons Grid ###
##########################

### Field data ###
Data <- read.csv("EX_DATA.csv",header = T,fileEncoding="UTF-8-BOM")

### MAP ###
Map <- read.csv("EX_MAP.csv",header = F,fileEncoding="UTF-8-BOM")
Map
rotate <- function(x) t(apply(x, 2, rev))
Map<-rotate(rotate(Map))
Map

### Grid ###
Shapefile <- fieldShape_render(mosaic = Test,
                        ncols = 11,
                        nrows = 20,
                        fieldData = Data,
                        PlotID = "PLOT",
                        fieldMap = Map)
fieldView(mosaic = Test,
          fieldShape = Shapefile,
          type = 2,
          alpha = 0.2)

### Saving Shapefile ###
# write_sf(Shapefile, "Shapefile.shp")
# Shapefile <- read_sf("Shapefile.shp") # Reading the saved shapefile option 02.

### If you downloaded the example in the descriptions above use this code to read the shapefile:
# unzip("Shapefile.zip")
# Shapefile <- read_sf("./Shapefile/Shapefile.shp") # Reading the saved shapefile.
# Shapefile$PlotID<-Shapefile$PlotName
# plot(Shapefile, border="red")

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

> Extracting data from all field images at once in a loop. The code below takes time to run. If you don't want to run it you can just **DOWNLOAD** the *'DataTotal.csv'* with extracted data [HERE](https://drive.google.com/file/d/1t6ZHtbldBouJ675i8XJ1nMQm7ReFeaEd/view?usp=sharing).

<br />

```r
################################################
### Extracting data for 14 flights in a loop ###
################################################

DataTotal<-NULL
for(i in 2:length(Field)){
  print(paste("===",Field[i],"==="))
  EX <- rast(paste("./5band/",Field[i],sep = ""))
  EX.1 <- fieldMask(EX,plot = F)
  EX.I <- fieldIndex(mosaic = EX,Red = 3,Green = 2,Blue = 1,RedEdge = 4,NIR = 5,
                     index = c("NGRDI","BGI","GLI","NDVI","NDRE","CIG","CIRE"),plot = F)
  EX.I<- fieldInfo_extra(mosaic = EX.I,
                   fieldShape = Shapefile)

  # EPH
  DSM0 <- rast(paste("./DSM/",Field_DSM[1],sep = ""))
  DSM1 <- rast(paste("./DSM/",Field_DSM[i],sep = ""))
  
  # Canopy Height Model (CHM):
  # DSM0 <- resample(DSM0, DSM1)
  CHM <- fieldHeight(dsm_before = DSM0, dsm_after = DSM1)
  CHM <- fieldMask(CHM$height, mask = Mask$mask, plot=F)
  CHM <- CHM$newMosaic
  
  # Extracting the estimate plant height average (EPH):
  EPH <- exactextractr::exact_extract(x = CHM,
                                           y = st_as_sf(EX.I),
                                           fun = "quantile",
                                           quantiles=c(0,0.1,0.25,0.5,0.75,0.9,1))
  
  EX.I <- cbind(EX.I, EPH)

  # Canopy Cover
  EX.AI<-fieldArea(mosaic = EX.1$newMosaic, 
                  fieldShape = EX.I)
  EX.I<-merge.data.frame(EX.I,EX.AI[,1:3],by="PlotID")

  DataTotal1<-data.frame(Date=Field[i],EX.I)
  DataTotal1<- DataTotal1[,!colnames(DataTotal1)%in%c("geometry.y","geometry.x","geometry")]
  
  # Data Table:
  DataTotal<-rbind(DataTotal,DataTotal1)
  
  # Making plots:
  fieldView(mosaic = EX,
            fieldShape = st_as_sf(EX.I),
            plotCol = "NDVI_mean",
            col_grid = c("blue", "grey", "red"),
            type = 2,
            alpha_grid = 0.6)
  
}

### Correcting Names ###
Data.names<-gsub("q00",'Height_0',colnames(DataTotal))
Data.names<-gsub("q10",'Height_10',Data.names)
Data.names<-gsub("q25",'Height_25',Data.names)
Data.names<-gsub("q50",'Height_50',Data.names)
Data.names<-gsub("q75",'Height_75',Data.names)
Data.names<-gsub("q90",'Height_90',Data.names)
Data.names<-gsub("q100",'Height_100',Data.names)
Data.names<-gsub("AreaPercentage",'Canopy',Data.names)
Data.names<-gsub("_mean",'',Data.names)
colnames(DataTotal)<-Data.names
DataTotal$DAP<-as.numeric(do.call(rbind,strsplit(DataTotal$Date,split = "_"))[,4])
head(DataTotal)

### Saving extracted data in a .CSV ###
# write.csv(DataTotal,"DataTotal.csv",row.names = F,col.names = T)
# DataTotal<-read.csv("DataTotal.csv",header = T)

```

<p align="center">
  <img src="https://raw.githubusercontent.com/filipematias23/images/master/readme/BF_8.jpg">
</p>

<br />

[Menu](#menu)

<div id="p1" />

## Agronomic Traits

> Each species has a key set of important agronomic traits with economical potential that are evaluated to characterize populations and apply selection. The following code is an example how to calculate adjusted means and heritability in a simple and fast way. As a reminder, this is just an example that must be adapted for each experiment to account for different field designs. For instance the same model was evaluated twice, the first with genotyped with random effect to calculate heritability (*package [lme4]( https://cran.r-project.org/web/packages/lme4/index.html)*) and the second as fixed effect (*function lm*) to calculate the adjusted means using the *package [emmeans]( https://cran.r-project.org/web/packages/emmeans/)*. The adjusted means will be used for further statistical analysis in this tutorial as well as principal component analysis and yield prediction. 

```r
##########################
### Agronomical Traits ###
##########################

### Field data collected manually ###
Data <- read.csv("EX_DATA.csv",header = T,fileEncoding="UTF-8-BOM")
Data$RANGE<-as.factor(Data$RANGE)
Data$ROW<-as.factor(Data$ROW)
Data$NAME<-as.factor(Data$NAME)

### Mixed model: getting adjusted means and heritability (H2) ###
Trait<-c("MAT_DAY","HT","DH","LODG","YLD")
H2.AG<-NULL
for(t in 1:length(Trait)){
  Data1<-droplevels(Data[!is.na(Data[,colnames(Data)==Trait[t]]),])
  # mod<-lmer(eval(parse(text = paste(Trait[t],"~RANGE+ROW+(1|NAME)",sep=""))),data = Data1)
  mod<-lmer(eval(parse(text = paste(Trait[t],"~(1|NAME)",sep=""))),data = Data1)
  Var1<-as.data.frame(VarCorr(mod))$vcov
  names(Var1)<-as.data.frame(VarCorr(mod))$grp
  H2.AG<-rbind(H2.AG,data.frame(Trait=Trait[t],
                          H2=round(c(Var1[1]/sum(Var1[1],
                                           Var1[2]/2 #2 replicates
                          )),3)
                          ))
  # mod<-lm(eval(parse(text = paste(Trait[t],"~RANGE+ROW+NAME",sep=""))),data = Data1)
  mod<-lm(eval(parse(text = paste(Trait[t],"~NAME",sep=""))),data = Data1)
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

### Agronomical traits heritability ###
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

<p align="center">
  <img src="https://raw.githubusercontent.com/filipematias23/images/master/readme/BF_9a.jpeg" width="80%" height="80%">
</p>

<br />

[Menu](#menu)

<div id="p2" />

## UAV Traits

> In the same way presented above for agronomic traits, the code below was adapted to evaluate many UAV traits and calculate their heritability and adjusted means over time in a loop (Days After Planting - DAP). **DOWNLOAD** the *'DataTotal.csv'* with extracted data [HERE](https://drive.google.com/file/d/1t6ZHtbldBouJ675i8XJ1nMQm7ReFeaEd/view?usp=sharing).

```r
##################
### UAV Traits ###
##################

### UAV data extracted above (DataTotal) ###
DataTotal<-read.csv("DataTotal.csv",header = T)
DataTotal$RANGE<-as.factor(DataTotal$RANGE)
DataTotal$ROW<-as.factor(DataTotal$ROW)
DataTotal$NAME<-as.factor(DataTotal$NAME)

### Preparing information for running the statistic models below in a loop ###
DAP<-unique(DataTotal$DAP)
Trait.UAV<-c("Blue","Green","Red","RedEdge","NIR", # Single Bands
             "NGRDI","BGI","GLI","NDVI","NDRE","CIG","CIRE", # Vegetation indices
             "Canopy", # Canopy cover
             "Height_0","Height_10","Height_25","Height_50","Height_75","Height_90","Height_100" #Estimated Plant Height
)

### Mixed model: getting adjusted means and heritability (H2) ###
H2.UAV<-NULL
Pheno.UAV<-list()
for(t1 in 1:length(DAP)){
  Data<-droplevels(DataTotal[DataTotal$DAP==DAP[t1],])
  for(t2 in 1:length(Trait.UAV)){
    # mod<-lmer(eval(parse(text = paste(Trait.UAV[t2]," ~ RANGE+ROW+(1|NAME)",sep=""))),data = Data)
    mod<-lmer(eval(parse(text = paste(Trait.UAV[t2]," ~ (1|NAME)",sep=""))),data = Data)
    H2.UAV<-rbind(H2.UAV,cbind(DAP=DAP[t1],
                               Trait=Trait.UAV[t2],
                               H2=round(as.data.frame(VarCorr(mod))$vcov[1]/sum(as.data.frame(VarCorr(mod))$vcov[1],
                                                                                as.data.frame(VarCorr(mod))$vcov[2]/2),3)))
    # mod<-lm(eval(parse(text = paste(Trait.UAV[t2],"~RANGE+ROW+NAME",sep=""))),data = Data)
    mod<-lm(eval(parse(text = paste(Trait.UAV[t2],"~NAME",sep=""))),data = Data)
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

### UAV traits heritability ###
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

<p align="center">
  <img src="https://raw.githubusercontent.com/filipematias23/images/master/readme/BF_10a.jpeg" width="80%" height="80%">
</p>

<br />

[Menu](#menu)

<div id="p3" />

## Area Under the Curve (AUC)

> There are many interesting ways to use the UAV data for plant breeding. The first way is evaluating single flights for specific time point aims (e.g., evaluating maturity, plant height, or disease resistance). These traits normally occur in a specific moment during the season and can be evaluated with 1 or few flights. On the other hand, traits such as plant development, biomass, and yield can be evaluated throughout the entire season using UAV data. In this case, one good strategy is combining all flights in the same analysis. The first option is using a vector with DAP as fixed effect in the model to capture the trait performance over time. The second option is described below which is calculating the area under the curve (AUV) for one trait over the time. This is an interesting combined UAV-multitrait, because some genotypes have slow growing abilities in the beginning of the season but can perform well on the other growing stages. At the same point, genotypes with great performance in early stages can reduce competitivity at the end of the cycle. Using AUC is an interesting way to observe and compare these different biological paths and use only one general trait (e.g., growing performance) for applying selection.

```r 
##################################
### Area Under the Curve (AUC) ###
##################################

### Choosing some days after planting (DAP) to investigate ###
unique(DataTotal$DAP)[order(unique(DataTotal$DAP))] # Options: c(16,19,25,31,34,41,46,49,52,55,67,73,74,80) 
DAP<-c(19,34,46,55,73,80)
Data<-DataTotal[DataTotal$DAP%in%DAP,]

### Extracted data visualization per DAP - NDVI ###
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

<p align="center">
  <img src="https://raw.githubusercontent.com/filipematias23/images/master/readme/BF_11.jpg" width="80%" height="80%">
</p>

<br />

```r
#########################
### AUC Visualization ###
#########################

### Choosing genotypes to highlight ###
unique(Data$NAME) # "FALLER","GLENN","NDVITPRO"
Data1<-Data[as.character(Data$NAME)%in%c("FALLER","GLENN","NDVITPRO"),]
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

<p align="center">
  <img src="https://raw.githubusercontent.com/filipematias23/images/master/readme/BF_Curve_New.jpeg" width="80%" height="80%">
</p>

<br />

```r
#######################
### Calculating AUC ###
#######################

### Choosing some UAV traits ###
Trait<-c("NGRDI","NDRE","CIRE","Canopy","Height_50","Height_90") 
DataAUC<-fieldAUC(data = Data,
                  trait = Trait,
                  keep.columns = colnames(Data)[2:27],
                  frame = "long")
DataAUC$AUC<-as.numeric(as.character(DataAUC$AUC))
DataAUC$TRAIT<-factor(DataAUC$TRAIT,levels = Trait)
DataAUC$NAME<-as.factor(DataAUC$NAME)
DataAUC$RANGE<-as.factor(DataAUC$RANGE)
DataAUC$ROW<-as.factor(DataAUC$ROW)
write.csv(DataAUC,"DataAUC.csv",row.names = F,col.names = T)
ggplot(DataAUC, aes(x = AUC)) +
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.4,position = 'identity', fill="gray") +
  facet_wrap(~TRAIT,scales = "free")+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) 
```

<p align="center">
  <img src="https://raw.githubusercontent.com/filipematias23/images/master/readme/BF_13.jpg" width="80%" height="80%">
</p>

<br />

```r
####################
### AUC analysis ###
####################

DataAUC<-read.csv("DataAUC.csv",header = T)
DataAUC$RANGE<-as.factor(DataAUC$RANGE)
DataAUC$ROW<-as.factor(DataAUC$ROW)
DataAUC$NAME<-as.factor(DataAUC$NAME)

### Mixed model: getting adjusted means and heritability (H2) ###
H2.AUC<-NULL
for(t in 1:length(Trait)){
  Data1<-droplevels(DataAUC[as.character(DataAUC$TRAIT)==Trait[t],])
  # mod<-lmer(AUC~RANGE+ROW+(1|NAME),data = Data1)
  mod<-lmer(AUC~(1|NAME),data = Data1)
  Var1<-as.data.frame(VarCorr(mod))$vcov
  names(Var1)<-as.data.frame(VarCorr(mod))$grp
  H2.AUC<-rbind(H2.AUC,data.frame(Trait=Trait[t],
                                  H2=round(c(Var1[1]/sum(Var1[1],
                                                         Var1[2]/2 #2 replicates
                                  )),3)
  ))
  # mod<-lm(AUC~RANGE+ROW+NAME,data = Data1)
  mod<-lm(AUC~NAME,data = Data1)
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

### AUC traits heritability ###
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
       x="", 
       fill="UAV_AUC") +
  geom_text(aes(label=paste(H2*100,"%")), position=position_dodge(width=0.9), vjust=-0.25)+
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

<p align="center">
  <img src="https://raw.githubusercontent.com/filipematias23/images/master/readme/BF_14a.jpeg" width="80%" height="80%">
</p>

<br />

[Menu](#menu)

<div id="p4" />

## Principal component analysis (PCA)

> Principal component analysis (PCA) is an exploratory data analysis commonly used for dimensionality reduction by projecting data point variability into the first few principal components. In the plant breeding perspective, PCA analysis helps to find multivariate patterns in the population by using lower-dimensional data visualization while preserving as much as possible data's variation. In the example below, agronomic traits and AUC-UAV-traits were used to characterize the wheat breeding population and observe the existence of potential genotypes compared with the check varieties. It’s possible to observe a great potential for selection in this population that has genotypes with higher multivariate performance compared with the named check varieties.  

```r
##########################################
### Principal component analysis (PCA) ###
##########################################

### Merging Agro and UAV adjusted means ###
Pheno.PCA<-merge(Pheno.AG,Pheno.AUC,by="NAME")
Pheno.PCA.1<-Pheno.PCA[,c("MAT_DAY","HT","DH","LODG","YLD",
                          "NDRE","Canopy","Height_90")]
rownames(Pheno.PCA.1)<-Pheno.PCA$NAME

### PCA ###
Pheno.PCA.2 <- prcomp(Pheno.PCA.1, center = TRUE, scale = TRUE)

### Highlighting checks ###
checks<-c("NDSW0932","NDSW14098","NDVITPRO","CHECK_01","CHECK_02","BARLOW","ELGIN-ND","FALLER","GLENN","MAX")
groups <- as.character(Pheno.PCA$NAME)
groups[groups%in%checks]<-"Checks"
groups[groups!="Checks"]<-"Lines"
groups.text <- as.character(Pheno.PCA$NAME)
groups.text[!groups.text%in%checks]<-""

### PCA visualization ###
fviz_pca_biplot(Pheno.PCA.2,
             col.ind = groups, # color by groups
             legend.title = "",
             palette = c("black", "gray70"),
             repel = TRUE,
             geom.ind = "point",mean.point=F,
             pointsize = 3,
             fill.ind = groups
             )+
  geom_text_repel(aes(label = groups.text),size = 3.5)+
  scale_shape_manual(values=c(17,21))+
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.text = element_text(color="black",size=18),
        legend.title = element_text(color="black",size=18),
        axis.text.y = element_text(color="black",size=14),
        # axis.title = element_text(color="black",size=18),
        axis.text.x = element_text(color="black",size=14),
        strip.text = element_text(color="black",size=18),
        strip.background = element_rect(fill="white"))
```

<p align="center">
  <img src="https://raw.githubusercontent.com/filipematias23/images/master/readme/BF_PCA_New2.jpeg" width="80%" height="80%">
</p>

<br />

[Menu](#menu)

<div id="p5" />

## Correlation Analysis (r)

>  Pearson correlation coefficient is a measure of linear correlation between two traits. This value gives to the breeders an idea how traits can be direct or indirect related. This measurement must be used with caution- more complex statistical analyses should be used to understand the real genetic/genomic connection among traits before applying selection. However, for preliminary data evaluation this analysis gives a good idea with what potential agronomical traits the UAV traits can be connected. There are some different ways to make this analysis using single flights to observe specific agronomical traits or AUC-UAV-traits for complex whole growth cycle traits as yield.  

```r
###################
### Correlation ###
###################

### Specific flight ###
# 31 DAP (MAT_DAY)
# 46 DAP (YLD)
# 49 DAP (DH)
# 74 DAP (LODG) and 73 DAP (1 day before heavy storm)

### 73 DAP ###
Pheno.UAV.2<-Pheno.UAV$`73`[,c("NAME","NGRDI", "NDRE", "CIRE","Canopy","Height_50","Height_90")]
Pheno.COR<-merge(Pheno.AG,Pheno.UAV.2,by="NAME")
Pheno.COR.1<-scale(Pheno.COR[,-1],scale = T)
rownames(Pheno.COR.1)<-Pheno.COR[,1]

### r (73 DAP) ###
r<-correlation(Pheno.COR.1)
r$correlation
round(r$pvalue,2)
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(r$correlation, 
         p.mat = r$pvalue,
         sig.level = 0.05, # Level of significance 5%
         method="color", col=col(200),  
         type="upper", order="alphabet",addCoef.col = "black", 
         tl.col="black", tl.srt=45, 
         insig = "blank", 
         diag=FALSE)
         
### 74 DAP ###
Pheno.UAV.2<-Pheno.UAV$`74`[,c("NAME","NGRDI", "NDRE", "CIRE","Canopy","Height_50","Height_90")]
Pheno.COR<-merge(Pheno.AG,Pheno.UAV.2,by="NAME")
Pheno.COR.1<-scale(Pheno.COR[,-1],scale = T)
rownames(Pheno.COR.1)<-Pheno.COR[,1]

### r (74 DAP) ###
r<-correlation(Pheno.COR.1)
r$correlation
round(r$pvalue,2)
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(r$correlation, 
         p.mat = r$pvalue,
         sig.level = 0.05, # Level of significance 5%
         method="color", col=col(200),  
         type="upper", order="alphabet",addCoef.col = "black", 
         tl.col="black", tl.srt=45, 
         insig = "blank", 
         diag=FALSE)

### AUC ###
Pheno.AUC<-Pheno.AUC[,c("NAME","NGRDI", "NDRE", "CIRE","Canopy","Height_50","Height_90")]
Pheno.COR<-merge(Pheno.AG,Pheno.AUC,by="NAME")
Pheno.COR.1<-scale(Pheno.COR[,-1],scale = T)
rownames(Pheno.COR.1)<-Pheno.COR[,1]

### r (AUC) ###
r<-correlation(Pheno.COR.1)
r$correlation
round(r$pvalue,2)
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(r$correlation, 
         p.mat = r$pvalue,
         sig.level = 0.1,
         method="color", col=col(200),  
         type="upper", order="alphabet",addCoef.col = "black", 
         tl.col="black", tl.srt=45, 
         insig = "blank", 
         diag=FALSE)
```

<p align="center">
  <img src="https://raw.githubusercontent.com/filipematias23/images/master/readme/BF_16_New.jpg">
</p>

<br />

[Menu](#menu)

<div id="p6" />

## Heading Day

> Heading date in wheat is associated with the timing of the floral transition and is an important agronomic trait that affects crop production. Heading and anthesis are critical periods for management strategies to cope with important diseases such as Fusarium Head Blight (FHB). We measure heading day as the Days After Planting when 50% of the plants in the plot has a complete spike exposed as observed in the left plot in the image below. In the example below, only the simple regression at 49 DAP showed some connection between heading date and estimate plant height. Also, heading date in this population is a high herdable trait (around 90% as shown in the H2 for agronomical traits above) compared with any H2 observed for UAV traits. In this case, using the simple regression approach was not a good strategy to evaluate heading date.  

<p align="center">
  <img src="https://raw.githubusercontent.com/filipematias23/images/master/readme/BF_HD.jpg">
</p>

<br />

```r
########################
### Heading Day (DH) ###
########################

### Choosing flights around the DH ###
DAP<-c(41,49,52)
Data<-droplevels(DataTotal[as.numeric(DataTotal$DAP)%in%DAP,])

### Choosing UAV traits to compare with DH ###
Data<-Data[,c("DAP","DH","Height_90")] # Other options: c("CIG","CIRE","Canopy","NGRDI","BGI","GLI","NDVI","NDRE","Height_50","Height_75","Height_90")
Data$DAP<-as.factor(Data$DAP)

### Preparing data to make plots ###
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

### Simple regression visualization ###
ggplot(data = Data.2,
       aes(x = Index.var,
           y = Trait.var,
           colour=Index)) +
  facet_grid(DAP~Trait, scales = "free_y")+
  geom_smooth(method=lm) +
  geom_point(size = 2) +
  scale_color_grey(start=0.4, end=0.7)+
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

<p align="center">
  <img src="https://raw.githubusercontent.com/filipematias23/images/master/readme/BF_17.jpg">
</p>

<br />

[Menu](#menu)

<div id="p7" />

## Maturity 

> Maturity is another important trait in plant breeding and has great potential for being measured using UAV traits because is related with greenness. For instance, vegetation indices are great tools to identify when the plants are physiological ready to be harvest. In the example below using a simple linear regression, vegetation indices were used to evaluate maturity. It is possible to identify a great moment for collecting data with UAVs (flying day) around 60 DAP, when is already possible to observe high correlation between these traits. The genotype is considered phenotypically matured when 50% of the plants’ peduncle in the plot are yellow as represented on the left site in the image below. 

<p align="center">
  <img src="https://raw.githubusercontent.com/filipematias23/images/master/readme/BF_M.jpg">
</p>

<br />

```r
##########################
### Maturity (MAT_DAY) ###
##########################

### Choosing flights around the MAT_DAY ###
DAP<-c(49,52,67,74,80)
Data<-DataTotal[as.numeric(DataTotal$DAP)%in%DAP,]

### Choosing UAV traits to compare with MAT_DAY ###
Data<-Data[,c("DAP","MAT_DAY","CIG","CIRE")] # Other options: c("CIG","CIRE","Canopy","NGRDI","BGI","GLI","NDVI","NDRE","Height_50","Height_75","Height_90")
Data$DAP<-as.factor(Data$DAP)

### Preparing data to make plots ###
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

### Simple regression visualization ###
ggplot(data = Data.2, 
       aes(x = Index.var,
           y = Trait.var,
           colour=Index)) + 
  facet_grid(DAP~Trait, scales = "free")+
  geom_smooth(method=lm) + 
  geom_point(size = 2) +
  scale_color_grey(start=0.4, end=0.7)+
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

<p align="center">
  <img src="https://raw.githubusercontent.com/filipematias23/images/master/readme/BF_18.jpg">
</p>

<br />

[Menu](#menu)

<div id="p8" />

## Lodging

> Lodging is when grain crops lose their vertical position. This makes them difficult to harvest, and can dramatically reduce yield. This trait is easily measured using the canopy height model (CHM) or estimated plant height (EPH) calculated using the digital surface model (DSM). DSM is one output from the orthomosaicking step when the images have georeferenced information or also called as gridded surface elevation. The idea is to collect images before the plants start emerging (Soil Base) and later in the growth cycle when the plants are already grown (more information in [HERE]( https://github.com/OpenDroneMap/FIELDimageR#P10)). In the example below, images were collected at 73 DAP before a heavy rain event. Another set images were collected in the next day at 74 DAP to investigate the ability of this approach on evaluating lodging. Genotypes with high lodging scores also were among the tallest plants at 73 DAP (Upper set of boxplots in the image below). The same genotypes had a reduction of EPH measured using our UAV platform at 74 DAP after the rain (Center set of boxplots). However, later in the cycle at 80 DAP some of these genotypes had the plasticity to recover and keep upright (Lower set of boxplots). This information can be used by the breeder to infer the ability of recovering of potential tall yield genotypes. 

<p align="center">
  <img src="https://raw.githubusercontent.com/filipematias23/images/master/readme/BF_L.jpg">
</p>

<br />

```r
###############
### Lodging ###
###############

### Flight before: 07/19/2021 (73 DAP) ###
### Flight after_1: 07/20/2021 (74 DAP) ###
### Flight after_2: 07/26/2021 (80 DAP) ###

### Choosing flights around the Lodging event ###
DAP<-c(73,74,80)
Data<-DataTotal[as.numeric(DataTotal$DAP)%in%DAP,]

### Choosing UAV traits to compare with LODG ###
Data<-Data[,c("DAP","LODG","Height_90")] # Other options: c("CIG","CIRE","Canopy","NGRDI","BGI","GLI","NDVI","NDRE","Height_50","Height_75","Height_90") 
Data$DAP<-as.factor(Data$DAP)

### Preparing data to make plots ###
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

### Simple regression visualization ###
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

<p align="center">
  <img src="https://raw.githubusercontent.com/filipematias23/images/master/readme/BF_19.jpg">
</p>

<br />

> In the example below, we compared the same area in field in two different moments of the wheat growth cycle to facilitate lodging visualization using UAV. The first flight was at 52 DAP during the flowering time (all genotypes standing) and the second flight at 74 DAP during the senescence time (some genotypes with lodging). For instance, a great visualization tool is comparing the *Digital Surface Model (DSM)* from both moments to evaluate the variation on average plant height. You can **DOWNLOAD** the example below [**HERE**]().

```r
###################################
### RGB - Lodging Visualization ###
###################################

# unzip("Lodging.zip")

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

<p align="center">
  <img src="https://raw.githubusercontent.com/filipematias23/images/master/readme/BF_20.jpg">
</p>

<br />

```r
###################################
### DSM - Lodging Visualization ###
###################################

dev.off()
par(mfrow=c(1,2))

# 49 DAP #
lodg.52DAP.dsm <- stack("./Lodging/DSM/10_DAP_52_2021_Casselton_YT_06-28_dem.tif")
plot(lodg.52DAP.dsm)

# 74 DAP #
lodg.74DAP.dsm <- stack("./Lodging/DSM/15_DAP_74_2021_Casselton_YT_07-20_dem.tif")
plot(lodg.74DAP.dsm)
```

<p align="center">
  <img src="https://raw.githubusercontent.com/filipematias23/images/master/readme/BF_21.jpg">
</p>

<br />

```r
#############################################
### Drawing Lines - Lodging Visualization ###
#############################################

# Observing EPH profile for 1 rows:
Draw <- fieldDraw(mosaic = lodg.74DAP.dsm,
                     ndraw = 1) # Try changing for 2 to draw 2 lines (Remember to press ESC after finishing the drawing)

# Profile plot:
dev.off()
par(mfrow=c(1,2))
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

<p align="center">
  <img src="https://raw.githubusercontent.com/filipematias23/images/master/readme/BF_22.jpg">
</p>

<br />

[Menu](#menu)

<div id="p9" />

## UAV Data on statistic applications for YIELD selection

> There are many ways to use UAV data to evaluate yield in plant breeding. Below we suggest three simple strategies that can be adapted to different crops. (A) The first suggestion is to use UAV-Traits as covariates in the model. (B) The second is to use AUC-UAV-Traits for indirect selection. (C) The third is to apply indirect selection using UAV-Traits from all flights, however using DAP as cofactor in the model. Please, feel free to contact us and suggest new ways of application. Using NDRE as covariate in the first strategy provided the lower Akaike information criterion (AIC) for both 55 DAP and AUC methods. NDRE-AUC and CANOPY-AUC provided the greater coincidence selection with yield ~ 60% using the indirect selection strategy "B". Strategy "C" had similar performance than strategy "B", however with some inversions on selection coincidence values between NDRE and CIRE.  

```r
##################################
### A. Covariates in the model ###
##################################

dev.off()

### 1) No covariate ###  
# Preparing the data:
Data <- read.csv("EX_DATA.csv",header = T,fileEncoding="UTF-8-BOM")
Data$RANGE<-as.factor(Data$RANGE)
Data$ROW<-as.factor(Data$ROW)
Data$NAME<-as.factor(Data$NAME)
# Mixed model:
# mod<-lmer(YLD~RANGE+ROW+(1|NAME),data = Data)
mod<-lmer(YLD~(1|NAME),data = Data)
# AIC (Comparing models):
(Yield.AIC<-AIC(mod))
# Residuals visualization:
qqPlot(residuals(mod))

# Making a loop:
Trait<-c("NGRDI","NDRE","CIRE","Canopy","Height_50","Height_90")
Data.AIC<-NULL
for(i in 1:length(Trait)){

### 2) Single flight co-variate ###
# Only one flight (e.g., 55 DAP)
DataTotal<-read.csv("DataTotal.csv",header = T)
DAP<-c(55) 
# Preparing the data:
Data<-DataTotal[DataTotal$DAP%in%DAP,]
Data$RANGE<-as.factor(Data$RANGE)
Data$ROW<-as.factor(Data$ROW)
Data$NAME<-as.factor(Data$NAME)
# Mixed model:
# mod<-lmer(eval(parse(text = paste("YLD~",Trait[i],"+RANGE+ROW+(1|NAME)",sep=""))),data = Data)
mod<-lmer(eval(parse(text = paste("YLD~",Trait[i],"+(1|NAME)",sep=""))),data = Data)
# AIC (Comparing models):
Data.AIC<-rbind(Data.AIC,cbind(Trait=Trait[i],AIC=AIC(mod), Model="55DAP")) 
# Residuals visualization:
qqPlot(residuals(mod)) 

### 3) AUC co-variate ### 
# Preparing the data:
DataAUC<-read.csv("DataAUC.csv",header = T)
Data<-DataAUC[as.character(DataAUC$TRAIT)%in%Trait[i],]
Data$RANGE<-as.factor(Data$RANGE)
Data$ROW<-as.factor(Data$ROW)
Data$NAME<-as.factor(Data$NAME)
# Mixed model:
# mod<-lmer(YLD~AUC+RANGE+ROW+(1|NAME),data = Data)
mod<-lmer(YLD~AUC+(1|NAME),data = Data)
# AIC (Comparing models):
Data.AIC<-rbind(Data.AIC,cbind(Trait=Trait[i],AIC=AIC(mod), Model="AUC")) 
# Residuals visualization:
qqPlot(residuals(mod))
}
Data.AIC<-as.data.frame(Data.AIC)
Data.AIC$AIC<-as.numeric(as.character(Data.AIC$AIC))

### AIC visualization ###
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

<p align="center">
  <img src="https://raw.githubusercontent.com/filipematias23/images/master/readme/BF_23a.jpeg"  width="80%" height="80%">
</p>

<br />

```r
################################################################
### B. UAV as a main trait in the model (Indirect Selection) ###
################################################################

### 1) Yield selection based on BLUP (observed data): ###
# Preparing the data:
Data <- read.csv("EX_DATA.csv",header = T,fileEncoding="UTF-8-BOM")
Data$RANGE<-as.factor(Data$RANGE)
Data$ROW<-as.factor(Data$ROW)
Data$NAME<-as.factor(Data$NAME)
# Mixed model:
# mod<-lmer(YLD~RANGE+ROW+(1|NAME),data = Data)
mod<-lmer(YLD~(1|NAME),data = Data)
# BLUPs:
BLUP.Pheno<-as.matrix(ranef(mod)$NAME)
# Ranking:
Sel.Pheno<-rownames(BLUP.Pheno)[order(BLUP.Pheno,decreasing = T)]

### 2) Yield selection based on BLUP (UAV AUC data): ###
DataAUC<-read.csv("DataAUC.csv",header = T)
Trait<-c("NGRDI","NDRE","CIRE","Canopy","Height_50","Height_90")

# Making a loop:
Data.SC<-NULL
for(i in 1:length(Trait)){
# Preparing the data:
  Data<-droplevels(DataAUC[as.character(DataAUC$TRAIT)%in%Trait[i],])
  Data$RANGE<-as.factor(Data$RANGE)
  Data$ROW<-as.factor(Data$ROW)
  Data$NAME<-as.factor(Data$NAME)
# Mixed model:
  # mod<-lmer(AUC~RANGE+ROW+(1|NAME),data = Data)
  mod<-lmer(AUC~(1|NAME),data = Data)
# BLUPs:
  BLUP.AUC<-as.matrix(ranef(mod)$NAME)
# Ranking:
  Sel.AUC<-rownames(BLUP.AUC)[order(BLUP.AUC,decreasing = T)]
# Selection coincidence (25%):
  n.sel<-round(length(Sel.Pheno)*0.25,0)
  Data.SC<-rbind(Data.SC,cbind(Trait=Trait[i],SC=sum(Sel.Pheno[1:n.sel]%in%Sel.AUC[1:n.sel])/n.sel))
}
Data.SC<-as.data.frame(Data.SC)
Data.SC$SC<-as.numeric(as.character(Data.SC$SC))

# Indirect selection coincidence:
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

<p align="center">
  <img src="https://raw.githubusercontent.com/filipematias23/images/master/readme/BF_24a.jpeg"  width="80%" height="80%">
</p>

<br />

```r
###############################################################################
### C. Indirect Selection using UAV traits and DAP as cofactor in the model ###
###############################################################################

### 1) Yield selection based on BLUP (observed data): ###
# Preparing the data:
Data <- read.csv("EX_DATA.csv",header = T,fileEncoding="UTF-8-BOM")
Data$RANGE<-as.factor(Data$RANGE)
Data$ROW<-as.factor(Data$ROW)
Data$NAME<-as.factor(Data$NAME)
# Mixed model:
# mod<-lmer(YLD~RANGE+ROW+(1|NAME),data = Data)
mod<-lmer(YLD~(1|NAME),data = Data)
# BLUPs:
BLUP.Pheno<-as.matrix(ranef(mod)$NAME)
# Ranking:
Sel.Pheno<-rownames(BLUP.Pheno)[order(BLUP.Pheno,decreasing = T)]

### 2) Yield selection based on BLUP (UAV data): ###
DataTotal<-read.csv("DataTotal.csv",header = T)
DAP<-c(31, 41, 49, 55, 67, 73)
DataTotal<-droplevels(DataTotal[DataTotal$DAP%in%DAP,])
# Preparing the data:
DataTotal$DAP<-as.factor(DataTotal$DAP)
DataTotal$RANGE<-as.factor(DataTotal$RANGE)
DataTotal$ROW<-as.factor(DataTotal$ROW)
DataTotal$NAME<-as.factor(DataTotal$NAME)
Trait<-c("NGRDI","NDRE","CIRE","Canopy","Height_50","Height_90")

# Making a loop:
Data.SC<-NULL
for(i in 1:length(Trait)){
  # Mixed model:
  # mod<-lmer(eval(parse(text = paste(Trait[i]," ~ DAP+RANGE+ROW+(1|NAME)",sep=""))),data = DataTotal)
  mod<-lmer(eval(parse(text = paste(Trait[i]," ~ DAP+(1|NAME)",sep=""))),data = DataTotal)
  # BLUPs:
  BLUP.UAV<-as.matrix(ranef(mod)$NAME)
  # Ranking:
  Sel.UAV<-rownames(BLUP.UAV)[order(BLUP.UAV,decreasing = T)]
  # Selection coincidence (25%):
  n.sel<-round(length(Sel.Pheno)*0.25,0)
  Data.SC<-rbind(Data.SC,cbind(Trait=Trait[i],SC=sum(Sel.Pheno[1:n.sel]%in%Sel.UAV[1:n.sel])/n.sel))
}
Data.SC<-as.data.frame(Data.SC)
Data.SC$SC<-as.numeric(as.character(Data.SC$SC))

# Indirect selection coincidence:
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

<p align="center">
  <img src="https://raw.githubusercontent.com/filipematias23/images/master/readme/BF_25a.jpeg"  width="80%" height="80%">
</p>

<br />

[Menu](#menu)

<div id="pC" />

### Forum for questions 

> This discussion group provides an online source of information about the **Bison-Fly Pipeline**. Report a bug and ask a question at: 
* [https://groups.google.com/g/Bison-Fly](https://groups.google.com/g/Bison-Fly) 

<br />

### Licenses

> The R/Bison-Fly package as a whole is distributed under [GPL-2 (GNU General Public License version 2)](https://www.gnu.org/licenses/gpl-2.0.en.html).

<br />

### Citation

> Filipe Inácio Matias, Andrew Green, Jennifer A Lachowiec, David LeBauer, Max Feldman. Bison-Fly: An open-source UAV pipeline for plant breeding data collection. The Plant Phenome Journal, 2022. [https://doi.org/10.1002/ppj2.20048](https://doi.org/10.1002/ppj2.20048)

<br />

### Author

> * [Filipe Inacio Matias](https://github.com/filipematias23)
> * [Andrew Green](https://www.ndsu.edu/agriculture/ag-home/directory/andrew-green)

<br />

### Support and collaboration 

> * [David LeBauer](https://github.com/dlebauer)
> * [Jennifer Lachowiec](https://www.montana.edu/lachowieclab/)
> * [Max Feldman](https://github.com/maxjfeldman)

<br />

### Acknowledgments

> * [Department of Plant Sciences at North Dakota State University](https://www.ndsu.edu/agriculture/academics/academic-units/plant-sciences)
> * [Drone2Phenome]()

<br />

[Menu](#menu)
