#Lehmann et al. 2018
#Complex responses of global insect pests to climate change
#================================================================================
#Ambient Temperatures extracted from the HADLEY model outputs 
#model=HadGEM2-CC; 
#time=Monthly
#Experiments: historical and rcp85
#sourced at:https://esgf-node.llnl.gov/projects/esgf-llnl/
#================================================================================

rm(list=ls())
#packages required:
library(maptools)
library(raster)   
library(rgeos)
library(ncdf4)
library(rgdal)

#--------------------------------------------------------------------------------------
#Prepare climate files for extraction
#--------------------------------------------------------------------------------------

#list names of downloaded datafiles
filename<-c("tas_Amon_HadGEM2-CC_historical_r3i1p1_195912-198411",
            "tas_Amon_HadGEM2-CC_rcp85_r3i1p1_200512-203011",
            "tas_Amon_HadGEM2-CC_rcp85_r3i1p1_205512-208011",
            "tas_Amon_HadGEM2-CC_rcp85_r3i1p1_205512-208011")

#list relevant information for the following loop
names<-c("historical", "current", "near_future", "future")
start_list<-c("X1960.01.15", "X2006.01.15", "X2056.01.15", "X2070.01.15")
end_list<-c("X1969.12.16", "X2015.12.16",  "X2065.12.16","X2079.12.16")

#loop through files, extract relevant years, apply global mask, and write cleaned
#climate data (as .nc) to file for later use
for (i in (1:length(filename))){
  pathname<-paste( "./path/",filename[i],".nc",sep="")
  nc <- nc_open(pathname) #call in data
  b<- brick(pathname, varname="tas") #as a rasterbrick
  # b #check ok
  start <- which((names(b) == start_list[i])==TRUE) #define start layer
  end <- which((names(b) == end_list[i])==TRUE)     #define end layer
  time<-seq(start,end,1)    #list layer  levels according to start and end
  b1<-subset(b, time)  #subset out relevant layers
  #plot(b1) #check ok
  #b1[1,1,] #check ok 
  #assign projection
  projectExtent(b1, "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0") 
  
  myshp <-rgdal::readOGR("~/ne_110m_admin_0_countries.shp") #call in global shape file
  e <- extent(myshp)  #define the extent according to the global map
  myraster.crop <- crop(b1, e, snap="out") #crop the layers to include land only with global extent
  crop <- setValues(myraster.crop, NA) #all sea areas get NA
  myshp.r <- rasterize(myshp, crop) #raster the new shape file
  xx <- raster::mask(x=myraster.crop, mask=myshp.r)#mask the two rasters
  
  b2<-rotate(b1) #rotate the map to accout for longitude definition of CMIP files (0-360)
  myraster.crop <- crop(b2, e, snap="out") #repeat the above on the western hemisphere
  crop <- setValues(myraster.crop, NA) #as above
  myshp.r <- rasterize(myshp, crop) #as above
  xx2 <- raster::mask(x=myraster.crop, mask=myshp.r) #as above
  
  r3<-merge(xx, xx2) #combine the two hemispheres, so now -180 to 180
  names(r3)<-names(xx2) #assign the correct name to the new raster
  #plot(r3) #check ok
  
  r3 <- r3-273 #convert from Kelvin to Celsius
  #Assign projection 
  projectExtent(r3, "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0") 
  #write raster output as new .nc file, only the landareas required
  writeRaster(r3, filename=paste( "./path/", names[i],".nc",sep=""), overwrite=TRUE)
}


#--------------------------------------------------------------------------------------
#Extract Tamb from the species sites:
#--------------------------------------------------------------------------------------

#Prepare the co-ordinates for the species locations
my_coords<-as.data.frame(read.csv(file = "./path/my_coordinates.csv", sep = ",", header=TRUE))

#set up new filenames for each timeperiod
filename <- c("historical", "current", "near_future", "future")

#loop through the four prepped climate files and extract data from each location
for (i in (1:4)){
  pathname<-paste("./path/",filename[i],".nc",sep="")
  nc <- nc_open(pathname) #call in the cleaned climate .nc file
  b<- brick(pathname, varname="variable") 
  plot(b$X1) 
  points(my_coords$lat~my_coords$lon, pch = 19, cex=0.5 ) #plot the points to check ok
  abline(h=0, col="red", lty=2) #plot equator to check ok

  points.sp <- na.omit(my_coords) #omit missing locations from lat/lon data
  coordinates(points.sp) <- ~ lon + lat #convert points to coordinates
  climate.points <- raster::extract(b, points.sp) #extract coordinates from climate raster
  print(climate.points) #check ok

  #transpose and clean output file
  output <- cbind(na.omit(my_coords), climate.points) #remove any missing values (not within global extent)
  output <- as.data.frame(t(climate.points)) #transpose data
  colnames(output) <- points.sp$Common.name #assign species name to headers
  output$month <- c(rep(1:12,10)) #include month column
  output$year <- c(rep(1:10,each=12)) #include year column
  head(output) #check data structure ok

  #select only the averages required months of growing season: 
  #if lat higher than 45, average over months 5-9, otherwise average
  # over all months of the year
  ave_temp = list() #generate an empty list into which average temperatures are saved
  for (j in (1:(ncol(output)-2))){
    if(points.sp$lat[j]>45){
     ave_temp[j] = mean(subset(output, (month>=5 & month<=9))[,j])  
   }else{
     ave_temp[j] = mean(output[,j])
   }
  }
  #combine temperature with the lat/lons and assign timestamp before returning to top of loop
  temperature <- unlist(ave_temp)
  final.output <-cbind(na.omit(my_coords, temperature))
  assign(filename[i], final.output) #assign name to final.output and continue to top of loop
}

#combine all decades into one dataframe and write to file
fin<-as.data.frame(cbind(historical, current$temperature, near_future$temperature, future$temperature))
colnames(fin)<-c("lat","lon","common_name","historical","current" ,"near_future","future")
head(fin)
write.csv(fin,"./path/output_filname.csv") 

#----------------------------------------------------------------------------------------------
#END_SCRIPT
