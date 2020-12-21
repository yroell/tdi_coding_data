setwd("/Users/yroell/Desktop/tdi/section1/")

library(XML)
library(rgdal)
library(gdalUtils)
library(raster)
library(sf)
library(dplyr)
library(plyr)
library(ranger)
library(randomForest)
library(caret)
library(elevatr)

# load in colorado county shapefile
colo = st_read("/Users/yroell/Desktop/tdi/section1/colorado/census_CENDEMG_co_3874117_01/census/census2015_demo_a_co.gdb", layer = "county_census2010dp1_a_co_2015")
colo_proj = st_transform(colo, crs = 26912)

coord_sys = st_crs(colo_proj)

bbox = st_bbox(colo_proj)

ulx = bbox$xmin
uly = bbox$ymax
lrx= bbox$xmax
lry = bbox$ymin
(bb <- c(ulx, uly, lrx, lry))

# download climate data
climate = getData("worldclim", var = "prec", res = 2.5)
precip = crop(climate, colo)
precip_proj = projectRaster(precip, crs = coord_sys[2]$proj4string)

# load in elevation data and create indices
dem = raster("/Users/yroell/Desktop/tdi/section1/terrain/output_srtm.tif")
ele = crop(dem, colo)
ele_proj = projectRaster(ele, precip_proj[[1]])
slope_ele = terrain(ele_proj, opt = "slope")
tpi_ele = terrain(ele_proj, opt = "TPI")
tri_ele = terrain(ele_proj, opt = "TRI")
rough_ele = terrain(ele_proj, opt = "roughness")

terrain_stack = stack(c(ele_proj, slope_ele, tpi_ele, tri_ele, rough_ele))

# download soil data
#url = "https://files.isric.org/soilgrids/latest/data/"

sg_url="/vsicurl/https://files.isric.org/soilgrids/latest/data/"

# get bulk density raster
voi_url = "bdod/bdod_0-5cm_mean.vrt"
lfile = "./soil/bdod_0-5.tif"
gdal_translate(paste0(sg_url, voi_url), lfile,
               tr = c(250,250),
               projwin = bb,
               projwin_srs = coord_sys[2]$proj4string,
               verbose= TRUE)
bulk_density = raster(lfile)/100

# get clay percentage raster
voi_url = "clay/clay_0-5cm_mean.vrt"
lfile = "./soil/clay_0-5.tif"
gdal_translate(paste0(sg_url, voi_url), lfile,
               tr = c(250,250),
               projwin = bb,
               projwin_srs = coord_sys[2]$proj4string,
               verbose= TRUE)
clay_percent = raster(lfile)/100

# get sand percentage raster
voi_url = "sand/sand_0-5cm_mean.vrt"
lfile = "./soil/sand_0-5.tif"
gdal_translate(paste0(sg_url, voi_url), lfile,
               tr = c(250,250),
               projwin = bb,
               projwin_srs = coord_sys[2]$proj4string,
               verbose= TRUE)
sand_percent = raster(lfile)/100

# get clay percentage raster
voi_url = "soc/soc_0-5cm_mean.vrt"
lfile = "./soil/soc_0-5.tif"
gdal_translate(paste0(sg_url, voi_url), lfile,
               tr = c(250,250),
               projwin = bb,
               projwin_srs = coord_sys[2]$proj4string,
               verbose= TRUE)
soc_percent = raster(lfile)/100

soil_stack = stack(c(bulk_density, clay_percent, sand_percent, soc_percent))
soil_stack_proj = projectRaster(soil_stack, precip_proj[[1]])

# stack all data together
covariate = stack(soil_stack_proj, terrain_stack, precip_proj)

# load in landslide data
landslide = st_read("/Users/yroell/Desktop/tdi/section1/landslide/Landslide_Locations%3A_Global_Landslide_Catalog__GLC_.shp")
landslide_colo = subset(landslide, admin_divi == "Colorado")
landslide_proj = st_transform(landslide_colo, crs = 26912)

# extract covariate values to target variable
sp_df = extract(covariate, landslide_proj, sp = TRUE)
df = as.data.frame(sp_df[, c(13, 37:57)])
df = df[complete.cases(df), ]

# create random forest model
df$landslid_2 = droplevels(df$landslid_2)
df$landslid_2 = revalue(df$landslid_2, c("small" = 1, "medium" = 2, "large" = 3, "very_large" = 3))
df$landslid_2 = factor(df$landslid_2, levels = c(1, 2, 3))
target = df$landslid_2
covar = df[, 2:22]

rf = randomForest(x = covar, y = target)
varImpPlot(rf, main = "Variable Importance")

# create landslide map
predict_landslide = predict(covariate, rf)
plot(predict_landslide)

# determine landslide per county
county_sp = extract(predict_landslide, colo_proj, fun = mean, na.rm = TRUE, sp = TRUE)
county_sp$risk = county_sp$layer >= mean(county_sp$layer)
county_sp$risk = as.factor(county_sp$risk)
county_sp$risk = revalue(county_sp$risk, c("FALSE" = "normal", "TRUE" = "high_risk"))
spplot(county_sp, "layer")
spplot(county_sp, "DP0010001")
spplot(county_sp, "risk")
county = as.data.frame(county_sp)
high_risk = subset(county, risk == "high_risk")
sum(high_risk$DP0010001) / sum(county$DP0010001)

# save output
writeRaster(covariate, filename = "./output/covariates.tif")
writeRaster(predict_landslide, filename = "./output/landslide_prediction.tif")
writeOGR(county_sp, "./output", "county_sp", driver = "ESRI Shapefile")
write.csv(county, "county.csv")

#setwd("/Users/yroell/Desktop/esri/")
#r = covariate
#nlayers(r)
#for(i in 1:nlayers(r)){
#  band<-r[[i]]
  #save raster in a separate file
#  writeRaster(band,paste(names(covariate)[i],'.tif', sep=''))
#}
