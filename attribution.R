

rm(list=ls())

# module load  nixpkgs/16.09
# module load gcc/7.3.0
# module load netcdf/4.6.1
# module load r/3.6.1
# module load gdal/3.0.1
# module load proj/6.3.0
# module load udunits/2.2.26



library(ncdf4)

precip_file <- nc_open("/scratch/xtan/LENS_monthly/ERA5/total_precipitation.nc")
runoff_file <- nc_open("/scratch/xtan/LENS_monthly/ERA5/Runoff.nc")
aet_file <- nc_open("/scratch/xtan/LENS_monthly/ERA5/AET.nc")
pet_file <- nc_open("/scratch/xtan/LENS_monthly/ERA5/PET.nc")
temp_file <- nc_open("/scratch/xtan/LENS_monthly/ERA5/2m_temperature.nc")
dew_file <- nc_open("/scratch/xtan/LENS_monthly/ERA5/2m_dewpoint_temperature.nc")
uwind_file <- nc_open("/scratch/xtan/LENS_monthly/ERA5/10m_U_wind.nc")
vwind_file <- nc_open("/scratch/xtan/LENS_monthly/ERA5/10m_V_wind.nc")


time <- ncvar_get(precip_file, varid = "time")
lon <- ncvar_get(precip_file, varid = "longitude")
lat <- ncvar_get(precip_file, varid = "latitude")

library(ncdf4)
q_file = nc_open("Q.nc")

longitude = ncvar_get(q_file, varid = "longitude")[241:721]
latitude = ncvar_get(q_file, varid = "latitude")[161:361]
time = ncvar_get(q_file, varid = "time")
expver = ncvar_get(q_file, varid = "expver")
x = ncdim_def("lon","degreesE",longitude)
y = ncdim_def("lat","degreesN",latitude)
exp = ncdim_def("expver","unit",expver)
t = ncdim_def("time","days since 1900-01-01",time)

uq1 = ncvar_get(q_file, varid = "p71.162", start = c(241,161,1,1), count = c(481,201,-1,-1) )
uq <- ncvar_def( "uq", "units", list(x, y , exp, t) )
nc <- nc_create( "Asia_moisture_uq.nc", list(uq) )
ncvar_put( nc, uq, uq1 )

uq1 = ncvar_get(q_file, varid = "p72.162", start = c(241,161,1,1), count = c(481,201,-1,-1))
vq <- ncvar_def( "vq", "units", list(x, y , exp, t) )
nc <- nc_create( "Asia_moisture_vq.nc", list(vq) )
ncvar_put( nc, vq, uq1 )

uq1 = ncvar_get(q_file, varid = "vimd", start = c(241,161,1,1), count = c(481,201,-1,-1))
dq <- ncvar_def( "dq", "units", list(x, y , exp, t) )
nc <- nc_create( "Asia_moisture_dq.nc", list(dq) )
ncvar_put( nc, dq, uq1 )

# south China lat: 20-30N (601:700) lon:100-123E (1001:1230)


precip <- ncvar_get(precip_file, start=c(1001, 601, 1), count=c(230, 100, -1), varid = "tp")
precip <- precip * 1000

AET <- ncvar_get(aet_file, start=c(1001, 601, 1), count=c(230, 100, -1), varid = "e")
AET <- AET * 1000 * (-1)

PET <- ncvar_get(pet_file, start=c(1001, 601, 1), count=c(230, 100, -1), varid = "pev")
PET <- PET * 1000 * (-1) 

runoff <- ncvar_get(runoff_file, start=c(1001, 601, 1), count=c(230, 100, -1), varid = "ro")
runoff <- runoff * 1000

temp <- ncvar_get(temp_file, start=c(1001, 601, 1), count=c(230, 100, -1), varid = "t2m")

dew_temp <- ncvar_get(dew_file, start=c(1001, 601, 1), count=c(230, 100, -1), varid = "d2m")

uwind <- ncvar_get(uwind_file, start=c(1001, 601, 1), count=c(230, 100, -1), varid = "u10")

vwind <- ncvar_get(vwind_file, start=c(1001, 601, 1), count=c(230, 100, -1), varid = "v10")

save(lon, lat, precip, AET, PET, runoff, temp, dew_temp, uwind, vwind, file = "south_china_ERA5.RData")

precip_area_mean <- apply(X=precip, MARGIN =c(3), FUN = mean, na.rm=T )


library(extRemes)
library(RcppRoll)
library(raster)
library(abind)
library(weathermetrics)
library(SPEI)

# dat <- qevd(seq(1e-8,1-1e-8,length.out=481)[1:480], 1, 0.5, 0.8)
load("precip_ERA_5.RData")
load("south_china_ERA5.RData")
lat <- lat[601:700]
lon <- lon[1001:1230]


dat1 <- precip[,,480]

# grid <- expand.grid(y=lon,x=lat)
grid <- expand.grid(y=seq(100,122.9,0.1),x=seq(30,20.1,-0.1))
dat <- data.frame(x=grid$y,y=grid$x,z=as.vector(dat1))
crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
dat <- raster::rasterFromXYZ(xyz=dat, crs = crs)
plot(dat)

#####################################################################

load("south_china_ERA5.RData")
duration <- c(3,4,5,6,7,8)
year <- 1981:2020
month <- 1:12
year_month <- expand.grid(x=month, y=year)
year_month <- data.frame (year = year_month$y, month= year_month$x)
grid <- expand.grid(y=seq(100,122.9,0.1),x=seq(30,20.1,-0.1))
for(i_d in 1:length(duration)){
  fill <- apply(X=precip, MARGIN = c(1,2), FUN = mean, na.rm=T )  
  fill <- array(fill, dim=c(1, 230, 100) )*duration[i_d]
  
  dat_roll <- apply(X=precip, MARGIN = c(1,2), FUN = roll_sum, n=duration[i_d] , na.rm=T ) 
  for(id in 1:(i_d+1)){
    dat_roll <- abind(fill, dat_roll, along = c(1))
  }
  
  
  fill <- apply(X=PET, MARGIN = c(1,2), FUN = mean, na.rm=T )  
  fill <- array(fill, dim=c(1, 230, 100) )*duration[i_d]
  
  pet_roll <- apply(X=PET, MARGIN = c(1,2), FUN = roll_sum, n=duration[i_d] , na.rm=T ) 
  for(id in 1:(i_d+1)){
    pet_roll <- abind(fill, pet_roll, along = c(1))
  }
  
  
  humidity <- dewpoint.to.humidity(dp= dew_temp-273.15, t=temp-273.15, temperature.metric = "celsius")
  fill <- apply(X=humidity, MARGIN = c(1,2), FUN = mean, na.rm=T )  
  fill <- array(fill, dim=c(1, 230, 100) )
  
  humidity_roll <- apply(X=humidity, MARGIN = c(1,2), FUN = roll_mean, n=duration[i_d] , na.rm=T ) 
  for(id in 1:(i_d+1)){
    humidity_roll <- abind(fill,  humidity_roll, along = c(1))
  }
  
  ######## calculation of SPEI
  
  
  for(i in 1:dim(precip)[1]){
    results <- array(NA, c(dim(precip)[2:3],4) )
    for(j in 1:dim(precip)[2]){

      ############# anormalies of variables
      # dat <- cbind(year_month,
      #              precip= anomalies ( cbind(year_month, dat=dat_roll[,i,j])),
      #              pet= anomalies ( cbind(year_month, dat=pet_roll[,i,j])),
      #              P_E= anomalies ( cbind(year_month, dat=dat_roll[,i,j]- pet_roll[,i,j] )) ,
      #              humidity = anomalies ( cbind(year_month, dat=humidity_roll[,i,j])) )

     # plot(dat$precip,type="l")
     # plot(dat$P_E,type="l")

      dat1 <- cbind(year_month, precip= precip[i,j,], pet= PET[i,j,],  humidity=humidity[i,j,] )
      if(length(dat1$precip [is.na(dat1$precip)] ) < 10 ){
        spi <- spi(dat1, scale = duration[i_d])$fitted[,3:5]
        P_E <- cbind(year_month, P_E= precip[i,j,]-PET[i,j,] )
        spei <- spei(data= P_E$P_E, scale = duration[i_d])$fitted
        spi <- cbind(spei, spi)
      }
      else{
        spi <- array(NA, dim = c(480,4) )
      }

      results[j,,] <- spi
    }
    save(results, file = paste0("spi_4_variables_",i_d,"_",i,".RData") )
  }
  
}

setwd("D:\\South_China\\spi")
load("south_china_ERA5.RData")
duration <- c(3,4,5,6,7,8)
grid <- expand.grid(y=seq(100,122.9,0.1),x=seq(30,20.1,-0.1))
for(i_d in 2:length(duration)){
  
  dat <- array(NA,dim = c(230,100,480,4))
  for(i in 1:dim(precip)[1]){
    load(file = paste0("spi_4_variables_",i_d,"_",i,".RData") )
    dat[i,,,] <- results
  }
  
  grid <- expand.grid(y=seq(100,122.9,0.1),x=seq(30,20.1,-0.1))
  dat1 <- data.frame(x=grid$y,y=grid$x,z=as.vector(dat[,,310,1])) ######### make the raster data for the 100 month of spei
  crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  dat1 <- raster::rasterFromXYZ(xyz=dat1, crs = crs)
  plot(dat1)
  
  library(trend)
  for(i_var in 1:4){
    spei <- dat[,,,i_var]
    trend_spei <- trend_sig_spei <- array(NA,dim=dim(spei)[1:2])
    for(i_lon in 1:nrow(spei)){
      for(i_lat in 1:ncol(spei)){
        time_series <- spei[i_lon, i_lat,][!is.na(spei[i_lon, i_lat,])]
        if(length(time_series)>0){
          test <- sens.slope(x=time_series)
          trend_spei[i_lon, i_lat] <- test$estimates
          trend_sig_spei[i_lon, i_lat] <- test$p.value
        }
        
      }
    }
    
    dat1 <- data.frame(x=grid$y,y=grid$x,z=as.vector(trend_spei)) 
    ######### make the raster data for the 100 month of spei
    crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
    dat1 <- raster::rasterFromXYZ(xyz=dat1, crs = crs)
    writeRaster(dat1,filename = paste0("trend_spi_4_variables_",i_d,"_",i_var,".nc"), overwrite=TRUE )
    
    dat2 <- data.frame(x=grid$y,y=grid$x,z=as.vector(trend_sig_spei)) 
    dat2 <- raster::rasterFromXYZ(xyz=dat2, crs = crs)
    # plot(dat1)
    writeRaster(dat2,filename = paste0("trend_sig_spi_4_variables_",i_d,"_",i_var,".nc"), overwrite=TRUE )
  }

}

library(raster)
############ variables with P-E, precip, pet,  humidity
var = c("P-E","Precipitation", "Potential Evapotranspiration","Humidity")
for(i_d in 1:6){
  for(i_var in 1:4){
    trend = raster(paste0("trend_spi_4_variables_",6,"_",i_var,".nc"))
    trend = as.data.frame(trend, xy=TRUE)
    sig = raster(paste0("trend_sig_spi_4_variables_",6,"_",i_var,".nc"))
    sig = as.data.frame(sig, xy=TRUE)
    if(i_var %in% 3){
      p = ggplot()+
        geom_raster(data = trend, mapping=aes(x=x, y=y, fill=z))+
        geom_contour(data = sig, mapping=aes(x=x, y=y, z=z),breaks=0.05) +
        labs(x = 'Longitude', y = 'Latitude', title = 'SPI_12', size = 10 ) +
        scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = .0, 
                             breaks = c(-0.004,-0.003,-0.002,-0.001,0,0.001,0.002,0.003,0.004),
                             guide = guide_colourbar(direction = "horizontal", barwidth = 20,
                                                     title.position = "bottom" ) ) + 
        theme(panel.background = element_rect(fill="white",color="black",size=0.3),
              panel.grid = element_blank(),legend.position='bottom',
              axis.title=element_text(size=10),axis.text=element_text(size=8),
              axis.ticks=element_line(size=0.1))
    }
    else{
      p = ggplot()+
        geom_raster(data = trend, mapping=aes(x=x, y=y, fill=z))+
        geom_contour(data = sig, mapping=aes(x=x, y=y, z=z),breaks=0.05) +
        labs(x = 'Longitude', y = 'Latitude', title = 'SPI_12', size = 10 ) +
        scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = .0,
                             breaks = c(-0.004,-0.003,-0.002,-0.001,0,0.001,0.002,0.003,0.004),
                             guide = guide_colourbar(direction = "horizontal", barwidth = 20,
                                                     title.position = "bottom" ) ) + 
        theme(panel.background = element_rect(fill="white",color="black",size=0.3),
              panel.grid = element_blank(),legend.position='bottom',
              axis.title=element_text(size=10),axis.text=element_text(size=8),
              axis.ticks=element_line(size=0.1))
    }
    ggsave(paste0("trend_4_variables_",i_d,"_",var[i_var],".pdf"),p,width = 5, height = 5)
  }
}

month = rep(1:12, 40)

dd = dim(precip)
precip_monthly_anomalies = array(NA,dim = dim(precip))
for(i in 1:dd[1]){
  for(j in 1:dd[2]){
    for(k in 1:dd[3]){
      precip_monthly_anomalies[i,j,k] = precip[i,j,k] - 
        mean(precip[i,j, seq(k%%12,dd[3],12 ) ], na.rm = T)
    }
  }
}
precip_monthly_anomalies_mean = apply(precip_monthly_anomalies, c(3), FUN = mean, na.rm=T)
write.csv(precip_monthly_anomalies_mean, file = "precip_monthly_anomalies_mean_ERA5.csv" )

precip_spring = precip_summer = precip_fall = precip_winter = array(NA, dim = c(dim(precip)[1:2],40))
for(i in 1:40){
  precip_spring[,,i] <- apply(precip[,,((i-1)*12+3) : ((i-1)*12+5) ], c(1,2), FUN = sum, na.rm=T)
  precip_summer[,,i] <- apply(precip[,,((i-1)*12+6) : ((i-1)*12+8) ], c(1,2), FUN = sum, na.rm=T)
  precip_fall[,,i] <- apply(precip[,,((i-1)*12+9) : ((i-1)*12+11) ], c(1,2), FUN = sum, na.rm=T)
}
precip_winter[,,1] <- apply(precip[,,c(1,2,480) ], c(1,2), FUN = sum, na.rm=T)
for(i in 1:39){
  precip_winter[,,i+1] <- apply(precip[,,c((i-1)*12+12,i*12+1,i*12+2) ], c(1,2), FUN = sum, na.rm=T)
}

precip_spring_climtology <- apply(X=precip_spring, MARGIN=c(1,2), FUN = mean, na.rm=T )
precip_summer_climtology <- apply(X=precip_summer, MARGIN=c(1,2), FUN = mean, na.rm=T )
precip_fall_climtology <- apply(X=precip_fall, MARGIN=c(1,2), FUN = mean, na.rm=T )
precip_winter_climtology <- apply(X=precip_winter, MARGIN=c(1,2), FUN = mean, na.rm=T )

precip_spring_anomalies <- precip_summer_anomalies <- precip_fall_anomalies <- precip_winter_anomalies <- precip_spring

for(i in 1:40){
  precip_spring_anomalies[,,i] <-  (precip_spring[,,i] - precip_spring_climtology)/precip_spring_climtology*100
  precip_summer_anomalies[,,i] <-  (precip_summer[,,i] - precip_summer_climtology)/precip_summer_climtology*100
  precip_fall_anomalies[,,i] <-  (precip_fall[,,i] - precip_fall_climtology)/precip_fall_climtology*100
  precip_winter_anomalies[,,i] <-  (precip_winter[,,i] - precip_winter_climtology)/precip_winter_climtology*100
}

plot_raster (dat=precip_summer_anomalies[,,40],  title = "Summer Precipitation ")
plot_raster (dat=precip_fall_anomalies[,,40], title = "Fall Precipitation")
plot_raster (dat=precip_spring_anomalies[,,40], title = "SpringPrecipitation ")

for(i in 1:40){
  precip_spring_anomalies[,,i] <-  (precip_spring[,,i] - precip_spring_climtology)*30
  precip_summer_anomalies[,,i] <-  (precip_summer[,,i] - precip_summer_climtology)*30
  precip_fall_anomalies[,,i] <-  (precip_fall[,,i] - precip_fall_climtology)*30
  precip_winter_anomalies[,,i] <-  (precip_winter[,,i] - precip_winter_climtology)*30
}

plot_raster (dat=precip_summer_anomalies[,,40],  title = "Summer Precipitation ")
plot_raster (dat=precip_fall_anomalies[,,40], title = "Fall Precipitation")
plot_raster (dat=precip_spring_anomalies[,,40], title = "SpringPrecipitation ")

region_mean_precip_summer <- apply(precip_summer[101:201, 21:100,], MARGIN = c(3), FUN = mean, na.rm=T )*30
region_mean_precip_spring <- apply(precip_spring[101:201, 21:100,], MARGIN = c(3), FUN = mean, na.rm=T )*30
region_mean_precip_fall <- apply(precip_fall[101:201, 21:100,], MARGIN = c(3), FUN = mean, na.rm=T )*30
region_mean_precip_winter <- apply(precip_winter[101:201, 21:100,], MARGIN = c(3), FUN = mean, na.rm=T )*30
save(year, region_mean_precip_summer,region_mean_precip_spring ,region_mean_precip_fall,region_mean_precip_winter,
     file = "ERA_South_China_precipitation.RData")

region_anomalies_precip_summer <- apply(precip_summer_anomalies[101:201, 21:100,], MARGIN = c(3), FUN = mean, na.rm=T )
region_anomalies_precip_spring <- apply(precip_spring_anomalies[101:201, 21:100,], MARGIN = c(3), FUN = mean, na.rm=T )
region_anomalies_precip_fall <- apply(precip_fall_anomalies[101:201, 21:100,], MARGIN = c(3), FUN = mean, na.rm=T )
region_anomalies_precip_winter <- apply(precip_winter_anomalies[101:201, 21:100,], MARGIN = c(3), FUN = mean, na.rm=T )
save(year, region_anomalies_precip_summer,region_anomalies_precip_spring ,
     region_anomalies_precip_fall,region_anomalies_precip_winter,
     file = "ERA_South_China_precipitation_anomalies.RData")

load(file = "ERA_South_China_precipitation_anomalies.RData")
annual_anomalies = region_anomalies_precip_spring  + region_anomalies_precip_summer +
                     region_anomalies_precip_fall + region_anomalies_precip_winter

############### frequency analysis ###########################

i = 1 
j  = 1

dat <- cbind(year_month, precip= dat_roll[,i,j], pet= pet_roll[,i,j], P_E= dat_roll[,i,j]- pet_roll[,i,j], humidity=humidity_roll[,i,j] )
spi <- spi(dat$precip, scale = 1)
plot(spi)
spei <- spei(dat$P_E, scale=1)
plot(spei)
spi_humidity <- spi(dat$humidity, scale=1)
plot(spi_humidity)

############## guangdong province lon 109:118 lat 20:24 ######################
GD <- c(91,181,61,100)
load("south_china_ERA5.RData")
precip <- precip[GD[1]:GD[2],GD[3]:GD[4],]
precip <- apply(X= precip, MARGIN =  c(3), FUN = mean, na.rm=T)
PET <- PET[GD[1]:GD[2],GD[3]:GD[4],]
PET <- apply(X= PET, MARGIN =  c(3), FUN = mean, na.rm=T)
humidity <- dewpoint.to.humidity(dp= dew_temp-273.15, t=temp-273.15, temperature.metric = "celsius")
humidity <- humidity[GD[1]:GD[2],GD[3]:GD[4],]
humidity <- apply(X= humidity, MARGIN =  c(3), FUN = mean, na.rm=T)


dat <- cbind(year_month, precip= precip, pet= PET, P_E= precip-PET, humidity=humidity )
spi <- spi(ts(dat$precip, start = c(1981,1), end = c(2020, 12), frequency = 12 ), scale = 6)
plot(spi)
spei <- spei(ts(dat$P_E, start = c(1981,1), end = c(2020, 12), frequency = 12 ), scale=6)
plot(spei)
spi_humidity <- spi(ts(dat$humidity, start = c(1981,1), end = c(2020, 12), frequency = 12 ), scale=6)
plot(spi_humidity)

