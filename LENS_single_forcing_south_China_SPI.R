

# 
# module load  nixpkgs/16.09
# module load gcc/7.3.0
# module load netcdf/4.6.1
# module load r/3.6.1
# module load gdal/3.0.1
# module load proj/6.3.0

R

rm(list=ls())

setwd("F:/data/LENS_single_forcing")

forcings = c("AER","BMB","GHG","hist","LULC")

models = list( c("001","002","003","004","005","007","008","009","010","011","012","013","014","015","016","017","018","019","020"),
               c("001","002","003","004","005","006","007","008","009","010","011","012","013","014","015"),
               c("001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020"), 
               c("002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020",
                 "021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","101","102","103","104","105"),
               c("001","002","003","004","005")      )


######################################################################################################
######################################  SPI frequency analyses #######################################
######################################################################################################


# module load udunits/2.2.26

library(ncdf4)
library(abind)
library(SPEI)
library(RcppRoll)
#library(standaRdized)
spi_calc = function(precip, ref.data){
  d = dim(precip)
  spi = precip * NA
  for(i_model in 1:d[1]){
    for(i_lon in 1:d[2]){
      spi[i_model, i_lon, ,] = t (spei(data=t(precip[i_model, i_lon,, ]), scale=12, kernel=list(type='rectangular',shift=0),
                                       distribution='Gamma', fit='ub-pwm', na.rm=T, 
                                       ref.data=t(ref.data[ i_lon, , ]))$fitted)
    }
  }
  return(spi)
}

for(i_forcing in 1){  ####### AER forcing 
  precip = array(NA, dim = c(19,9,10,1932))
  
  # spi = precip = array(NA, dim = c(19,9,10,1932))
  for(i_model in 1:length(models[[i_forcing]])){
    
    years = "192001-200512"
    
    file = nc_open( paste0( "/scratch/xtan/LENS_single_forcing/PRECL/b.e11.B20TRLENS_RCP85.f09_g16.xaer.",
                            models[[i_forcing]][i_model],".cam.h0.PRECL.", years,".nc" ))
    
    #  spi[lon,lat,time]
    ## stduy region (part of south China) lat 20:28 (117:126) , lon, 110-120 (88:96)
    # precip_south_china = dat(lat|197:212, lon|131:144,  time|:, height|0)    ; extract the study region data
    
    dat = ncvar_get( file, varid="PRECL", start = c(89,96,1), count = c(9,10,-1) )
    
    file  = nc_open( paste0( "/scratch/xtan/LENS_single_forcing/PRECC/b.e11.B20TRLENS_RCP85.f09_g16.xaer.",
                             models[[i_forcing]][i_model], ".cam.h0.PRECC.", years, ".nc"))
    
    dat1 = ncvar_get( file, varid="PRECC", start = c(89,96,1), count = c(9,10,-1) )
    
    precip1 = dat + dat1
    
    years = "200601-208012"
    
    file = nc_open( paste0( "/scratch/xtan/LENS_single_forcing/PRECL/b.e11.B20TRLENS_RCP85.f09_g16.xaer.",
                            models[[i_forcing]][i_model],".cam.h0.PRECL.", years,".nc" ))
    
    #  spi[lon,lat,time]
    ## stduy region (part of south China) lat 20:28 (117:126) , lon, 110-120 (88:96)
    # precip_south_china = dat(lat|197:212, lon|131:144,  time|:, height|0)    ; extract the study region data
    
    dat = ncvar_get( file, varid="PRECL", start = c(89,96,1), count = c(9,10,-1) )
    
    file  = nc_open( paste0( "/scratch/xtan/LENS_single_forcing/PRECC/b.e11.B20TRLENS_RCP85.f09_g16.xaer.",
                             models[[i_forcing]][i_model], ".cam.h0.PRECC.", years, ".nc"))
    
    dat1 = ncvar_get( file, varid="PRECC", start = c(89,96,1), count = c(9,10,-1) )
    
    precip2 = dat + dat1
    
    precip[i_model,,,] = abind(precip1, precip2, along = 3)
    
    # file = nc_open(filename = paste0("/scratch/xtan/LENS_single_forcing/spi/aer/AER_spi_12_",models[[i_forcing]][i_model],".nc"))
    # time <- ncvar_get(file, varid = "time")
    # lon <- ncvar_get(file, varid = "lon")
    # lat <- ncvar_get(file, varid = "lat")
    
    # spi[i_model,,,] = ncvar_get(file, varid="spi", start = c(89,96,1), count = c(9,10,-1) )
  }
  
  # save(precip, spi, file = paste0("precip_spi_aer.RData"))
  save(precip, file = paste0("precip_spi_aer.RData"))
}

for(i_forcing in 2){  ####### BMB forcing 
  spi = precip = array(NA, dim = c(length(models[[i_forcing]]),9,10,1320))
  years = "192001-202912"
  
  for(i_model in 1:length(models[[i_forcing]])){
    file = nc_open( paste0( "/scratch/xtan/LENS_single_forcing/PRECL/b.e11.B20TRLENS_RCP85.f09_g16.xbmb.",
                            models[[i_forcing]][i_model],".cam.h0.PRECL.", years,".nc" ))
    
    #  spi[lon,lat,time]
    ## stduy region (part of south China) lat 20:28 (117:126) , lon, 110-120 (88:96)
    # precip_south_china = dat(lat|197:212, lon|131:144,  time|:, height|0)    ; extract the study region data
    
    dat = ncvar_get( file, varid="PRECL", start = c(89,96,1), count = c(9,10,-1) )
    
    file  = nc_open( paste0( "/scratch/xtan/LENS_single_forcing/PRECC/b.e11.B20TRLENS_RCP85.f09_g16.xbmb.",
                             models[[i_forcing]][i_model], ".cam.h0.PRECC.", years, ".nc"))
    
    dat1 = ncvar_get( file, varid="PRECC", start = c(89,96,1), count = c(9,10,-1) )
    
    precip[i_model,,,] = dat + dat1
    
    file = nc_open(filename = paste0("/scratch/xtan/LENS_single_forcing/spi/bmb/BMB_spi_12_",models[[i_forcing]][i_model],".nc"))
    # time <- ncvar_get(file, varid = "time")
    # lon <- ncvar_get(file, varid = "lon")
    # lat <- ncvar_get(file, varid = "lat")
    
    spi[i_model,,,] = ncvar_get(file, varid="spi", start = c(89,96,1), count = c(9,10,-1) )
  }
  
  save(precip, spi, file = paste0("precip_spi_bmb.RData"))
}

for(i_forcing in 3){  ####### GHG forcing 
 precip = array(NA, dim = c(length(models[[i_forcing]]),9,10,1932))

  for(i_model in 1:length(models[[i_forcing]])){
    years = "192001-200512"
    
    file = nc_open( paste0( "/scratch/xtan/LENS_single_forcing/PRECL/b.e11.B20TRLENS_RCP85.f09_g16.xghg.",
                            models[[i_forcing]][i_model],".cam.h0.PRECL.", years,".nc" ))
    
    #  spi[lon,lat,time]
    ## stduy region (part of south China) lat 20:28 (117:126) , lon, 110-120 (88:96)
    # precip_south_china = dat(lat|197:212, lon|131:144,  time|:, height|0)    ; extract the study region data
    
    dat = ncvar_get( file, varid="PRECL", start = c(89,96,1), count = c(9,10,-1) )
    
    file  = nc_open( paste0( "/scratch/xtan/LENS_single_forcing/PRECC/b.e11.B20TRLENS_RCP85.f09_g16.xghg.",
                             models[[i_forcing]][i_model], ".cam.h0.PRECC.", years, ".nc"))
    
    dat1 = ncvar_get( file, varid="PRECC", start = c(89,96,1), count = c(9,10,-1) )
    
    precip1 = dat + dat1
    
    years = "200601-208012"
    
    file = nc_open( paste0( "/scratch/xtan/LENS_single_forcing/PRECL/b.e11.B20TRLENS_RCP85.f09_g16.xghg.",
                            models[[i_forcing]][i_model],".cam.h0.PRECL.", years,".nc" ))
    
    #  spi[lon,lat,time]
    ## stduy region (part of south China) lat 20:28 (117:126) , lon, 110-120 (88:96)
    # precip_south_china = dat(lat|197:212, lon|131:144,  time|:, height|0)    ; extract the study region data
    
    dat = ncvar_get( file, varid="PRECL", start = c(89,96,1), count = c(9,10,-1) )
    
    file  = nc_open( paste0( "/scratch/xtan/LENS_single_forcing/PRECC/b.e11.B20TRLENS_RCP85.f09_g16.xghg.",
                             models[[i_forcing]][i_model], ".cam.h0.PRECC.", years, ".nc"))
    
    dat1 = ncvar_get( file, varid="PRECC", start = c(89,96,1), count = c(9,10,-1) )
    
    precip2 = dat + dat1
    
    precip[i_model,,,] = abind(precip1, precip2, along = 3)
    
    # file = nc_open(filename = paste0("/scratch/xtan/LENS_single_forcing/spi/GHG/GHG_spi_12_",models[[i_forcing]][i_model],".nc"))
    # # time <- ncvar_get(file, varid = "time")
    # # lon <- ncvar_get(file, varid = "lon")
    # # lat <- ncvar_get(file, varid = "lat")
    # 
    # spi[i_model,,,] = ncvar_get(file, varid="spi", start = c(89,96,1), count = c(9,10,-1) )
  }
  
  save(precip, file = paste0("precip_spi_GHG.RData"))
}

for(i_forcing in 4){  ####### HIST forcing 
  
  precip = array(NA, dim = c(length(models[[i_forcing]]),9,10,2172))
  years = "192001-200512"
  
  for(i_model in c(1:15,17:32)){

    years = "192001-200512"
    
    file = nc_open( paste0( "/scratch/xtan/LENS_single_forcing/PRECL/b.e11.B20TRC5CNBDRD.f09_g16.",
                            models[[i_forcing]][i_model],".cam.h0.PRECL.", years,".nc" ))
    
    #  spi[lon,lat,time]
    ## stduy region (part of south China) lat 20:28 (117:126) , lon, 110-120 (88:96)
    # precip_south_china = dat(lat|197:212, lon|131:144,  time|:, height|0)    ; extract the study region data
    
    dat = ncvar_get( file, varid="PRECL", start = c(89,96,1), count = c(9,10,-1) )
    
    file  = nc_open( paste0( "/scratch/xtan/LENS_single_forcing/PRECC/b.e11.B20TRC5CNBDRD.f09_g16.",
                             models[[i_forcing]][i_model], ".cam.h0.PRECC.", years, ".nc"))
    
    dat1 = ncvar_get( file, varid="PRECC", start = c(89,96,1), count = c(9,10,-1) )
    
    precip1 = dat + dat1
    
    years = "200601-208012"
    
    file = nc_open( paste0( "/scratch/xtan/LENS_single_forcing/PRECL/b.e11.BRCP85C5CNBDRD.f09_g16.",
                            models[[i_forcing]][i_model],".cam.h0.PRECL.", years,".nc" ))
    
    #  spi[lon,lat,time]
    ## stduy region (part of south China) lat 20:28 (117:126) , lon, 110-120 (88:96)
    # precip_south_china = dat(lat|197:212, lon|131:144,  time|:, height|0)    ; extract the study region data
    
    dat = ncvar_get( file, varid="PRECL", start = c(89,96,1), count = c(9,10,-1) )
    
    file  = nc_open( paste0( "/scratch/xtan/LENS_single_forcing/PRECC/b.e11.BRCP85C5CNBDRD.f09_g16.",
                             models[[i_forcing]][i_model], ".cam.h0.PRECC.", years, ".nc"))
    
    dat1 = ncvar_get( file, varid="PRECC", start = c(89,96,1), count = c(9,10,-1) )
    
    precip2 = dat + dat1
    
    years = "208101-210012"
    
    file = nc_open( paste0( "/scratch/xtan/LENS_single_forcing/PRECL/b.e11.BRCP85C5CNBDRD.f09_g16.",
                            models[[i_forcing]][i_model],".cam.h0.PRECL.", years,".nc" ))
    
    #  spi[lon,lat,time]
    ## stduy region (part of south China) lat 20:28 (117:126) , lon, 110-120 (88:96)
    # precip_south_china = dat(lat|197:212, lon|131:144,  time|:, height|0)    ; extract the study region data
    
    dat = ncvar_get( file, varid="PRECL", start = c(89,96,1), count = c(9,10,-1) )
    
    file  = nc_open( paste0( "/scratch/xtan/LENS_single_forcing/PRECC/b.e11.BRCP85C5CNBDRD.f09_g16.",
                             models[[i_forcing]][i_model], ".cam.h0.PRECC.", years, ".nc"))
    
    dat1 = ncvar_get( file, varid="PRECC", start = c(89,96,1), count = c(9,10,-1) )
    
    precip3 = dat + dat1
    
    precip[i_model,,,] = abind(precip1, precip2, precip3, along = 3)
    
    # file = nc_open(filename = paste0("/scratch/xtan/LENS_single_forcing/spi/historical/HIST_spi_12_",
    #                                  models[[i_forcing]][i_model],".nc"))
    # # time <- ncvar_get(file, varid = "time")
    # # lon <- ncvar_get(file, varid = "lon")
    # # lat <- ncvar_get(file, varid = "lat")
    # 
    # spi[i_model,,,] = ncvar_get(file, varid="spi", start = c(89,96,1), count = c(9,10,-1) )
  }
  
  for(i_model in 33:length(models[[i_forcing]])){
    
    years = "192001-200512"
    
    file = nc_open( paste0( "/scratch/xtan/LENS_single_forcing/PRECL/b.e11.B20TRC5CNBDRD.f09_g16.",
                            models[[i_forcing]][i_model],".cam.h0.PRECL.", years,".nc" ))
    
    #  spi[lon,lat,time]
    ## stduy region (part of south China) lat 20:28 (117:126) , lon, 110-120 (88:96)
    # precip_south_china = dat(lat|197:212, lon|131:144,  time|:, height|0)    ; extract the study region data
    
    dat = ncvar_get( file, varid="PRECL", start = c(89,96,1), count = c(9,10,-1) )
    
    file  = nc_open( paste0( "/scratch/xtan/LENS_single_forcing/PRECC/b.e11.B20TRC5CNBDRD.f09_g16.",
                             models[[i_forcing]][i_model], ".cam.h0.PRECC.", years, ".nc"))
    
    dat1 = ncvar_get( file, varid="PRECC", start = c(89,96,1), count = c(9,10,-1) )
    
    precip1 = dat + dat1
    
    years = "200601-210012"
    
    file = nc_open( paste0( "/scratch/xtan/LENS_single_forcing/PRECL/b.e11.BRCP85C5CNBDRD.f09_g16.",
                            models[[i_forcing]][i_model],".cam.h0.PRECL.", years,".nc" ))
    
    #  spi[lon,lat,time]
    ## stduy region (part of south China) lat 20:28 (117:126) , lon, 110-120 (88:96)
    # precip_south_china = dat(lat|197:212, lon|131:144,  time|:, height|0)    ; extract the study region data
    
    dat = ncvar_get( file, varid="PRECL", start = c(89,96,1), count = c(9,10,-1) )
    
    file  = nc_open( paste0( "/scratch/xtan/LENS_single_forcing/PRECC/b.e11.BRCP85C5CNBDRD.f09_g16.",
                             models[[i_forcing]][i_model], ".cam.h0.PRECC.", years, ".nc"))
    
    dat1 = ncvar_get( file, varid="PRECC", start = c(89,96,1), count = c(9,10,-1) )
    
    precip2 = dat + dat1
    
    precip[i_model,,,] = abind(precip1, precip2, along = 3)
    
    # file = nc_open(filename = paste0("/scratch/xtan/LENS_single_forcing/spi/historical/HIST_spi_12_",
    #                                  models[[i_forcing]][i_model],".nc"))
    # # time <- ncvar_get(file, varid = "time")
    # # lon <- ncvar_get(file, varid = "lon")
    # # lat <- ncvar_get(file, varid = "lat")
    # 
    # spi[i_model,,,] = ncvar_get(file, varid="spi", start = c(89,96,1), count = c(9,10,-1) )
  }
  
  save(precip, file = paste0("precip_spi_hist.RData"))
}

for(i_forcing in 5){  ####### HIST forcing 
   precip = array(NA, dim = c(15,9,10,1320))
  years = "192001-202912"
  
  for(i_model in 1:length(models[[i_forcing]])){
    
    file = nc_open( paste0( "/scratch/xtan/LENS_single_forcing/PRECL/b.e11.B20TRLENS_RCP85.f09_g16.xlulc.",
                            models[[i_forcing]][i_model],".cam.h0.PRECL.", years,".nc" ))
    
    #  spi[lon,lat,time]
    ## stduy region (part of south China) lat 20:28 (117:126) , lon, 110-120 (88:96)
    # precip_south_china = dat(lat|197:212, lon|131:144,  time|:, height|0)    ; extract the study region data
    
    dat = ncvar_get( file, varid="PRECL", start = c(89,96,1), count = c(9,10,-1) )
    
    file  = nc_open( paste0( "/scratch/xtan/LENS_single_forcing/PRECC/b.e11.B20TRLENS_RCP85.f09_g16.xlulc.",
                             models[[i_forcing]][i_model], ".cam.h0.PRECC.", years, ".nc"))
    
    dat1 = ncvar_get( file, varid="PRECC", start = c(89,96,1), count = c(9,10,-1) )
    
    precip[i_model,,,] = dat + dat1
    
    # file = nc_open(filename = paste0("/scratch/xtan/LENS_single_forcing/spi/lucl/LULC_spi_12_",models[[i_forcing]][i_model],".nc"))
    # # time <- ncvar_get(file, varid = "time")
    # # lon <- ncvar_get(file, varid = "lon")
    # # lat <- ncvar_get(file, varid = "lat")
    # 
    # spi[i_model,,,] = ncvar_get(file, varid="spi", start = c(89,96,1), count = c(9,10,-1) )
  }
  
  save(precip, file = paste0("precip_spi_lulc.RData"))
}


#######################################################################################################
load(file = paste0("precip_spi_hist.RData"))

precip =  precip * 86400 * 1000 * 30
ref.data = precip[,,,493:852] ### reference precipitation 1960-1990
ref.data = apply(ref.data, MARGIN = c(2,3,4), FUN = mean, na.rm=T)
spi_hist = spi_calc(precip=precip, ref.data= ref.data)
spi_hist = spi_hist[c(1:15,17:39),,,]  # lack of data in model 16
# spi_hist = apply(X=spi_hist, MARGIN = c(1,4), FUN = mean, na.rm=T )
# spi_hist = spi_hist[c(1:15,17:39),]
filter = 12
spi_hist_roll = apply(spi_hist, FUN = roll_mean, MARGIN = c(1,2,3), n=filter, fill=NA )
spi_hist_quantile = apply(spi_hist_roll, FUN = quantile, MARGIN = c(1), 
                          probs=c(0.1,0.5,0.9), na.rm=T )
plot(spi_hist_quantile[2,])

load(file = paste0("precip_spi_aer.RData"))
precip =  precip * 86400 * 1000 * 30
spi_aer = spi_calc(precip=precip, ref.data= ref.data)
filter = 12
spi_aer_roll = apply(spi_aer, FUN = roll_mean, MARGIN = c(1,2,3), n=filter, fill=NA )
spi_aer_quantile = apply(spi_aer_roll, FUN = quantile, MARGIN = c(1), 
                         probs=c(0.1,0.5,0.9), na.rm=T )
plot(spi_aer_quantile[2,])

load(file = paste0("precip_spi_GHG.RData"))
precip =  precip * 86400 * 1000 * 30
spi_ghg = spi_calc(precip=precip, ref.data= ref.data)
filter = 12
spi_ghg_roll = apply(spi_ghg, FUN = roll_mean, MARGIN = c(1,2,3), n=filter, fill=NA )
spi_ghg_quantile = apply(spi_ghg_roll, FUN = quantile, MARGIN = c(1), 
                         probs=c(0.1,0.5,0.9), na.rm=T )
plot(spi_ghg_quantile[2,])

load(file = paste0("precip_spi_bmb.RData"))
precip =  precip * 86400 * 1000 * 30
spi_bmb = spi_calc(precip=precip, ref.data= ref.data)
filter = 12
spi_bmb_roll = apply(spi_bmb, FUN = roll_mean, MARGIN = c(1,2,3), n=filter, fill=NA )
spi_bmb_quantile = apply(spi_bmb_roll, FUN = quantile, MARGIN = c(1), 
                         probs=c(0.1,0.5,0.9), na.rm=T )
plot(spi_bmb_quantile[2,])

load(file = paste0("precip_spi_lulc.RData"))
precip =  precip * 86400 * 1000 * 30
spi_lulc = spi_calc(precip=precip, ref.data= ref.data)
filter = 12
spi_lulc_roll = apply(spi_lulc, FUN = roll_mean, MARGIN = c(1,2,3), n=filter, fill=NA )
spi_lulc_quantile = apply(spi_lulc_roll, FUN = quantile, MARGIN = c(1), 
                         probs=c(0.1,0.5,0.9), na.rm=T )
plot(spi_lulc_quantile[2,])

save(spi_lulc, spi_lulc_quantile,
     spi_hist, spi_hist_quantile, 
     spi_aer, spi_aer_quantile,
     spi_bmb, spi_bmb_quantile, 
     spi_ghg, spi_ghg_quantile, file = "spi_ref_1960-1990.RData")

dat = NULL
dat = rbind(dat, cbind( year=seq(from=1920, to=2100+11/12, by=1/12 ), 
                         ymean = spi_hist_quantile[2,] , 
                         ymin = spi_hist_quantile[1,] ,  
                         ymax = spi_hist_quantile[3,] ,  
                         forcing = rep("Hist", 2172 ) ) ) 

dat = rbind(dat, cbind( year=seq(from=1920, to=2080+11/12, by=1/12 ), 
                        ymean = spi_aer_quantile[2,] , 
                        ymin = spi_aer_quantile[1,] ,  
                        ymax = spi_aer_quantile[3,] ,  
                        forcing = rep("AER", 1932 ) ) ) 

dat = rbind(dat, cbind( year=seq(from=1920, to=2080+11/12, by=1/12 ), 
                        ymean = spi_ghg_quantile[2,] , 
                        ymin = spi_ghg_quantile[1,] ,  
                        ymax = spi_ghg_quantile[3,] ,  
                        forcing = rep("ghg", 1932 ) ) ) 

dat = rbind(dat, cbind( year=seq(from=1920, to=2029+11/12, by=1/12 ), 
                        ymean = spi_bmb_quantile[2,] , 
                        ymin = spi_bmb_quantile[1,] ,  
                        ymax = spi_bmb_quantile[3,] ,  
                        forcing = rep("BMB", 1320 ) ) ) 

dat = rbind(dat, cbind( year=seq(from=1920, to=2029+11/12, by=1/12 ), 
                        ymean = spi_lulc_quantile[2,] , 
                        ymin = spi_lulc_quantile[1,] ,  
                        ymax = spi_lulc_quantile[3,] ,  
                        forcing = rep("LULC", 1320 ) ) ) 

dat = as.data.frame(dat)
dat$year = as.numeric(dat$year)
dat$ymean = as.numeric(dat$ymean)
dat$ymax = as.numeric(dat$ymax)
dat$ymin = as.numeric(dat$ymin)

p = ggplot(dat) +
  geom_ribbon(  mapping = aes(x=year, ymax = ymax, ymin = ymin, fill=forcing), alpha=0.3, size=0.1 ) +
  geom_line( mapping = aes(x=year, y = ymean, col=forcing), size=0.5 ) +
  theme(panel.background = element_rect(fill="white",color="black",size=0.3),
        panel.grid = element_blank(),legend.position='bottom',
        axis.title=element_text(size=10),axis.text=element_text(size=8),
        axis.ticks=element_line(size=0.1))+
  #stat_ecdf(alpha=0.55, aes(fwi)) + 
  geom_hline(aes(yintercept=0), colour="black",size=0.2) +
  scale_x_continuous(breaks= seq(1920,2100,10) ) +
  scale_y_continuous(breaks=seq(-2,2,0.4)) 
ggsave(filename = "forcings_spi_12_1.pdf",p,width = 8,height = 5)


dat = NULL
for(i_model in 1:19){
  aaa = as.vector(spi_aer[i_model,,,])
  dat = rbind (dat, cbind ( aaa, rep(i_model, length(aaa)), rep("AER", length(aaa) ) ) )  
} 
for(i_model in 1:38){
  aaa = as.vector(spi_hist[i_model,,,])
  dat = rbind (dat, cbind ( aaa, rep(i_model, length(aaa)), rep("HIST", length(aaa) ) ) )  
} 
for(i_model in 1:20){
  aaa = as.vector(spi_ghg[i_model,,,])
  dat = rbind (dat, cbind ( aaa, rep(i_model, length(aaa)), rep("GHG", length(aaa) ) ) )  
} 
for(i_model in 1:5){
  aaa = as.vector(spi_bmb[i_model,,,])
  dat = rbind (dat, cbind ( aaa, rep(i_model, length(aaa)), rep("BMB", length(aaa) ) ) )  
} 
for(i_model in 1:5){
  aaa = as.vector(spi_lulc[i_model,,,])
  dat = rbind (dat, cbind ( aaa, rep(i_model, length(aaa)), rep("ULC", length(aaa) ) ) )  
} 
colnames(dat) = c("SPI", "model", "forcing")
dat = as.data.frame(dat)
dat$SPI = as.numeric(dat$SPI)

p = ggplot(dat) +
  geom_density ( mapping = aes(x= SPI, col= forcing, group=model ), size=0.05, alpha = 0.5) +
  geom_density ( mapping = aes(x= SPI, col= forcing), size=0.5 ) +
  theme(panel.background = element_rect(fill="white",color="black",size=0.3),
        panel.grid = element_blank(),legend.position='bottom',
        axis.title=element_text(size=14),axis.text=element_text(size=14),
        axis.ticks=element_line(size=0.2)) +
  xlim (-3.0, 3.0)
  # scale_y_continuous(trans = 'log10')  +
  # scale_x_continuous(trans = 'log10')
ggsave(filename = "SPI_forcing_density.pdf",p,width = 4,height = 8)

############## risk ratio analyses ###############

# ## XGHG
# 
# threshold = -1.5
# 
# library(fmsb)
# rr = matrix(NA, 3, 161)
#   for(i_year in 1:161){
#     aaa= as.vector(spi_ghg[,,,((i_year-1)*12+1):(i_year*12)])
#     aaa = aaa[!is.na(aaa)]
#     occ_ghg = length(aaa[aaa < threshold])
#     bbb= as.vector(spi_hist[,,,((i_year-1)*12+1):(i_year*12)])
#     bbb = bbb[!is.na(bbb)]
#     occ_hist = length(bbb[bbb < threshold])
#     rr1 = riskratio(X=occ_hist, Y=occ_ghg,
#                     m1=length(bbb), m2=length(aaa), conf.level=0.95, p.calc.by.independence=TRUE)
#     rr[1, ] = rr1$estimate 
#     rr[2:3, ] = rr1$conf.int
# }
# rr_ghg = rr
# 
# ## XGHG
# library(fmsb)
# rr = matrix(NA, 3, 161)
# for(i_year in 1:161){
#   aaa= as.vector(spi_ghg[,,,((i_year-1)*12+1):(i_year*12)])
#   aaa = aaa[!is.na(aaa)]
#   occ_ghg = length(aaa[aaa < threshold])
#   bbb= as.vector(spi_hist[,,,((i_year-1)*12+1):(i_year*12)])
#   bbb = bbb[!is.na(bbb)]
#   occ_hist = length(bbb[bbb < threshold ])
#   rr1 = riskratio(X=occ_hist, Y=occ_ghg,
#                   m1=length(bbb), m2=length(aaa), conf.level=0.95, p.calc.by.independence=TRUE)
#   rr[1, i_year] = rr1$estimate 
#   rr[2:3, i_year ] = rr1$conf.int
# }
# rr_ghg = rr
# 
# ## XAER
# library(fmsb)
# rr = matrix(NA, 3, 161)
# for(i_year in 1:161){
#   aaa= as.vector(spi_aer[,,,((i_year-1)*12+1):(i_year*12)])
#   aaa = aaa[!is.na(aaa)]
#   occ_aer = length(aaa[aaa < threshold ])
#   bbb= as.vector(spi_hist[,,,((i_year-1)*12+1):(i_year*12)])
#   bbb = bbb[!is.na(bbb)]
#   occ_hist = length(bbb[bbb < threshold ])
#   rr1 = riskratio(X=occ_hist, Y=occ_aer,
#                   m1=length(bbb), m2=length(aaa), conf.level=0.95, p.calc.by.independence=TRUE)
#   rr[1, i_year] = rr1$estimate 
#   rr[2:3, i_year ] = rr1$conf.int
# }
# rr_aer = rr
# 
# ## XBMB
# library(fmsb)
# rr = matrix(NA, 3, 110)
# for(i_year in 1:110){
#   aaa= as.vector(spi_bmb[,,,((i_year-1)*12+1):(i_year*12)])
#   aaa = aaa[!is.na(aaa)]
#   occ_bmb = length(aaa[aaa < threshold])
#   bbb= as.vector(spi_hist[,,,((i_year-1)*12+1):(i_year*12)])
#   bbb = bbb[!is.na(bbb)]
#   occ_hist = length(bbb[bbb < threshold ])
#   rr1 = riskratio(X=occ_hist, Y=occ_bmb,
#                   m1=length(bbb), m2=length(aaa), conf.level=0.95, p.calc.by.independence=TRUE)
#   rr[1, i_year] = rr1$estimate 
#   rr[2:3, i_year ] = rr1$conf.int
# }
# rr_bmb = rr
# 
# ## XLULC
# library(fmsb)
# rr = matrix(NA, 3, 110)
# for(i_year in 1:110){
#   aaa= as.vector(spi_lulc[,,,((i_year-1)*12+1):(i_year*12)])
#   aaa = aaa[!is.na(aaa)]
#   occ_lulc = length(aaa[aaa < threshold ])
#   bbb= as.vector(spi_hist[,,,((i_year-1)*12+1):(i_year*12)])
#   bbb = bbb[!is.na(bbb)]
#   occ_hist = length(bbb[bbb < threshold ])
#   rr1 = riskratio(X=occ_hist, Y=occ_lulc,
#                   m1=length(bbb), m2=length(aaa), conf.level=0.95, p.calc.by.independence=TRUE)
#   rr[1, i_year] = rr1$estimate 
#   rr[2:3, i_year ] = rr1$conf.int
# }
# rr_lulc = rr
# 
# save(rr_ghg, rr_aer, rr_bmb, rr_lulc, file = "risk_ratio_LENS.RData")
# 
# for(i in 1:3){
#   rr_ghg[i,] [is.na(rr_ghg[i,]) | is.infinite(rr_ghg[i,])  ] = 
#     mean (rr_ghg[i,] [ ! (is.na(rr_ghg[i,]) | is.infinite(rr_ghg[i,]))  ])
#   rr_aer[i,] [is.na(rr_aer[i,]) | is.infinite(rr_aer[i,])  ] = 
#     mean (rr_aer[i,] [ ! (is.na(rr_aer[i,]) | is.infinite(rr_aer[i,]))  ])
#   rr_bmb[i,] [is.na(rr_bmb[i,]) | is.infinite(rr_bmb[i,])  ] = 
#     mean (rr_bmb[i,] [ ! (is.na(rr_bmb[i,]) | is.infinite(rr_bmb[i,]))  ])
#   rr_lulc[i,] [is.na(rr_lulc[i,]) | is.infinite(rr_lulc[i,])  ] = 
#     mean (rr_lulc[i,] [ ! (is.na(rr_lulc[i,]) | is.infinite(rr_lulc[i,]))  ])
# }

dat = NULL

dat = rbind(dat, cbind( year=1920:2080, 
                        ymean =rr_aer[2,] , 
                        ymin = rr_aer[1,] ,  
                        ymax = rr_aer[3,] ,  
                        forcing = rep("AER", 161 ) ) ) 

dat = rbind(dat, cbind( year=1920:2080, 
                        ymean =rr_ghg[2,] , 
                        ymin = rr_ghg[1,] ,  
                        ymax = rr_ghg[3,] ,  
                        forcing = rep("GHG", 161 ) ) ) 

dat = rbind(dat, cbind( year=1920:2029, 
                        ymean =rr_bmb[2,] , 
                        ymin = rr_bmb[1,] ,  
                        ymax = rr_bmb[3,] ,  
                        forcing = rep("BMB", 110 ) ) ) 

dat = rbind(dat, cbind( year=1920:2029, 
                        ymean =rr_lulc[2,] , 
                        ymin = rr_lulc[1,] ,  
                        ymax = rr_lulc[3,] ,  
                        forcing = rep("LULC", 110 ) ) ) 

dat = as.data.frame(dat)
dat$year = as.numeric(dat$year)
dat$ymean = as.numeric(dat$ymean)
dat$ymax = as.numeric(dat$ymax)
dat$ymin = as.numeric(dat$ymin)

p = ggplot(dat) +
  geom_ribbon(  mapping = aes(x=year, ymax = ymax, ymin = ymin, fill=forcing), alpha=0.3, size=0.1 ) +
  geom_line( mapping = aes(x=year, y = ymean, col=forcing), size=0.5 ) +
  theme(panel.background = element_rect(fill="white",color="black",size=0.3),
        panel.grid = element_blank(),legend.position='bottom',
        axis.title=element_text(size=10),axis.text=element_text(size=8),
        axis.ticks=element_line(size=0.1))+
  #stat_ecdf(alpha=0.55, aes(fwi)) + 
  ylim(0,5) +
  geom_hline(aes(yintercept=0), colour="black",size=0.2) +
  scale_x_continuous(breaks= seq(1920,2100,10) ) 
  # scale_y_continuous(breaks=seq(-5,5,0.4)) 
ggsave(filename = "risk_ratio.pdf",p,width = 8,height = 5)

rr = function(spi_hist, spi_nat, threshold, wet=F){
  library(fmsb)
  dd1 = dim(spi_hist)
  dd2 = dim(spi_nat)
  dd = min(dd1[4], dd2[4])/12
  dd_model = min( dd1[1], dd2[1] )
  rr = matrix(NA, dd_model, dd)
  rr_mean = NULL
  for(i_model in 1:dd_model){
    for(i_year in 1:dd){
      if(wet){
        aaa= as.vector(spi_hist[i_model,,,((i_year-1)*12+1):(i_year*12)])
        aaa = aaa[!is.na(aaa)]
        occ_nat = length(aaa[aaa > threshold])
        bbb= as.vector(spi_nat[i_model,,,((i_year-1)*12+1):(i_year*12)])
        bbb = bbb[!is.na(bbb)]
        occ_hist = length(bbb[bbb > threshold])
        rr1 = riskratio(X=occ_hist, Y=occ_nat,
                        m1=length(bbb), m2=length(aaa), conf.level=0.95, p.calc.by.independence=TRUE)
        rr[i_model, i_year] = rr1$estimate 
      }
      else{
        aaa= as.vector(spi_hist[i_model,,,((i_year-1)*12+1):(i_year*12)])
        aaa = aaa[!is.na(aaa)]
        occ_nat = length(aaa[aaa < threshold])
        bbb= as.vector(spi_nat[i_model,,,((i_year-1)*12+1):(i_year*12)])
        bbb = bbb[!is.na(bbb)]
        occ_hist = length(bbb[bbb < threshold])
        rr1 = riskratio(X=occ_hist, Y=occ_nat,
                        m1=length(bbb), m2=length(aaa), conf.level=0.95, p.calc.by.independence=TRUE)
        rr[i_model, i_year] = rr1$estimate 
      }

    }
    if(wet){
      aaa= as.vector(spi_hist[i_model,,,1:(dd*12)])
      aaa = aaa[!is.na(aaa)]
      occ_nat = length(aaa[aaa > threshold])
      bbb= as.vector(spi_nat[i_model,,,1:(dd*12)])
      bbb = bbb[!is.na(bbb)]
      occ_hist = length(bbb[bbb > threshold])
      rr1 = riskratio(X=occ_hist, Y=occ_nat,
                      m1=length(bbb), m2=length(aaa), conf.level=0.95, p.calc.by.independence=TRUE)
      
      rr_mean [i_model] = rr1$estimate
    }
    else{
      aaa= as.vector(spi_hist[i_model,,,1:(dd*12)])
      aaa = aaa[!is.na(aaa)]
      occ_nat = length(aaa[aaa < threshold])
      bbb= as.vector(spi_nat[i_model,,,1:(dd*12)])
      bbb = bbb[!is.na(bbb)]
      occ_hist = length(bbb[bbb < threshold])
      rr1 = riskratio(X=occ_hist, Y=occ_nat,
                      m1=length(bbb), m2=length(aaa), conf.level=0.95, p.calc.by.independence=TRUE)
      
      rr_mean [i_model] = rr1$estimate
    }
  }
  return( risk_ration = list (rr=rr, rr_mean=rr_mean) )
}

load (file = "spi_ref_1960-1990.RData")

rr_ghg_10 = rr(spi_hist, spi_ghg, threshold = -1.0)
rr_ghg_15 = rr(spi_hist, spi_ghg, threshold = -1.5)
rr_ghg_20 = rr(spi_hist, spi_ghg, threshold = -2.0)
rr_ghg_25 = rr(spi_hist, spi_ghg, threshold = -2.5)
rr_ghg_30 = rr(spi_hist, spi_ghg, threshold = -3.0)

rr_ghg_wet_10 = rr(spi_hist, spi_ghg, threshold = 1.0, wet = T)
rr_ghg_wet_15 = rr(spi_hist, spi_ghg, threshold = 1.5, wet = T)
rr_ghg_wet_20 = rr(spi_hist, spi_ghg, threshold = 2.0, wet = T)
rr_ghg_wet_25 = rr(spi_hist, spi_ghg, threshold = 2.5, wet = T)
rr_ghg_wet_30 = rr(spi_hist, spi_ghg, threshold = 3.0, wet = T)

rr_aer_10 = rr(spi_hist, spi_aer, threshold = -1.0)
rr_aer_15 = rr(spi_hist, spi_aer, threshold = -1.5)
rr_aer_20 = rr(spi_hist, spi_aer, threshold = -2.0)
rr_aer_25 = rr(spi_hist, spi_aer, threshold = -2.5)
rr_aer_30 = rr(spi_hist, spi_aer, threshold = -3.0)

rr_aer_wet_10 = rr(spi_hist, spi_aer, threshold = 1.0, wet = T)
rr_aer_wet_15 = rr(spi_hist, spi_aer, threshold = 1.5, wet = T)
rr_aer_wet_20 = rr(spi_hist, spi_aer, threshold = 2.0, wet = T)
rr_aer_wet_25 = rr(spi_hist, spi_aer, threshold = 2.5, wet = T)
rr_aer_wet_30 = rr(spi_hist, spi_aer, threshold = 3.0, wet = T)

spi_bmb = spi_bmb[1:5,,,]
rr_bmb_10 = rr(spi_hist, spi_bmb, threshold = -1.0)
rr_bmb_15 = rr(spi_hist, spi_bmb, threshold = -1.5)
rr_bmb_20 = rr(spi_hist, spi_bmb, threshold = -2.0)
rr_bmb_25 = rr(spi_hist, spi_bmb, threshold = -2.5)
rr_bmb_30 = rr(spi_hist, spi_bmb, threshold = -3.0)

rr_bmb_wet_10 = rr(spi_hist, spi_bmb, threshold = 1.0, wet = T)
rr_bmb_wet_15 = rr(spi_hist, spi_bmb, threshold = 1.5, wet = T)
rr_bmb_wet_20 = rr(spi_hist, spi_bmb, threshold = 2.0, wet = T)
rr_bmb_wet_25 = rr(spi_hist, spi_bmb, threshold = 2.5, wet = T)
rr_bmb_wet_30 = rr(spi_hist, spi_bmb, threshold = 3.0, wet = T)

spi_lulc = spi_lulc[1:5,,,]
rr_lulc_10 = rr(spi_hist, spi_lulc, threshold = -1.0)
rr_lulc_15 = rr(spi_hist, spi_lulc, threshold = -1.5)
rr_lulc_20 = rr(spi_hist, spi_lulc, threshold = -2.0)
rr_lulc_25 = rr(spi_hist, spi_lulc, threshold = -2.5)
rr_lulc_30 = rr(spi_hist, spi_lulc, threshold = -3.0)

rr_lulc_wet_10 = rr(spi_hist, spi_lulc, threshold = 1.0, wet = T)
rr_lulc_wet_15 = rr(spi_hist, spi_lulc, threshold = 1.5, wet = T)
rr_lulc_wet_20 = rr(spi_hist, spi_lulc, threshold = 2.0, wet = T)
rr_lulc_wet_25 = rr(spi_hist, spi_lulc, threshold = 2.5, wet = T)
rr_lulc_wet_30 = rr(spi_hist, spi_lulc, threshold = 3.0, wet = T)

dat = NULL
dat = rbind(dat, cbind (rr_10$rr_mean, rep("10", length(rr_10$rr_mean)),
                        rep("drought",length(rr_10$rr_mean)), rep("HadGEM3-GA6_HIST/NAT",length(rr_10$rr_mean)) ) ) 
dat = rbind(dat, cbind (rr_15$rr_mean, rep("15", length(rr_15$rr_mean)), 
                        rep("drought",length(rr_15$rr_mean)), rep("HadGEM3-GA6_HIST/NAT",length(rr_15$rr_mean)) ) )
dat = rbind(dat, cbind (rr_20$rr_mean, rep("20", length(rr_20$rr_mean)), 
                        rep("drought",length(rr_20$rr_mean)), rep("HadGEM3-GA6_HIST/NAT",length(rr_20$rr_mean)) ) )
dat = rbind(dat, cbind (rr_25$rr_mean, rep("25", length(rr_25$rr_mean)),
                        rep("drought",length(rr_25$rr_mean)), rep("HadGEM3-GA6_HIST/NAT",length(rr_25$rr_mean)) ) )
dat = rbind(dat, cbind (rr_30$rr_mean, rep("30", length(rr_30$rr_mean)), 
                        rep("drought",length(rr_30$rr_mean)), rep("HadGEM3-GA6_HIST/NAT",length(rr_30$rr_mean)) ) )

dat = rbind(dat, cbind (rr_wet_10$rr_mean, rep("10", length(rr_wet_10$rr_mean)),
                        rep("wet",length(rr_wet_10$rr_mean)), rep("HadGEM3-GA6_HIST/NAT",length(rr_wet_10$rr_mean)) ) )
dat = rbind(dat, cbind (rr_wet_15$rr_mean, rep("15", length(rr_wet_15$rr_mean)), 
                        rep("wet",length(rr_wet_15$rr_mean)), rep("HadGEM3-GA6_HIST/NAT",length(rr_wet_15$rr_mean)) ) )
dat = rbind(dat, cbind (rr_wet_20$rr_mean, rep("20", length(rr_wet_20$rr_mean)), 
                        rep("wet",length(rr_wet_20$rr_mean)), rep("HadGEM3-GA6_HIST/NAT",length(rr_wet_20$rr_mean)) ) )
dat = rbind(dat, cbind (rr_wet_25$rr_mean, rep("25", length(rr_wet_25$rr_mean)), 
                        rep("wet",length(rr_wet_25$rr_mean)) , rep("HadGEM3-GA6_HIST/NAT",length(rr_wet_25$rr_mean)) ) )
dat = rbind(dat, cbind (rr_wet_30$rr_mean, rep("30", length(rr_wet_30$rr_mean)), 
                        rep("wet",length(rr_wet_30$rr_mean)), rep("HadGEM3-GA6_HIST/NAT",length(rr_wet_30$rr_mean)) ) )

dat = rbind(dat, cbind (rr_ghg_10$rr_mean, rep("10", length(rr_ghg_10$rr_mean)),
                        rep("drought",length(rr_ghg_10$rr_mean)), rep("LENS_HIST/no_GHG",length(rr_ghg_10$rr_mean)) ) )
dat = rbind(dat, cbind (rr_ghg_15$rr_mean, rep("15", length(rr_ghg_15$rr_mean)), 
                        rep("drought",length(rr_ghg_15$rr_mean)), rep("LENS_HIST/no_GHG",length(rr_ghg_15$rr_mean)) ) )
dat = rbind(dat, cbind (rr_ghg_20$rr_mean, rep("20", length(rr_ghg_20$rr_mean)),
                        rep("drought",length(rr_ghg_20$rr_mean)), rep("LENS_HIST/no_GHG",length(rr_ghg_20$rr_mean)) ) )
dat = rbind(dat, cbind (rr_ghg_25$rr_mean, rep("25", length(rr_ghg_25$rr_mean)), 
                        rep("drought",length(rr_ghg_25$rr_mean)), rep("LENS_HIST/no_GHG",length(rr_ghg_25$rr_mean)) ) )
dat = rbind(dat, cbind (rr_ghg_30$rr_mean, rep("30", length(rr_ghg_30$rr_mean)),
                        rep("drought",length(rr_ghg_30$rr_mean)), rep("LENS_HIST/no_GHG",length(rr_ghg_30$rr_mean)) ) )

dat = rbind(dat, cbind (rr_ghg_wet_10$rr_mean, rep("10", length(rr_ghg_wet_10$rr_mean)), 
                        rep("wet",length(rr_ghg_wet_10$rr_mean)), rep("LENS_HIST/no_GHG",length(rr_ghg_wet_10$rr_mean)) ) )
dat = rbind(dat, cbind (rr_ghg_wet_15$rr_mean, rep("15", length(rr_ghg_wet_15$rr_mean)), 
                        rep("wet",length(rr_ghg_wet_15$rr_mean)), rep("LENS_HIST/no_GHG",length(rr_ghg_wet_15$rr_mean)) ) )
dat = rbind(dat, cbind (rr_ghg_wet_20$rr_mean, rep("20", length(rr_ghg_wet_20$rr_mean)), 
                        rep("wet",length(rr_ghg_wet_20$rr_mean)), rep("LENS_HIST/no_GHG",length(rr_ghg_wet_20$rr_mean)) ) )
dat = rbind(dat, cbind (rr_ghg_wet_25$rr_mean, rep("25", length(rr_ghg_wet_25$rr_mean)),
                        rep("wet",length(rr_ghg_wet_25$rr_mean)), rep("LENS_HIST/no_GHG",length(rr_ghg_wet_25$rr_mean)) ) )
dat = rbind(dat, cbind (rr_ghg_wet_30$rr_mean, rep("30", length(rr_ghg_wet_30$rr_mean)),
                        rep("wet",length(rr_ghg_wet_30$rr_mean)), rep("LENS_HIST/no_GHG",length(rr_ghg_wet_30$rr_mean)) ) )

dat = rbind(dat, cbind (rr_aer_10$rr_mean, rep("10", length(rr_aer_10$rr_mean)),
                        rep("drought",length(rr_aer_10$rr_mean)), rep("LENS_HIST/no_AER",length(rr_aer_10$rr_mean)) ) )
dat = rbind(dat, cbind (rr_aer_15$rr_mean, rep("15", length(rr_aer_15$rr_mean)),
                        rep("drought",length(rr_aer_15$rr_mean)), rep("LENS_HIST/no_AER",length(rr_aer_15$rr_mean)) ) )
dat = rbind(dat, cbind (rr_aer_20$rr_mean, rep("20", length(rr_aer_20$rr_mean)), 
                        rep("drought",length(rr_aer_20$rr_mean)), rep("LENS_HIST/no_AER",length(rr_aer_20$rr_mean)) ) )
dat = rbind(dat, cbind (rr_aer_25$rr_mean, rep("25", length(rr_aer_25$rr_mean)),
                        rep("drought",length(rr_aer_25$rr_mean)), rep("LENS_HIST/no_AER",length(rr_aer_25$rr_mean)) ) )
dat = rbind(dat, cbind (rr_aer_30$rr_mean, rep("30", length(rr_aer_30$rr_mean)),
                        rep("drought",length(rr_aer_30$rr_mean)), rep("LENS_HIST/no_AER",length(rr_aer_30$rr_mean)) ) )

dat = rbind(dat, cbind (rr_aer_wet_10$rr_mean, rep("10", length(rr_aer_wet_10$rr_mean)), 
                        rep("wet",length(rr_aer_wet_10$rr_mean)), rep("LENS_HIST/no_AER",length(rr_aer_wet_10$rr_mean)) ) )
dat = rbind(dat, cbind (rr_aer_wet_15$rr_mean, rep("15", length(rr_aer_wet_15$rr_mean)),
                        rep("wet",length(rr_aer_wet_15$rr_mean)), rep("LENS_HIST/no_AER",length(rr_aer_wet_15$rr_mean)) ) )
dat = rbind(dat, cbind (rr_aer_wet_20$rr_mean, rep("20", length(rr_aer_wet_20$rr_mean)),
                        rep("wet",length(rr_aer_wet_20$rr_mean)), rep("LENS_HIST/no_AER",length(rr_aer_wet_20$rr_mean)) ) )
dat = rbind(dat, cbind (rr_aer_wet_25$rr_mean, rep("25", length(rr_aer_wet_25$rr_mean)),
                        rep("wet",length(rr_aer_wet_25$rr_mean)), rep("LENS_HIST/no_AER",length(rr_aer_wet_25$rr_mean)) ) )
dat = rbind(dat, cbind (rr_aer_wet_30$rr_mean, rep("30", length(rr_aer_wet_30$rr_mean)),
                        rep("wet",length(rr_aer_wet_30$rr_mean)), rep("LENS_HIST/no_AER",length(rr_aer_wet_30$rr_mean)) ) )

dat = rbind(dat, cbind (rr_bmb_10$rr_mean, rep("10", length(rr_bmb_10$rr_mean)), 
                        rep("drought",length(rr_bmb_10$rr_mean)), rep("LENS_HIST/no_BMB",length(rr_bmb_10$rr_mean)) ) )
dat = rbind(dat, cbind (rr_bmb_15$rr_mean, rep("15", length(rr_bmb_15$rr_mean)),
                        rep("drought",length(rr_bmb_15$rr_mean)), rep("LENS_HIST/no_BMB",length(rr_bmb_15$rr_mean)) ) )
dat = rbind(dat, cbind (rr_bmb_20$rr_mean, rep("20", length(rr_bmb_20$rr_mean)),
                        rep("drought",length(rr_bmb_20$rr_mean)), rep("LENS_HIST/no_BMB",length(rr_bmb_20$rr_mean)) ) )
dat = rbind(dat, cbind (rr_bmb_25$rr_mean, rep("25", length(rr_bmb_25$rr_mean)),
                        rep("drought",length(rr_bmb_25$rr_mean)), rep("LENS_HIST/no_BMB",length(rr_bmb_25$rr_mean)) ) )
dat = rbind(dat, cbind (rr_bmb_30$rr_mean, rep("30", length(rr_bmb_30$rr_mean)),
                        rep("drought",length(rr_bmb_30$rr_mean)), rep("LENS_HIST/no_BMB",length(rr_bmb_30$rr_mean)) ) )

dat = rbind(dat, cbind (rr_bmb_wet_10$rr_mean, rep("10", length(rr_bmb_wet_10$rr_mean)), 
                        rep("wet",length(rr_bmb_wet_10$rr_mean)), rep("LENS_HIST/no_BMB",length(rr_bmb_wet_10$rr_mean)) ) )
dat = rbind(dat, cbind (rr_bmb_wet_15$rr_mean, rep("15", length(rr_bmb_wet_15$rr_mean)),
                        rep("wet",length(rr_bmb_wet_15$rr_mean)), rep("LENS_HIST/no_BMB",length(rr_bmb_wet_15$rr_mean)) ) )
dat = rbind(dat, cbind (rr_bmb_wet_20$rr_mean, rep("20", length(rr_bmb_wet_20$rr_mean)),
                        rep("wet",length(rr_bmb_wet_20$rr_mean)), rep("LENS_HIST/no_BMB",length(rr_bmb_wet_20$rr_mean)) ) )
dat = rbind(dat, cbind (rr_bmb_wet_25$rr_mean, rep("25", length(rr_bmb_wet_25$rr_mean)),
                        rep("wet",length(rr_bmb_wet_25$rr_mean)), rep("LENS_HIST/no_BMB",length(rr_bmb_wet_25$rr_mean)) ) )
dat = rbind(dat, cbind (rr_bmb_wet_30$rr_mean, rep("30", length(rr_bmb_wet_30$rr_mean)), 
                        rep("wet",length(rr_bmb_wet_30$rr_mean)), rep("LENS_HIST/no_BMB",length(rr_bmb_wet_30$rr_mean)) ) )

dat = rbind(dat, cbind (rr_lulc_10$rr_mean, rep("10", length(rr_lulc_10$rr_mean)), 
                        rep("drought",length(rr_lulc_10$rr_mean)), rep("LENS_HIST/no_LULC",length(rr_lulc_10$rr_mean)) ) )
dat = rbind(dat, cbind (rr_lulc_15$rr_mean, rep("15", length(rr_lulc_15$rr_mean)),
                        rep("drought",length(rr_lulc_15$rr_mean)), rep("LENS_HIST/no_LULC",length(rr_lulc_15$rr_mean)) ) )
dat = rbind(dat, cbind (rr_lulc_20$rr_mean, rep("20", length(rr_lulc_20$rr_mean)),
                        rep("drought",length(rr_lulc_20$rr_mean)), rep("LENS_HIST/no_LULC",length(rr_lulc_20$rr_mean)) ) )
dat = rbind(dat, cbind (rr_lulc_25$rr_mean, rep("25", length(rr_lulc_25$rr_mean)),
                        rep("drought",length(rr_lulc_25$rr_mean)) , rep("LENS_HIST/no_LULC",length(rr_lulc_25$rr_mean))) )
dat = rbind(dat, cbind (rr_lulc_30$rr_mean, rep("30", length(rr_lulc_30$rr_mean)),
                        rep("drought",length(rr_lulc_30$rr_mean)) , rep("LENS_HIST/no_LULC",length(rr_lulc_30$rr_mean))) )

dat = rbind(dat, cbind (rr_lulc_wet_10$rr_mean, rep("10", length(rr_lulc_wet_10$rr_mean)),
                        rep("wet",length(rr_lulc_wet_10$rr_mean)), rep("LENS_HIST/no_LULC",length(rr_lulc_wet_10$rr_mean)) ) )
dat = rbind(dat, cbind (rr_lulc_wet_15$rr_mean, rep("15", length(rr_lulc_wet_15$rr_mean)),
                        rep("wet",length(rr_lulc_wet_15$rr_mean)), rep("LENS_HIST/no_LULC",length(rr_lulc_wet_15$rr_mean)) ) )
dat = rbind(dat, cbind (rr_lulc_wet_20$rr_mean, rep("20", length(rr_lulc_wet_20$rr_mean)),
                        rep("wet",length(rr_lulc_wet_20$rr_mean)), rep("LENS_HIST/no_LULC",length(rr_lulc_wet_20$rr_mean)) ) )
dat = rbind(dat, cbind (rr_lulc_wet_25$rr_mean, rep("25", length(rr_lulc_wet_25$rr_mean)),
                        rep("wet",length(rr_lulc_wet_25$rr_mean)), rep("LENS_HIST/no_LULC",length(rr_lulc_wet_25$rr_mean)) ) )
dat = rbind(dat, cbind (rr_lulc_wet_30$rr_mean, rep("30", length(rr_lulc_wet_30$rr_mean)),
                        rep("wet",length(rr_lulc_wet_30$rr_mean)), rep("LENS_HIST/no_LULC",length(rr_lulc_wet_30$rr_mean)) ) )

dat = as.data.frame(dat)
colnames(dat) = c("rr","threshold", "regime", "forcing" )
dat$rr = as.numeric(dat$rr)

save(dat, file = "risk_ratio.RData")

load(file = "risk_ratio.RData")
dat_drought = dat[dat$regime %in% "drought",]
dat_wet = dat[dat$regime %in% "wet",]
p = ggplot(dat_drought) +
  geom_boxplot(mapping = aes(x= forcing, y=rr, colour = threshold ))  +
  geom_hline(yintercept = 1, size=0.3) +
  theme(panel.background = element_rect(fill="white",color="black",size=0.3),
        panel.grid = element_blank(),legend.position='bottom',
        axis.title=element_text(size=10),axis.text=element_text(size=8),
        axis.ticks=element_line(size=0.1)) +
  ylim(0,3)
 
ggsave(filename = "risk_ratio_drought.pdf",p,width = 6,height = 4)

p = ggplot(dat_wet) +
  geom_boxplot(mapping = aes(x= forcing, y=rr, colour = threshold ))  +
  geom_hline(yintercept = 1, size=0.3) +
  theme(panel.background = element_rect(fill="white",color="black",size=0.3),
        panel.grid = element_blank(),legend.position='bottom',
        axis.title=element_text(size=10),axis.text=element_text(size=8),
        axis.ticks=element_line(size=0.1)) +
  ylim(0,3)
ggsave(filename = "risk_ratio_wet.pdf",p,width = 6,height = 4)

#############################################################################################
################## function for retrieval drought characteristics ###########################
#############################################################################################

#### new function
#### run_theory() 
#### applied the run theory to a time series ()
## time_serie: a numeric vector with no NA values
## threshold: a numeric value in which the features (below) of run theory is measured

library(magrittr)

run_theory <- function(time_serie, threshold = -.5, wet = F){
  if(!wet){
    dataBase <- data.frame(time_serie = time_serie) %>%
      transform(masked = ifelse(time_serie >= threshold, 1, 0)) %>%
      transform(index = cumsum(masked), index_rev = cumsum(abs(masked-1)))
  }
  else{
    dataBase <- data.frame(time_serie = time_serie) %>%
      transform(masked = ifelse(time_serie <= threshold, 1, 0)) %>%
      transform(index = cumsum(masked), index_rev = cumsum(abs(masked-1)))
  }
  
  ####  Duration, Severity, Intesity ####  
  
  dataBase[dataBase$masked == 0, ] %>%
    by(., .$index, function(z){
      
      data.frame(D = dim(z)[1],
                 S = abs(sum(z$time_serie)),
                 I = abs(sum(z$time_serie))/dim(z)[1],
                 date_ini = row.names(z)[1],
                 date_fin = row.names(z)[nrow(z)]  )
      
    }) %>% do.call(rbind, .) -> df1
  
  
  ####  Interarrival ####  
  
  dataBase[dataBase$masked != 0, ] %>%
    by(., .$index_rev, function(z){
      
      data.frame(Int = dim(z)[1])
      
    }) %>% unlist() -> Int
  
  # first condition 
  
  # if ((dataBase$masked[1]) == 1){
  #   
  #   Int <- Int[-1] 
  #   
  # } else {
  #   
  #   Int <- Int
  #   
  # }
  # 
  # # second condition
  # 
  # if (dataBase$masked[length(dataBase$masked)] == 1) {
  #   
  #   Int <- Int + df1$D
  #   
  # } else { 
  #   
  #   n <- c(Int, 0 ) + df1$D
  #   Int <- n[-length(n)]
  #   
  # }
  # 
  
  return(res = cbind(Duration = as.numeric(df1$D),
                     Severity = as.numeric(df1$S),
                     Intesity = as.numeric(df1$I),
                     Date_Ini_Ev = as.character(df1$date_ini),
                     Date_Fin_Ev = as.character(df1$date_fin),
                     Interarrival = as.numeric(Int)) )
  
}

wet_dry_period = function(spi){
  n = dim(spi)
  wet = dry = NULL
  for(i_model in 1:n[1]){
    for(i_lon in 1:n[2]){
      for(i_lat in 1:n[3]){
        
        dat = spi[i_model,i_lon, i_lat, 12:n[4]]
        
        wet_features = run_theory(time_serie = dat,
                                  threshold = 0.5, wet = T)
        ll = nrow(wet_features)
        wet <- rbind(wet, cbind (wet_features, model = rep (i_model, ll ), lon = rep(i_lon, ll), lat = rep(i_lat, ll) ) )
        
        dry_features = run_theory(time_serie = dat,
                                  threshold = -.5, wet = F)
        ll = nrow(dry_features)
        dry <- rbind(dry, cbind (dry_features, model = rep (i_model, ll ), lon = rep(i_lon, ll), lat = rep(i_lat, ll) ) )
      }
    }
  }
  return(list(wet= wet, dry=dry))
}

wet = wet_dry_period(spi=spi_hist[c(1:15,17:39),,,])
save(wet, file = "HIST_drought_regimes.RData")
wet = wet_dry_period(spi=spi_aer)
save(wet,  file = "AER_drought_regimes.RData")
wet = wet_dry_period(spi=spi_ghg)
save(wet,  file = "GHG_drought_regimes.RData")
wet = wet_dry_period(spi=spi_bmb)
save(wet,  file = "BMB_drought_regimes.RData")
wet = wet_dry_period(spi=spi_lulc[1:5,,,])
save(wet,  file = "LULC_drought_regimes.RData")

###################################################################################################
###################################################################################################
###################################################################################################

load(file = "HIST_drought_regimes.RData")
dat = cbind(wet$dry, forcing = rep( "HIST", nrow(wet$dry) ) )
load(file = "AER_drought_regimes.RData")
dat = rbind(dat, cbind ( wet$dry , forcing = rep( "AER", nrow(wet$dry) ))  )
load(file = "GHG_drought_regimes.RData")
dat = rbind(dat, cbind ( wet$dry , forcing = rep( "GHG", nrow(wet$dry) ))  )
load(file = "BMB_drought_regimes.RData")
dat = rbind(dat, cbind ( wet$dry , forcing = rep( "BMB", nrow(wet$dry) ))  )
load(file = "LULC_drought_regimes.RData")
dat = rbind(dat, cbind ( wet$dry , forcing = rep( "LULC", nrow(wet$dry) ))  )

dat = data.frame(dat)
dat$Duration = as.numeric(dat$Duration) 
dat$Severity = as.numeric(dat$Severity)
dat$Intesity = as.numeric(dat$Intesity)
dat$Date_Ini_Ev = as.numeric(dat$Date_Ini_Ev)
dat$Date_Fin_Ev = as.numeric(dat$Date_Fin_Ev)

library(ggplot2)
p = ggplot(dat) +
  geom_density( mapping = aes(x= Duration, col= forcing, group=model), size=0.1 ) +
  geom_density ( mapping = aes(x= Duration, col= forcing), size=0.5 ) +
  theme(panel.background = element_rect(fill="white",color="black",size=0.3),
        panel.grid = element_blank(),legend.position='bottom',
        axis.title=element_text(size=10),axis.text=element_text(size=8),
        axis.ticks=element_line(size=0.1)) +
  xlim(0, 30)
ggsave(filename = "duration_density_single_forcing.pdf",p,width = 2,height = 3)

library(ggplot2)
p = ggplot(dat) +
  geom_density( mapping = aes(x= Severity, col= forcing, group=model ), size=0.1) +
  geom_density ( mapping = aes(x= Severity, col= forcing), size=0.5 ) +
  theme(panel.background = element_rect(fill="white",color="black",size=0.3),
        panel.grid = element_blank(),legend.position='bottom',
        axis.title=element_text(size=10),axis.text=element_text(size=8),
        axis.ticks=element_line(size=0.1)) +
  xlim (0, 50) 
ggsave(filename = "Severity_density_single_forcing.pdf",p,width = 2,height = 3)

library(ggplot2)
p = ggplot(dat) +
  geom_density( mapping = aes(x= Intesity, col= forcing, group=model ), size=0.1) +
  geom_density ( mapping = aes(x= Intesity, col= forcing), size=0.5 ) +
  theme(panel.background = element_rect(fill="white",color="black",size=0.3),
        panel.grid = element_blank(),legend.position='bottom',
        axis.title=element_text(size=10),axis.text=element_text(size=8),
        axis.ticks=element_line(size=0.1)) +
  xlim (0, 3) 
ggsave(filename = "Intensity_density_single_forcing.pdf",p,width = 2,height = 3)





