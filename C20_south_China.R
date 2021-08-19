
rm(list=ls())

setwd("F:/data/c20")

forcings = c("hist","NAT")

models = 1:15

climatology = array( read.table("climatology_c20.txt")[[1]],  dim=c(12,14,16) ) 
climatology = climatology * 86400 * 30


for(i_forcing in 1){  ####### hist forcing 
  anomalies = array(NA, dim =c( 15, 648,14,16) )
  for(i_model in 1:15 ){
    ####### data dimension lat 10 lon 9 time 1932
    dat = read.table(file = paste0("precipitation_anomalies_",forcings[i_forcing],"_",i_model, ".txt" ) )
    anomalies [i_model, , ,] = array(data = dat[[1]], dim=c( 648,14,16)) * 86400  * 30 
  }
  anomalies_hist = anomalies
}


for(i_forcing in 2){  ####### NAT forcing 
  anomalies = array(NA, dim =c( 15, 648,14,16) )
  for(i_model in 1:15 ){
    ####### data dimension lat 10 lon 9 time 1932
    dat = read.table(file = paste0("precipitation_anomalies_",forcings[i_forcing],"_",i_model, ".txt" ) )
    anomalies [i_model, , ,] = array(data = dat[[1]], dim=c( 648,14,16)) * 86400 * 30 
  }
  anomalies_NAT = anomalies
}


library(RcppRoll)

filter = 12*1

anomalies_hist_mean = apply(anomalies_hist, FUN = mean, MARGIN = c(1, 2) )
anomalies_hist_mean_roll = apply(anomalies_hist_mean, FUN = roll_mean, MARGIN = c(1), n=filter, fill=NA )
anomalies_hist_mean_quantile = apply(anomalies_hist_mean_roll, FUN = quantile, MARGIN = c(1), 
                                    probs=c(0.1,0.5,0.9), na.rm=T )
plot(anomalies_hist_mean_quantile[2,])

anomalies_NAT_mean = apply(anomalies_NAT, FUN = mean, MARGIN = c(1, 2) )
anomalies_NAT_mean_roll = apply(anomalies_NAT_mean, FUN = roll_mean, MARGIN = c(1), n=filter, fill=NA )
anomalies_NAT_mean_quantile = apply(anomalies_NAT_mean_roll, FUN = quantile, MARGIN = c(1), 
                                    probs=c(0.1,0.5,0.9), na.rm=T )
plot(anomalies_NAT_mean_quantile[2,])

climatology_mean = apply(climatology, FUN = mean, MARGIN = c(1) )

library(abind)
anomalies_mean_quantile = rbind(anomalies_hist_mean_quantile, anomalies_NAT_mean_quantile )

dat = NULL
for(i_forcing in 1:2){  ####### HIST forcing 
    dat = rbind( dat, cbind( year=seq(from=1960, to=2013+11/12, by=1/12 ), 
                             ymean = as.numeric( anomalies_mean_quantile[(i_forcing-1)*3+2,]) , 
                             ymin = as.numeric( anomalies_mean_quantile[(i_forcing-1)*3+1,]) ,  
                             ymax = as.numeric( anomalies_mean_quantile[(i_forcing-1)*3+3,]) ,  
                             foring = rep(forcings[i_forcing], 648 ) ) )
}
dat = as.data.frame(dat)
dat$year = as.numeric(dat$year)
dat$ymean = as.numeric(dat$ymean)
dat$ymax = as.numeric(dat$ymax)
dat$ymin = as.numeric(dat$ymin)

####### station data 
dat_station = read.csv("dat_station.csv")
dat_station_monthly = cbind( read.csv("monthly_anomalies_south_China.csv"), 
                             year=seq(from=1959, to=2021+02/12, by=1/12 ) )
dat_station_monthly$x = dat_station_monthly$x/3
dat_era5_monthly = cbind( read.csv("precip_monthly_anomalies_mean_ERA5.csv"),
                          year= seq(from=1981, to=2020+11/12, by=1/12 ) )
dat_era5_monthly$x = dat_era5_monthly$x * 15
####### ERA5 data
load(file = "D:/South_China/SPI/ERA_South_China_precipitation_anomalies.RData")
anomalies = region_anomalies_precip_spring  + region_anomalies_precip_summer +
  region_anomalies_precip_fall + region_anomalies_precip_winter
dat_era5 = data.frame (year= 1981:2020, anomalies = anomalies)

library(ggplot2)
year=seq(from=1960, to=2013+11/12, by=1/12 )

p = ggplot(dat) +
  geom_ribbon(  mapping = aes(x=year, ymax = ymax, ymin = ymin, fill=foring), alpha=0.3, size=0.1 ) +
  geom_line( mapping = aes(x=year, y = ymean, col=foring), size=0.2 ) +
  geom_line( data = dat_station_monthly, mapping = aes(x=year, y = x), col="blue", size=0.2, show.legend = T ) +
  geom_line( data = dat_era5_monthly, mapping = aes(x=year, y = x), col=1, size=0.2, show.legend = T ) +
  theme(panel.background = element_rect(fill="white",color="black",size=0.3),
        panel.grid = element_blank(),legend.position='bottom',
        axis.title=element_text(size=10),axis.text=element_text(size=8),
        axis.ticks=element_line(size=0.1))+
  #stat_ecdf(alpha=0.55, aes(fwi)) + 
  geom_hline(aes(yintercept=0), colour="black",size=0.2) +
  scale_x_continuous(breaks=seq(1950,2020,10)) +
  scale_y_continuous(breaks=seq(-60,600,10))
ggsave(filename = "forcings_precipitation_monthly_1.pdf",p,width = 8,height = 5)
  # geom_path(data = china, aes(x=long, y=lat, group = group), color = 1, show.legend = F) +
  # geom_point(data = data_plot, aes(x = lon, y = lat, color = spi), size = 2, alpha = 1) +
  # geom_text_repel(data = mat.cities, aes(x = long, y = lat, label = names), family = "STHeiti") +
  # geom_text(data = mat.cities, aes(x = long, y = lat, label = names), family = "STHeiti") +
  # labs(x = 'Longitude', y = 'Latitude', title = 'SPI_12', size = 10 ) +
  # ylim(18,55) +
  # theme_bw() +
  # scale_color_gradient2(low = "red", mid = "white", high = "blue", midpoint = .0,
  #                       guide = guide_colourbar(direction = "horizontal",title.position = "bottom" ) ) + 
  # theme(text = element_text(family = "STHeiti"),
  #       #        panel.border = element_blank(),
  #       plot.title = element_text(hjust = 0.5),
  #       legend.position = "bottom" )

############## annual data figure ##############

dat_annual = dat[1:(1/12*nrow(dat)), ]
for(i in 1:nrow(dat_annual) ){
  dat_annual[i,] = c(dat[(i-1)*12+1, 1], colSums( dat[((i-1)*12+1):(i*12), 2:4 ]), dat[(i-1)*12+1, 5] )
}
dat_annual$year = as.numeric(dat_annual$year)
dat_annual$ymean = as.numeric(dat_annual$ymean)
dat_annual$ymax = as.numeric(dat_annual$ymax)
dat_annual$ymin = as.numeric(dat_annual$ymin)

p = ggplot(dat_annual) +
  geom_ribbon(  mapping = aes(x=year, ymax = ymax, ymin = ymin, fill=foring), alpha=0.3, size=0.1 ) +
  geom_line( mapping = aes(x=year, y = ymean, col=foring), size=0.5 ) +
  geom_line( data = dat_station, mapping = aes(x=year, y = sc), col="blue", size=0.5, show.legend = T ) +
  geom_line( data = dat_era5, mapping = aes(x=year, y = anomalies), col=1, size=0.5, show.legend = T ) +
  theme(panel.background = element_rect(fill="white",color="black",size=0.3),
        panel.grid = element_blank(),legend.position='bottom',
        axis.title=element_text(size=10),axis.text=element_text(size=8),
        axis.ticks=element_line(size=0.1))+
  #stat_ecdf(alpha=0.55, aes(fwi)) + 
  geom_hline(aes(yintercept=0), colour="black",size=0.2) +
  scale_x_continuous(breaks=seq(1950,2020,10)) +
  scale_y_continuous(breaks=seq(-800,800,200))
ggsave(filename = "forcings_precipitation_annual_1.pdf",p,width = 8,height = 5)


######################################################################################################
######################################  SPI frequency analyses #######################################
######################################################################################################

# module load  nixpkgs/16.09
# module load gcc/7.3.0
# module load netcdf/4.6.1
# module load r/3.6.1
# module load gdal/3.0.1
# module load proj/6.3.0
# module load udunits/2.2.26

library(ncdf4)

spi = wetrun = dryrun = array(NA, dim = c(15,14,16,648))
for(i_model in 1:15){
  file = nc_open(filename = paste0("Hist_spi_12_r1i1p",i_model,".nc"))
  time <- ncvar_get(file, varid = "time")
  lon <- ncvar_get(file, varid = "lon")
  lat <- ncvar_get(file, varid = "lat")
  
#  spi[lon,lat,height,time]
##  stduy region (part of south China) lat 20:28 (197:212) , lon, 110-120 (131:144)
  
  # precip_south_china = dat(lat|197:212, lon|131:144,  time|:, height|0)    ; extract the study region data
  
  spi[i_model,,,] = ncvar_get(file, varid="spi", start = c(132,198,1,1), count = c(14,16,-1,-1),)
  
  file = nc_open(filename = paste0("Hist_wetrun_12_r1i1p",i_model,".nc"))
  wetrun[i_model,,,] = ncvar_get(file, varid="wetrun", start = c(132,198,1,1), count = c(14,16,-1,-1),)
  
  file = nc_open(filename = paste0("Hist_dryrun_12_r1i1p",i_model,".nc"))
  dryrun[i_model,,,] = ncvar_get(file, varid="dryrun", start = c(132,198,1,1), count = c(14,16,-1,-1),)
}

save(spi, wetrun, dryrun, file = "Hist_regimes.RData")

spi = wetrun = dryrun = array(NA, dim = c(15,14,16,648))
for(i_model in 1:15){
  file = nc_open(filename = paste0("/scratch/xtan/Nat-Hist/Nat-Hist_spi_12_r1i1p",i_model,".nc"))
  time <- ncvar_get(file, varid = "time")
  lon <- ncvar_get(file, varid = "lon")
  lat <- ncvar_get(file, varid = "lat")

  #  spi[lon,lat,height,time]
  ##  stduy region (part of south China) lat 20:28 (197:212) , lon, 110-120 (131:144)

  # precip_south_china = dat(lat|197:212, lon|131:144,  time|:, height|0)    ; extract the study region data

  spi[i_model,,,] = ncvar_get(file, varid="spi", start = c(132,198,1,1), count = c(14,16,-1,-1),)

  file = nc_open(filename = paste0("/scratch/xtan/Nat-Hist/Nat-Hist_wetrun_12_r1i1p",i_model,".nc"))
  wetrun[i_model,,,] = ncvar_get(file, varid="wetrun", start = c(132,198,1,1), count = c(14,16,-1,-1),)

  file = nc_open(filename = paste0("/scratch/xtan/Nat-Hist/Nat-Hist_dryrun_12_r1i1p",i_model,".nc"))
  dryrun[i_model,,,] = ncvar_get(file, varid="dryrun", start = c(132,198,1,1), count = c(14,16,-1,-1),)
}

save(spi, wetrun, dryrun, file = "NAT_regimes.RData")

#######################################################################################################
load(file = "Hist_regimes.RData")
dat = NULL
for(i_model in 1:15){
  aaa = as.vector(spi[i_model,,,])
  dat = rbind (dat, cbind ( aaa, rep(i_model, length(aaa)), rep("HIST", length(aaa) ) ) )  
} 
load(file = "NAT_regimes.RData")
for(i_model in 1:15){
  aaa = as.vector(spi[i_model,,,])
  dat = rbind (dat, cbind ( aaa, rep(i_model, length(aaa)), rep("NAT", length(aaa) ) ) )  
}
colnames(dat) = c("SPI", "model", "forcing")
dat = as.data.frame(dat)
dat$SPI = as.numeric(dat$SPI)

library(ggplot2)
p = ggplot(dat) +
  geom_density( mapping = aes(x= SPI, col= forcing, group=model ), size=0.1 ) +
  geom_density ( mapping = aes(x= SPI, col= forcing), size=0.5 ) +
  xlim (-3, 3) +
  # geom_hline(aes(yintercept=0), colour="black",size=0.2) + 
  theme(panel.background = element_rect(fill="white",color="black",size=0.3),
        panel.grid = element_blank(),legend.position='bottom',
        axis.title=element_text(size=10),axis.text=element_text(size=8),
        axis.ticks=element_line(size=0.1)) 
ggsave(filename = "SPI_density_C20.pdf",p,width = 6,height = 4)


load(file = "Hist_regimes.RData")
spi_hist = spi
load(file = "NAT_regimes.RData")
spi_nat = spi

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

rr_10 = rr(spi_hist, spi_nat, threshold = -1.0)
rr_15 = rr(spi_hist, spi_nat, threshold = -1.5)
rr_20 = rr(spi_hist, spi_nat, threshold = -2.0)
rr_25 = rr(spi_hist, spi_nat, threshold = -2.5)
rr_30 = rr(spi_hist, spi_nat, threshold = -3.0)

rr_wet_10 = rr(spi_hist, spi_nat, threshold = 1.0, wet = T)
rr_wet_15 = rr(spi_hist, spi_nat, threshold = 1.5, wet = T)
rr_wet_20 = rr(spi_hist, spi_nat, threshold = 2.0, wet = T)
rr_wet_25 = rr(spi_hist, spi_nat, threshold = 2.5, wet = T)
rr_wet_30 = rr(spi_hist, spi_nat, threshold = 3.0, wet = T)

##########################################################################################
############################## SPI frequency analyses ####################################
##########################################################################################

load(file = "Hist_regimes.RData")
spi_mean = apply(spi, FUN = mean, MARGIN = c(1, 4) )
spi_max = apply(spi, FUN = max, MARGIN = c( 4), na.rm=T )
spi_min = apply(spi, FUN = min, MARGIN = c( 4), na.rm=T )
spi_regional = apply(spi, FUN = mean, MARGIN = c( 4), na.rm=T )
hist(spi)

dat = NULL
dat = rbind(dat, cbind( year=seq(from=1960, to=2013+11/12, by=1/12 ), 
                         ymean = spi_regional , 
                         ymin = spi_min ,  
                         ymax = spi_max ,  
                         foring = rep("Hist", 648 ) ) ) 

load(file = "NAT_regimes.RData")
spi_mean = apply(spi, FUN = mean, MARGIN = c(1, 4) )
spi_max = apply(spi, FUN = max, MARGIN = c( 4), na.rm=T )
spi_min = apply(spi, FUN = min, MARGIN = c( 4), na.rm=T )
spi_regional = apply(spi, FUN = mean, MARGIN = c( 4), na.rm=T )
dat = rbind(dat, cbind( year=seq(from=1960, to=2013+11/12, by=1/12 ), 
                        ymean = spi_regional , 
                        ymin = spi_min ,  
                        ymax = spi_max ,  
                        foring = rep("NAT", 648 ) ) ) 

dat = as.data.frame(dat)
dat$year = as.numeric(dat$year)
dat$ymean = as.numeric(dat$ymean)
dat$ymax = as.numeric(dat$ymax)
dat$ymin = as.numeric(dat$ymin)

p = ggplot(dat) +
  geom_ribbon(  mapping = aes(x=year, ymax = ymax, ymin = ymin, fill=foring), alpha=0.3, size=0.1 ) +
  geom_line( mapping = aes(x=year, y = ymean, col=foring), size=0.5 ) +
  theme(panel.background = element_rect(fill="white",color="black",size=0.3),
        panel.grid = element_blank(),legend.position='bottom',
        axis.title=element_text(size=10),axis.text=element_text(size=8),
        axis.ticks=element_line(size=0.1))+
  #stat_ecdf(alpha=0.55, aes(fwi)) + 
  geom_hline(aes(yintercept=0), colour="black",size=0.2) +
  scale_x_continuous(breaks= seq(1960,2020,10) ) +
  scale_y_continuous(breaks= seq(-3,3,0.5) )
ggsave(filename = "forcings_spi_12.pdf",p,width = 8,height = 5)


#############################################################################################
################## function for retrieval drought characteristics ###########################
#############################################################################################

#### new function
#### run_theory() 
#### applied the run theory to a time series ()
## time_serie: a numeric vector with no NA values
## threshold: a numeric value in which the features (below) of run theory is measured
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

library(magrittr)

load(file = "Hist_regimes.RData")

n = dim(spi)
wet = dry = NULL
for(i_model in 1:n[1]){
  for(i_lon in 1:n[2]){
    for(i_lat in 1:n[3]){
      
      dat = spi[i_model,i_lon, i_lat, 12:648]
      
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

save(wet, dry, file = "HIST_drought_regimes.RData")


load(file = "NAT_regimes.RData")

n = dim(spi)
wet = dry = NULL
for(i_model in 1:n[1]){
  for(i_lon in 1:n[2]){
    for(i_lat in 1:n[3]){
      
      dat = spi[i_model,i_lon, i_lat, 12:648]
      
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
save(wet, dry, file = "NAT_drought_regimes.RData")

###################################################################################################
###################################################################################################
###################################################################################################

load(file = "HIST_drought_regimes.RData")
dat = cbind(dry, forcing = rep( "HIST", nrow(dry) ) )
load(file = "NAT_drought_regimes.RData")
dat = rbind(dat, cbind ( dry , forcing = rep( "NAT", nrow(dry) ))  )

dat = data.frame(dat)
dat$Duration = as.numeric(dat$Duration) 
dat$Severity = as.numeric(dat$Severity)
dat$Intesity = as.numeric(dat$Intesity)
dat$Date_Ini_Ev = as.numeric(dat$Date_Ini_Ev)
dat$Date_Fin_Ev = as.numeric(dat$Date_Fin_Ev)

library(ggplot2)
p = ggplot(dat) +
  geom_density( mapping = aes(x= Duration, col= forcing, group=model ), size=0.1) +
  geom_density ( mapping = aes(x= Duration, col= forcing), size=0.5 ) +
  xlim (0, 30) +
  # geom_hline(aes(yintercept=0), colour="black",size=0.2) + 
  theme(panel.background = element_rect(fill="white",color="black",size=0.3),
        panel.grid = element_blank(),legend.position='bottom',
        axis.title=element_text(size=10),axis.text=element_text(size=8),
        axis.ticks=element_line(size=0.1)) 
ggsave(filename = "duration_density_C20.pdf",p,width = 1.5,height = 3)

library(ggplot2)
p = ggplot(dat) +
  geom_density( mapping = aes(x= Severity, col= forcing, group=model ), size=0.1) +
  geom_density ( mapping = aes(x= Severity, col= forcing), size=0.5 ) +
  xlim (0, 50) +
  # geom_hline(aes(yintercept=0), colour="black",size=0.2) + 
  theme(panel.background = element_rect(fill="white",color="black",size=0.3),
        panel.grid = element_blank(),legend.position='bottom',
        axis.title=element_text(size=10),axis.text=element_text(size=8),
        axis.ticks=element_line(size=0.1)) 
ggsave(filename = "Severity_density_C20.pdf",p,width = 1.5,height = 3)

library(ggplot2)
p = ggplot(dat) +
  geom_density( mapping = aes(x= Intesity, col= forcing, group=model ), size=0.1) +
  geom_density ( mapping = aes(x= Intesity, col= forcing), size=0.5 ) +
  xlim (0, 3) +
  # geom_hline(aes(yintercept=0), colour="black",size=0.2) + 
  theme(panel.background = element_rect(fill="white",color="black",size=0.3),
        panel.grid = element_blank(),legend.position='bottom',
        axis.title=element_text(size=10),axis.text=element_text(size=8),
        axis.ticks=element_line(size=0.1)) 
ggsave(filename = "Intensity_density_C20.pdf",p,width = 2,height = 3)





