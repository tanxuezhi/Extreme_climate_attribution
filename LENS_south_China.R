
rm(list=ls())

setwd("F:/data/LENS_single_forcing")

forcings = c("AER","BMB","GHG","hist","LULC")

models = list( c("001","002","003","004","005","007","008","009","010","011","012","013","014","015","016","017","018","019","020"),
               c("001","002","003","004","005","006","007","008","009","010","011","012","013","014","015"),
               c("001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020"), 
               c("002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020",
                 "021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","101","102","103","104","105"),
               c("001","002","003","004","005")      )

climatology = array( read.table("climatology.txt")[[1]],  dim=c(12,9,10) ) 
climatology = climatology * 86400 * 1000 * 30

for(i_forcing in 1){  ####### AER forcing 
  anomalies = array(NA, dim =c( length(models[[i_forcing]]),1932,9,10) )
  for(i_model in 1:length(models[[i_forcing]])){
    ####### data dimension lat 10 lon 9 time 1932
    dat = read.table(file = paste0("precipitation_anomalies_",forcings[i_forcing],"_",models[[i_forcing]][i_model], ".txt" ) )
    anomalies [i_model, , ,] = array(data = dat[[1]], dim=c(1932,9,10)) * 86400 * 1000 * 30 
  }
  anomalies_AER = anomalies
}


for(i_forcing in 2){  ####### BMB forcing 
  anomalies = array(NA, dim =c( length(models[[i_forcing]]), 1320,9,10) )
  for(i_model in 1:length(models[[i_forcing]])){
    ####### data dimension lat 10 lon 9 time 1932
    dat = read.table(file = paste0("precipitation_anomalies_",forcings[i_forcing],"_",models[[i_forcing]][i_model], ".txt" ) )
    anomalies [i_model, , ,] = array(data = dat[[1]], dim=c(1320,9,10)) * 86400 * 1000 * 30 
  }
  anomalies_BMB = anomalies
}

for(i_forcing in 3){  ####### GHG forcing 
  anomalies = array(NA, dim =c( length(models[[i_forcing]]), 1932,9,10) )
  for(i_model in 1:length(models[[i_forcing]])){
    ####### data dimension lat 10 lon 9 time 1932
    dat = read.table(file = paste0("precipitation_anomalies_",forcings[i_forcing],"_",models[[i_forcing]][i_model], ".txt" ) )
    anomalies [i_model, , ,] = array(data = dat[[1]], dim=c(1932,9,10)) * 86400 * 1000 * 30 
  }
  anomalies_GHG = anomalies
}

for(i_forcing in 4){  ####### HIST forcing 
  anomalies = array(NA, dim =c( length(models[[i_forcing]]), 2172,9,10) )
  for(i_model in 1:length(models[[i_forcing]])){
    ####### data dimension lat 10 lon 9 time 1932
    dat = read.table(file = paste0("precipitation_anomalies_",forcings[i_forcing],"_",models[[i_forcing]][i_model], ".txt" ) )
    anomalies [i_model, , ,] = array(data = dat[[1]], dim=c(2172,9,10)) * 86400 * 1000 * 30 
  }
  anomalies_HIST = anomalies
}

for(i_forcing in 5){  ####### LULC forcing 
  anomalies = array(NA, dim =c( length(models[[i_forcing]]), 1320,9,10) )
  for(i_model in 1:length(models[[i_forcing]])){
    ####### data dimension lat 10 lon 9 time 1932
    dat = read.table(file = paste0("precipitation_anomalies_",forcings[i_forcing],"_",models[[i_forcing]][i_model], ".txt" ) )
    anomalies [i_model, , ,] = array(data = dat[[1]], dim=c(1320,9,10)) * 86400 * 1000 * 30 
  }
  anomalies_LULC = anomalies
}

library(RcppRoll)

filter = 12*10

anomalies_AER_mean = apply(anomalies_AER, FUN = mean, MARGIN = c(1, 2) )
anomalies_AER_mean_roll = apply(anomalies_AER_mean, FUN = roll_mean, MARGIN = c(1), n=filter, fill=NA )
anomalies_AER_mean_quantile = apply(anomalies_AER_mean_roll, FUN = quantile, MARGIN = c(1), 
                                    probs=c(0.1,0.5,0.9), na.rm=T )
plot(anomalies_AER_mean_quantile[2,])

anomalies_BMB_mean = apply(anomalies_BMB, FUN = mean, MARGIN = c(1, 2) )
anomalies_BMB_mean_roll = apply(anomalies_BMB_mean, FUN = roll_mean, MARGIN = c(1), n=filter, fill=NA )
anomalies_BMB_mean_quantile = apply(anomalies_BMB_mean_roll, FUN = quantile, MARGIN = c(1), 
                                    probs=c(0.1,0.5,0.9), na.rm=T )
plot(anomalies_BMB_mean_quantile[2,])

anomalies_GHG_mean = apply(anomalies_GHG, FUN = mean, MARGIN = c(1, 2) )
anomalies_GHG_mean_roll = apply(anomalies_GHG_mean, FUN = roll_mean, MARGIN = c(1), n=filter, fill=NA )
anomalies_GHG_mean_quantile = apply(anomalies_GHG_mean_roll, FUN = quantile, MARGIN = c(1), 
                                    probs=c(0.1,0.5,0.9), na.rm=T )
plot(anomalies_GHG_mean_quantile[2,])

anomalies_HIST_mean = apply(anomalies_HIST, FUN = mean, MARGIN = c(1, 2) )
anomalies_HIST_mean_roll = apply(anomalies_HIST_mean, FUN = roll_mean, MARGIN = c(1), n=filter, fill=NA )
anomalies_HIST_mean_quantile = apply(anomalies_HIST_mean_roll, FUN = quantile, MARGIN = c(1), 
                                    probs=c(0.1,0.5,0.9), na.rm=T )
plot(anomalies_HIST_mean_quantile[2,])

anomalies_LULC_mean = apply(anomalies_LULC, FUN = mean, MARGIN = c(1, 2) )
anomalies_LULC_mean_roll = apply(anomalies_LULC_mean, FUN = roll_mean, MARGIN = c(1), n=filter, fill=NA )
anomalies_LULC_mean_quantile = apply(anomalies_LULC_mean_roll, FUN = quantile, MARGIN = c(1), 
                                    probs=c(0.1,0.5,0.9), na.rm=T )
plot(anomalies_LULC_mean_quantile[2,])

climatology_mean = apply(climatology, FUN = mean, MARGIN = c(1) )

library(abind)
anomalies_mean_quantile = rbind(anomalies_AER_mean_quantile, anomalies_GHG_mean_quantile)
forcings = c("AER","GHG","BMB","LULC", "HIST")
dat = NULL
for(i_forcing in 1:2){  ####### HIST forcing 
    dat = rbind( dat, cbind( year=seq(from=1920, to=2080+11/12, by=1/12 ), 
                             ymean = as.numeric( anomalies_mean_quantile[(i_forcing-1)*3+2,]) , 
                             ymin = as.numeric( anomalies_mean_quantile[(i_forcing-1)*3+1,]) ,  
                             ymax = as.numeric( anomalies_mean_quantile[(i_forcing-1)*3+3,]) ,  
                             foring = rep(forcings[i_forcing], 1932 ) ) )
}
anomalies_mean_quantile = rbind(anomalies_BMB_mean_quantile, anomalies_LULC_mean_quantile)
for(i_forcing in 3:4){  ####### HIST forcing 
  dat = rbind( dat, cbind( year=seq(from=1920, to=2029+11/12, by=1/12 ), 
                           ymean = as.numeric( anomalies_mean_quantile[(i_forcing-3)*3+2,]) , 
                           ymin = as.numeric( anomalies_mean_quantile[(i_forcing-3)*3+1,]) ,  
                           ymax = as.numeric( anomalies_mean_quantile[(i_forcing-3)*3+3,]) ,  
                           foring = rep(forcings[i_forcing], 1320 ) ) )
}
dat = rbind(dat, cbind( year=seq(from=1920, to=2100+11/12, by=1/12 ), 
                        ymean = as.numeric( anomalies_HIST_mean_quantile[2,]) , 
                        ymin = as.numeric( anomalies_HIST_mean_quantile[1,] ) ,  
                        ymax = as.numeric( anomalies_HIST_mean_quantile[3,] ) ,  
                        foring = rep("HIST", 2172 ) ) )
dat = as.data.frame(dat)
dat$year = as.numeric(dat$year)
dat$ymean = as.numeric(dat$ymean)
dat$ymax = as.numeric(dat$ymax)
dat$ymin = as.numeric(dat$ymin)

library(ggplot2)

p = ggplot(dat) +
  geom_ribbon(  mapping = aes(x=year, ymax = ymax, ymin = ymin, fill=foring), alpha=0.3, size=0.1 ) +
  geom_line( mapping = aes(x=year, y = ymean, col=foring), size=0.5 ) +
  theme(panel.background = element_rect(fill="white",color="black",size=0.3),
        panel.grid = element_blank(),legend.position='bottom',
        axis.title=element_text(size=10),axis.text=element_text(size=8),
        axis.ticks=element_line(size=0.1))+
  #stat_ecdf(alpha=0.55, aes(fwi)) + 
  geom_hline(aes(yintercept=0), colour="black",size=0.2) +
  scale_x_continuous(breaks=c(1925,1950, 1975, 2000, 2025, 2050, 2075,2100)) +
  scale_y_continuous(breaks=seq(-20,30,5))
ggsave(filename = "forcings_precipitation_monthly.pdf",p,width = 5,height = 5)
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
  theme(panel.background = element_rect(fill="white",color="black",size=0.3),
        panel.grid = element_blank(),legend.position='bottom',
        axis.title=element_text(size=10),axis.text=element_text(size=8),
        axis.ticks=element_line(size=0.1))+
  #stat_ecdf(alpha=0.55, aes(fwi)) + 
  geom_hline(aes(yintercept=0), colour="black",size=0.2) +
  scale_x_continuous(breaks=c(1925,1950, 1975, 2000, 2025, 2050, 2075,2100)) +
  scale_y_continuous(breaks=seq(-200,300,50))
ggsave(filename = "forcings_precipitation_annual.pdf",p,width = 8,height = 5)






