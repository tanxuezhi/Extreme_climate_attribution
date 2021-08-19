

begin

;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

; uq_his=hus_his*ua_his

; uq_rcp=hus_rcp*ua_rcp

; vq_his=hus_his*va_his

; vq_rcp=hus_rcp*va_rcp

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;; C20;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

month = (/3,4,5,6,9,12/)

; scenario historical

	models = (/"002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","101","102","103","104","105"/)
	years = (/"192001-200512"/)

 do model=0, dimsizes(models)-1
	file_precip_0 = addfile("/scratch/xtan/LENS_single_forcing/PRECL/b.e11.B20TRC5CNBDRD.f09_g16."+models(model)+".cam.h0.PRECL."+years(0)+".nc","r")
	dat = file_precip_0->PRECL
	; print(d)
	file_precip_0 = addfile("/scratch/xtan/LENS_single_forcing/PRECC/b.e11.B20TRC5CNBDRD.f09_g16."+models(model)+".cam.h0.PRECC."+years(0)+".nc","r")
	dat0 = file_precip_0->PRECC

	lat := file_precip_0->lat
	lon := file_precip_0->lon

	dat1 = dat + dat0  
	copy_VarMeta(dat,dat1)

        delete(dat)
	delete(dat0)

;; stduy region (part of south China) lat 20:28 (117:126) , lon, 110-120 (88:96)

  precip_south_china = dat1(lat|117:126, lon|88:96,  time|:)    ; extract the study region data

  clim = new( (/39,10,9,12/),float)

  clim(model,:,:,:) = clmMonLLT (precip_south_china(lat|:, lon|:, time|492:851) ) ; calculate the study region precipitation climatology

  tp_anomalies  = calcMonAnomLLT(precip_south_china, clim(model,:,:,:))  ;Calculate Anomalies from Means

  asciiwrite("precipitation_anomalies_hist_"+models(model)+".txt",tp_anomlies)

;; stduy region (part of south China) lat 20:28 (117:126) , lon, 110-120 (88:96)

	do i_month= 0,dimsizes(month)-1

	spi = dim_spi_n (dat1, month(i_month), False, 0)
	copy_VarMeta(dat1,spi)
	ncdf = addfile("BMB_spi_"+month(i_month)+"_"+models(model)+".nc" ,"c")  ; open output netCDF file
	ncdf->spi = spi

	wet = where(spi.gt.1.0, 1, 0)
	copy_VarMeta(dat1,wet)
	ncdf = addfile("BMB_wet_"+month(i_month)+"_"+models(model)+".nc" ,"c")  ; open output netCDF file
	ncdf->wet = wet

	dry = where(spi.lt.-1.0, 1,0)
	copy_VarMeta(dat1,dry)
	ncdf = addfile("BMB_dry_"+month(i_month)+"_"+models(model)+".nc" ,"c")  ; open output netCDF file
	ncdf->dry = dry

	dryrun = dim_numrun_n (dry, 1, 0) 
	copy_VarMeta(dat1,dryrun)
	ncdf = addfile("BMB_dryrun_"+month(i_month)+"_"+models(model)+".nc" ,"c")  ; open output netCDF file
	ncdf->dryrun = dryrun

	wetrun = dim_numrun_n (wet, 1, 0) 
	copy_VarMeta(dat1,wetrun)
	ncdf = addfile("BMB_wetrun_"+month(i_month)+"_"+models(model)+".nc" ,"c")  ; open output netCDF file
	ncdf->wetrun = wetrun

	delete(spi)
	delete(wet)
	delete(dry)
	delete(wetrun)
	delete(dryrun)

        end do;
	delete(dat1)	
	
end do;

  clim_mean = dim_avg_n (clim, 0)
  asciiwrite("climatology.txt",clim_mean)

; scenario rcp85 
	models = (/"002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","101","102","103","104","105"/)
	
 do model=0, dimsizes(models)-1
 	
	if model .gt. 32 then
	years := (/"200601-210012"/)
	file_precip_0 = addfile("/scratch/xtan/LENS_single_forcing/PRECL/b.e11.B20TRC5CNBDRD.f09_g16."+models(model)+".cam.h0.PRECL."+years(0)+".nc","r")
	dat = file_precip_0->PRECL
	; print(d)
	file_precip_0 = addfile("/scratch/xtan/LENS_single_forcing/PRECC/b.e11.B20TRC5CNBDRD.f09_g16."+models(model)+".cam.h0.PRECC."+years(0)+".nc","r")
	dat0 = file_precip_0->PRECC

	lat := file_precip_0->lat
	lon := file_precip_0->lon

	dat1 = dat + dat0  
	copy_VarMeta(dat,dat1)

        delete(dat)
	delete(dat0)

	
	else
	years := (/"200601-208012","208101-210012"/)
	file_precip_0 = addfile("/scratch/xtan/LENS_single_forcing/PRECL/b.e11.B20TRC5CNBDRD.f09_g16."+models(model)+".cam.h0.PRECL."+years(0)+".nc","r")
	dat1 = file_precip_0->PRECL
	time1 = file_precip_0->time
	file_precip_1 = addfile("/scratch/xtan/LENS_single_forcing/PRECL/b.e11.B20TRC5CNBDRD.f09_g16."+models(model)+".cam.h0.PRECL."+years(1)+".nc","r")
	dat2 = file_precip_1->PRECL
	dat = array_append_record (dat1, dat2, 0)
	delete(dat1)
	delete(dat2)
	d = getvardims(file_precip_0)
	; print(d)
	file_precip_0 = addfile("/scratch/xtan/LENS_single_forcing/PRECC/b.e11.B20TRC5CNBDRD.f09_g16."+models(model)+".cam.h0.PRECC."+years(0)+".nc","r")
	dat1 = file_precip_0->PRECC
	file_precip_1 = addfile("/scratch/xtan/LENS_single_forcing/PRECC/b.e11.B20TRC5CNBDRD.f09_g16."+models(model)+".cam.h0.PRECC."+years(1)+".nc","r")
	dat2 = file_precip_1->PRECC
	time2= file_precip_1->time
	dat0 = array_append_record (dat1, dat2, 0)
	delete(dat1)
	delete(dat2)
	lat := file_precip_0->lat
	lon := file_precip_0->lon

	dat1 = dat + dat0  
        delete(dat)
	delete(dat0)

	dat1!1 = "lat"
	dat1!2 = "lon"
	dat1!0 = "time"
	time := array_append_record (time1, time2, 0)
	dat1@time = time
	dat1@lat = lat
	dat1@lon = lon

	end


;; stduy region (part of south China) lat 20:28 (117:126) , lon, 110-120 (88:96)

  precip_south_china = dat1(lat|117:126, lon|88:96,  time|:)
  tp_anomalies  = calcMonAnomLLT(precip_south_china, clim(model,:,:,:)  ) ;Calculate Anomalies from Means
  asciiwrite("precipitation_anomalies_rcp85_"+models(model)+".txt",tp_anomlies)

	do i_month= 0,dimsizes(month)-1

	spi = dim_spi_n (dat1, month(i_month), False, 0)
	copy_VarMeta(dat1,spi)
	ncdf = addfile("GHG_spi_"+month(i_month)+"_"+models(model)+".nc" ,"c")  ; open output netCDF file
	ncdf->spi = spi

	wet = where(spi.gt.1.0, 1, 0)
	copy_VarMeta(dat1,wet)
	ncdf = addfile("GHG_wet_"+month(i_month)+"_"+models(model)+".nc" ,"c")  ; open output netCDF file
	ncdf->wet = wet

	dry = where(spi.lt.-1.0, 1,0)
	copy_VarMeta(dat1,dry)
	ncdf = addfile("GHG_dry_"+month(i_month)+"_"+models(model)+".nc" ,"c")  ; open output netCDF file
	ncdf->dry = dry

	dryrun = dim_numrun_n (dry, 1, 0) 
	copy_VarMeta(dat1,dryrun)
	ncdf = addfile("GHG_dryrun_"+month(i_month)+"_"+models(model)+".nc" ,"c")  ; open output netCDF file
	ncdf->dryrun = dryrun

	wetrun = dim_numrun_n (wet, 1, 0) 
	copy_VarMeta(dat1,wetrun)
	ncdf = addfile("GHG_wetrun_"+month(i_month)+"_"+models(model)+".nc" ,"c")  ; open output netCDF file
	ncdf->wetrun = wetrun

	delete(spi)
	delete(wet)
	delete(dry)
	delete(wetrun)
	delete(dryrun)

end do;

	delete(dat1)

end do;


; scenario xghg 

	models = (/"001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020"/)
	years = (/"192001-200512","200601-208012"/)


 do model=0, dimsizes(models)-1
	file_precip_0 = addfile("/scratch/xtan/LENS_single_forcing/PRECL/b.e11.B20TRLENS_RCP85.f09_g16.xghg."+models(model)+".cam.h0.PRECL."+years(0)+".nc","r")
	dat1 = file_precip_0->PRECL
	time1 = file_precip_0->time
	file_precip_1 = addfile("/scratch/xtan/LENS_single_forcing/PRECL/b.e11.B20TRLENS_RCP85.f09_g16.xghg."+models(model)+".cam.h0.PRECL."+years(1)+".nc","r")
	dat2 = file_precip_1->PRECL
	dat = array_append_record (dat1, dat2, 0)
	delete(dat1)
	delete(dat2)
	d = getvardims(file_precip_0)
	; print(d)
	file_precip_0 = addfile("/scratch/xtan/LENS_single_forcing/PRECC/b.e11.B20TRLENS_RCP85.f09_g16.xghg."+models(model)+".cam.h0.PRECC."+years(0)+".nc","r")
	dat1 = file_precip_0->PRECC
	file_precip_1 = addfile("/scratch/xtan/LENS_single_forcing/PRECC/b.e11.B20TRLENS_RCP85.f09_g16.xghg."+models(model)+".cam.h0.PRECC."+years(1)+".nc","r")
	dat2 = file_precip_1->PRECC
	time2= file_precip_1->time
	dat0 = array_append_record (dat1, dat2, 0)
	delete(dat1)
	delete(dat2)

	lat := file_precip_0->lat
	lon := file_precip_0->lon


	dat1 = dat + dat0  
        delete(dat)
	delete(dat0)

	dat1!1 = "lat"
	dat1!2 = "lon"
	dat1!0 = "time"
	time := array_append_record (time1, time2, 0)
	dat1@time = time
	dat1@lat = lat
	dat1@lon = lon

;; stduy region (part of south China) lat 20:28 (117:126) , lon, 110-120 (88:96)

  precip_south_china = dat1(lat|117:126, lon|88:96,  time|:)
  tp_anomalies  = calcMonAnomLLT(precip_south_china, clim_mean  ) ;Calculate Anomalies from Means
  asciiwrite("precipitation_anomalies_GHG_"+models(model)+".txt",tp_anomlies)

	do i_month= 0,dimsizes(month)-1

	spi = dim_spi_n (dat1, month(i_month), False, 0)
	copy_VarMeta(dat1,spi)
	ncdf = addfile("GHG_spi_"+month(i_month)+"_"+models(model)+".nc" ,"c")  ; open output netCDF file
	ncdf->spi = spi

	wet = where(spi.gt.1.0, 1, 0)
	copy_VarMeta(dat1,wet)
	ncdf = addfile("GHG_wet_"+month(i_month)+"_"+models(model)+".nc" ,"c")  ; open output netCDF file
	ncdf->wet = wet

	dry = where(spi.lt.-1.0, 1,0)
	copy_VarMeta(dat1,dry)
	ncdf = addfile("GHG_dry_"+month(i_month)+"_"+models(model)+".nc" ,"c")  ; open output netCDF file
	ncdf->dry = dry

	dryrun = dim_numrun_n (dry, 1, 0) 
	copy_VarMeta(dat1,dryrun)
	ncdf = addfile("GHG_dryrun_"+month(i_month)+"_"+models(model)+".nc" ,"c")  ; open output netCDF file
	ncdf->dryrun = dryrun

	wetrun = dim_numrun_n (wet, 1, 0) 
	copy_VarMeta(dat1,wetrun)
	ncdf = addfile("GHG_wetrun_"+month(i_month)+"_"+models(model)+".nc" ,"c")  ; open output netCDF file
	ncdf->wetrun = wetrun

	delete(spi)
	delete(wet)
	delete(dry)
	delete(wetrun)
	delete(dryrun)

end do;

	delete(dat1)

end do;

; scenario xaer 

	models = (/"001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020"/)
	years = (/"192001-200512","200601-208012"/)
 do model=0, dimsizes(models)-1
	file_precip_0 = addfile("/scratch/xtan/LENS_single_forcing/PRECL/b.e11.B20TRLENS_RCP85.f09_g16.xaer."+models(model)+".cam.h0.PRECL."+years(0)+".nc","r")
	dat1 = file_precip_0->PRECL
	time1 = file_precip_0->time
	file_precip_1 = addfile("/scratch/xtan/LENS_single_forcing/PRECL/b.e11.B20TRLENS_RCP85.f09_g16.xaer."+models(model)+".cam.h0.PRECL."+years(1)+".nc","r")
	dat2 = file_precip_1->PRECL
	dat = array_append_record (dat1, dat2, 0)
	delete(dat1)
	delete(dat2)
	d = getvardims(file_precip_0)
	; print(d)
	file_precip_0 = addfile("/scratch/xtan/LENS_single_forcing/PRECC/b.e11.B20TRLENS_RCP85.f09_g16.xaer."+models(model)+".cam.h0.PRECC."+years(0)+".nc","r")
	dat1 = file_precip_0->PRECC
	file_precip_1 = addfile("/scratch/xtan/LENS_single_forcing/PRECC/b.e11.B20TRLENS_RCP85.f09_g16.xaer."+models(model)+".cam.h0.PRECC."+years(1)+".nc","r")
	dat2 = file_precip_1->PRECC
	time2= file_precip_1->time
	dat0 = array_append_record (dat1, dat2, 0)
	delete(dat1)
	delete(dat2)

	lat := file_precip_0->lat
	lon := file_precip_0->lon

	dat1 = dat + dat0  
        delete(dat)
	delete(dat0)

	dat1!1 = "lat"
	dat1!2 = "lon"
	dat1!0 = "time"
	time := array_append_record (time1, time2, 0)
	dat1@time = time
	dat1@lat = lat
	dat1@lon = lon


;; stduy region (part of south China) lat 20:28 (117:126) , lon, 110-120 (88:96)

  precip_south_china = dat1(lat|117:126, lon|88:96,  time|:)
  tp_anomalies  = calcMonAnomLLT(precip_south_china, clim_mean  ) ;Calculate Anomalies from Means
  asciiwrite("precipitation_anomalies_AER_"+models(model)+".txt",tp_anomlies)

	do i_month= 0,dimsizes(month)-1

	spi = dim_spi_n (dat1, month(i_month), False, 0)
	copy_VarMeta(dat1,spi)
	ncdf = addfile("AER_spi_"+month(i_month)+"_"+models(model)+".nc" ,"c")  ; open output netCDF file
	ncdf->spi = spi

	wet = where(spi.gt.1.0, 1, 0)
	copy_VarMeta(dat1,wet)
	ncdf = addfile("AER_wet_"+month(i_month)+"_"+models(model)+".nc" ,"c")  ; open output netCDF file
	ncdf->wet = wet

	dry = where(spi.lt.-1.0, 1,0)
	copy_VarMeta(dat1,dry)
	ncdf = addfile("AER_dry_"+month(i_month)+"_"+models(model)+".nc" ,"c")  ; open output netCDF file
	ncdf->dry = dry

	dryrun = dim_numrun_n (dry, 1, 0) 
	copy_VarMeta(dat1,dryrun)
	ncdf = addfile("AER_dryrun_"+month(i_month)+"_"+models(model)+".nc" ,"c")  ; open output netCDF file
	ncdf->dryrun = dryrun

	wetrun = dim_numrun_n (wet, 1, 0) 
	copy_VarMeta(dat1,wetrun)
	ncdf = addfile("AER_wetrun_"+month(i_month)+"_"+models(model)+".nc" ,"c")  ; open output netCDF file
	ncdf->wetrun = wetrun

	delete(spi)
	delete(wet)
	delete(dry)
	delete(wetrun)
	delete(dryrun)

end do;

	delete(dat1)

end do;

; scenario xbmb 

	models := (/"001","002","003","004","005","006","007","008","009","010","011","012","013","014","015"/)
	years = (/"192001-202912"/)

 do model=0, dimsizes(models)-1
	file_precip_0 = addfile("/scratch/xtan/LENS_single_forcing/PRECL/b.e11.B20TRLENS_RCP85.f09_g16.xbmb."+models(model)+".cam.h0.PRECL."+years(0)+".nc","r")
	dat = file_precip_0->PRECL
	; print(d)
	file_precip_0 = addfile("/scratch/xtan/LENS_single_forcing/PRECC/b.e11.B20TRLENS_RCP85.f09_g16.xbmb."+models(model)+".cam.h0.PRECC."+years(0)+".nc","r")
	dat0 = file_precip_0->PRECC

	lat := file_precip_0->lat
	lon := file_precip_0->lon

	dat1 = dat + dat0  
	copy_VarMeta(dat,dat1)

        delete(dat)
	delete(dat0)


;; stduy region (part of south China) lat 20:28 (117:126) , lon, 110-120 (88:96)
  precip_south_china = dat1(lat|117:126, lon|88:96,  time|:)
  tp_anomalies  = calcMonAnomLLT(precip_south_china, clim_mean  ) ;Calculate Anomalies from Means
  asciiwrite("precipitation_anomalies_BMB_"+models(model)+".txt",tp_anomlies)


	do i_month= 0,dimsizes(month)-1

	spi = dim_spi_n (dat1, month(i_month), False, 0)
	copy_VarMeta(dat1,spi)
	ncdf = addfile("BMB_spi_"+month(i_month)+"_"+models(model)+".nc" ,"c")  ; open output netCDF file
	ncdf->spi = spi

	wet = where(spi.gt.1.0, 1, 0)
	copy_VarMeta(dat1,wet)
	ncdf = addfile("BMB_wet_"+month(i_month)+"_"+models(model)+".nc" ,"c")  ; open output netCDF file
	ncdf->wet = wet

	dry = where(spi.lt.-1.0, 1,0)
	copy_VarMeta(dat1,dry)
	ncdf = addfile("BMB_dry_"+month(i_month)+"_"+models(model)+".nc" ,"c")  ; open output netCDF file
	ncdf->dry = dry

	dryrun = dim_numrun_n (dry, 1, 0) 
	copy_VarMeta(dat1,dryrun)
	ncdf = addfile("BMB_dryrun_"+month(i_month)+"_"+models(model)+".nc" ,"c")  ; open output netCDF file
	ncdf->dryrun = dryrun

	wetrun = dim_numrun_n (wet, 1, 0) 
	copy_VarMeta(dat1,wetrun)
	ncdf = addfile("BMB_wetrun_"+month(i_month)+"_"+models(model)+".nc" ,"c")  ; open output netCDF file
	ncdf->wetrun = wetrun

	delete(spi)
	delete(wet)
	delete(dry)
	delete(wetrun)
	delete(dryrun)

end do;

	delete(dat1)

end do;

; scenario xlulc 

	models = (/"001","002","003","004","005"/)
	years = (/"192001-202912"/)

 do model=0, dimsizes(models)-1
	file_precip_0 = addfile("/scratch/xtan/LENS_single_forcing/PRECL/b.e11.B20TRLENS_RCP85.f09_g16.xlulc."+models(model)+".cam.h0.PRECL."+years(0)+".nc","r")
	dat = file_precip_0->PRECL
	; print(d)
	file_precip_0 = addfile("/scratch/xtan/LENS_single_forcing/PRECC/b.e11.B20TRLENS_RCP85.f09_g16.xlulc."+models(model)+".cam.h0.PRECC."+years(0)+".nc","r")
	dat0 = file_precip_0->PRECC

	lat := file_precip_0->lat
	lon := file_precip_0->lon

	dat1 = dat + dat0  
	copy_VarMeta(dat,dat1)

        delete(dat)
	delete(dat0)

;; stduy region (part of south China) lat 20:28 (117:126) , lon, 110-120 (88:96)
  precip_south_china = dat1(lat|117:126, lon|88:96,  time|:)
  tp_anomalies  = calcMonAnomLLT(precip_south_china, clim_mean  ) ;Calculate Anomalies from Means
  asciiwrite("precipitation_anomalies_LULC_"+models(model)+".txt",tp_anomlies)


	do i_month= 0,dimsizes(month)-1

	spi = dim_spi_n (dat1, month(i_month), False, 0)
	copy_VarMeta(dat1,spi)
	ncdf = addfile("LULC_spi_"+month(i_month)+"_"+models(model)+".nc" ,"c")  ; open output netCDF file
	ncdf->spi = spi

	wet = where(spi.gt.1.0, 1, 0)
	copy_VarMeta(dat1,wet)
	ncdf = addfile("LULC_wet_"+month(i_month)+"_"+models(model)+".nc" ,"c")  ; open output netCDF file
	ncdf->wet = wet

	dry = where(spi.lt.-1.0, 1,0)
	copy_VarMeta(dat1,dry)
	ncdf = addfile("LULC_dry_"+month(i_month)+"_"+models(model)+".nc" ,"c")  ; open output netCDF file
	ncdf->dry = dry

	dryrun = dim_numrun_n (dry, 1, 0) 
	copy_VarMeta(dat1,dryrun)
	ncdf = addfile("LULC_dryrun_"+month(i_month)+"_"+models(model)+".nc" ,"c")  ; open output netCDF file
	ncdf->dryrun = dryrun

	wetrun = dim_numrun_n (wet, 1, 0) 
	copy_VarMeta(dat1,wetrun)
	ncdf = addfile("LULC_wetrun_"+month(i_month)+"_"+models(model)+".nc" ,"c")  ; open output netCDF file
	ncdf->wetrun = wetrun

	delete(spi)
	delete(wet)
	delete(dry)
	delete(wetrun)
	delete(dryrun)

end do;

	delete(dat1)

end do;


end;

