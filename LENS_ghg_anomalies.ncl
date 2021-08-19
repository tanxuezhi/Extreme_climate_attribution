

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
  clim_mean = asciiread("climatology.txt",(/10,9,12/),"float")
  precip_south_china = dat1(lat|117:126, lon|88:96,  time|:)
  tp_anomalies  = calcMonAnomLLT(precip_south_china, clim_mean  ) ;Calculate Anomalies from Means
  asciiwrite("precipitation_anomalies_GHG_"+models(model)+".txt",tp_anomalies)

delete(dat1)

end do;

end;

