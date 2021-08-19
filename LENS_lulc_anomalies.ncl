

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
  clim_mean = asciiread("climatology.txt",(/10,9,12/),"float")
  precip_south_china = dat1(lat|117:126, lon|88:96,  time|:)
  tp_anomalies  = calcMonAnomLLT(precip_south_china, clim_mean  ) ;Calculate Anomalies from Means
  asciiwrite("precipitation_anomalies_LULC_"+models(model)+".txt",tp_anomalies)


	delete(dat1)

end do;


end;

