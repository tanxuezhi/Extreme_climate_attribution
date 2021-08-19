

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

years = (/"196001-196912","197001-197912","198001-198912","199001-199912","200001-200912","201001-201312"/)
month = (/3,4,5,6,9,12/)

 clim = new( (/15,16,14,12/),float)

 do model=1, 15
	file_precip_0 = addfile("/scratch/xtan/C20/mon/atmos/pr/r1i1p"+model+"/pr_Amon_HadGEM3-A-N216_All-Hist_est1_v1-0_r1i1p"+model+"_"+years(0)+".nc","r")
	dat1 = file_precip_0->pr
	file_precip_1 = addfile("/scratch/xtan/C20/mon/atmos/pr/r1i1p"+model+"/pr_Amon_HadGEM3-A-N216_All-Hist_est1_v1-0_r1i1p"+model+"_"+years(1)+".nc","r")
	dat2 = file_precip_1->pr
	dat = array_append_record (dat1, dat2, 0)
	delete(dat1)
	delete(dat2)
	file_precip_2 = addfile("/scratch/xtan/C20/mon/atmos/pr/r1i1p"+model+"/pr_Amon_HadGEM3-A-N216_All-Hist_est1_v1-0_r1i1p"+model+"_"+years(2)+".nc","r")
	dat2 = file_precip_2->pr
	dat1 = array_append_record (dat, dat2, 0)
	delete(dat)
	delete(dat2)
	file_precip_3 = addfile("/scratch/xtan/C20/mon/atmos/pr/r1i1p"+model+"/pr_Amon_HadGEM3-A-N216_All-Hist_est1_v1-0_r1i1p"+model+"_"+years(3)+".nc","r")
	dat2 = file_precip_3->pr
	dat = array_append_record (dat1, dat2, 0)
	delete(dat1)
	delete(dat2)
	file_precip_4 = addfile("/scratch/xtan/C20/mon/atmos/pr/r1i1p"+model+"/pr_Amon_HadGEM3-A-N216_All-Hist_est1_v1-0_r1i1p"+model+"_"+years(4)+".nc","r")
	dat2 = file_precip_4->pr
	dat1 = array_append_record (dat, dat2, 0)
	delete(dat)
	delete(dat2)
	file_precip_5 = addfile("/scratch/xtan/C20/mon/atmos/pr/r1i1p"+model+"/pr_Amon_HadGEM3-A-N216_All-Hist_est1_v1-0_r1i1p"+model+"_"+years(5)+".nc","r")
	dat2 = file_precip_5->pr
	dat = array_append_record (dat1, dat2, 0)
	delete(dat1)
	delete(dat2)
	d = getvardims(file_precip_0)
	; print(d)

	lat = file_precip_0->lat
	lon = file_precip_0->lon

;; stduy region (part of south China) lat 20:28 (197:212) , lon, 110-120 (131:144)
  
  precip_south_china = dat(lat|197:212, lon|131:144,  time|:, height|0)    ; extract the study region data

  clim(model-1,:,:,:) = clmMonLLT (precip_south_china(lat|:, lon|:, time|0:359) ) ; calculate the study region precipitation climatology

  tp_anomalies  = calcMonAnomLLT(precip_south_china, clim(model-1,:,:,:))  ;Calculate Anomalies from Means

  asciiwrite("precipitation_anomalies_hist_"+model+".txt",tp_anomalies)


	delete(dat)

end do;

  clim_mean = dim_avg_n (clim, 0)
  asciiwrite("climatology_c20.txt",clim_mean)


end;

