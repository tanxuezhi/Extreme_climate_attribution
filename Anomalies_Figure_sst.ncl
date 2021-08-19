

begin

;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

; uq_his=hus_his*ua_his

; uq_rcp=hus_rcp*ua_rcp

; vq_his=hus_his*va_his

; vq_rcp=hus_rcp*va_rcp

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;; Nat-Hist;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

 months = (/"Jan.","Feb.","Mar.","Apr.","May","Jun.","Jul.","Aug.","Sept.","Oct.","Nov.","Dec."/)

	file_precip_0 = addfile("total_precipitation.nc","r")
	lat = file_precip_0->latitude(400:900)
	lon = file_precip_0->longitude(600:1800)
	dat = file_precip_0->tp(time|:, latitude|400:900,longitude|600:1800)

;	file_precip_1 = addfile("total_precipitation_50-78.nc","r")
;	dat2 = file_precip_1->tp(time|:, latitude|400:900,longitude|600:1800)
;	dat1 = array_append_record (dat2, dat, 0)

	d = getvardims(file_precip_0)
	; print(d)

	file_precip_0 = addfile("SST.nc","r")
	sst =  file_precip_0->sst(time|0:503, expver|0, latitude|:,longitude|:)

	  dimice=dimsizes(dat)
	  ntim2=dimice(0)
	  nlat2=dimice(1)
 	  mlon2=dimice(2)

;====================================
; Calculate the JFM Clim.
;====================================
  temp   :=  dat(latitude|:,longitude|:,time|:)  
	temp!0 = "lat"
	temp!1 = "lon"
	temp!2 = "time"

;=====================================================================
; Calculate and Remove the Long-Term  Monthly Means, + compute JFM Avg.
;=====================================================================

  tp_clm = clmMonLLT (short2flt(temp))             ;Calculate LTMM's
  tp_anomalies  = calcMonAnomLLT(short2flt(temp),tp_clm) ;Calculate Anomalies from Means
  tp_anomalies = tp_anomalies * 1000 * 30
;  region_mean = dim_avg (tp_clm)

  temp   :=  sst(latitude|:,longitude|:,time|:)  
	temp!0 = "lat"
	temp!1 = "lon"
	temp!2 = "time"

  sst_clm = clmMonLLT (short2flt(temp))             ;Calculate LTMM's
  sst_anomalies  = calcMonAnomLLT(short2flt(temp),sst_clm) ;Calculate Anomalies from Means

;  region_mean = dim_avg (tp_clm)


;;  study region of south china  lat: 20:28 (220:300) & lon: 110:121(500:610)
  precip_south_china = tp_anomalies(lat|220:300, lon|500:610,  time|:)
  region_mean_anomalies =  dim_avg_n(precip_south_china, (/0,1/) )

  precip_sst_reg = regCoef(region_mean_anomalies, sst_anomalies(lat|:, lon|:, time|24:503 ) )
  
  tval = onedtond(precip_sst_reg@tval , dimsizes(precip_sst_reg) )
  df   = onedtond(precip_sst_reg@nptxy, dimsizes(precip_sst_reg)) - 2
  b = tval    ; b must be same size as tval (and df)
  b = 0.5
  prob = betainc(df/(df+tval^2),df/2.0,b)       ; prob(nlat,nlon)
  copy_VarCoords(sst_anomalies(lat|:, lon|:, time|0 ), prob  )

  precip_sst_reg = precip_sst_reg *100
  copy_VarCoords(sst_anomalies(lat|:, lon|:, time|0 ), precip_sst_reg )

  precip_sst_reg_lag3 = regCoef(region_mean_anomalies, sst_anomalies(lat|:, lon|:, time|21:500 ) )

  tval = onedtond(precip_sst_reg_lag3@tval , dimsizes(precip_sst_reg_lag3) )
  df   = onedtond(precip_sst_reg_lag3@nptxy, dimsizes(precip_sst_reg_lag3)) - 2
  b = tval    ; b must be same size as tval (and df)
  b = 0.5
  prob_lag3 = betainc(df/(df+tval^2),df/2.0,b)       ; prob(nlat,nlon)
  copy_VarCoords(sst_anomalies(lat|:, lon|:, time|0 ), prob_lag3  )

  precip_sst_reg_lag3 = precip_sst_reg_lag3 *100
  copy_VarCoords(sst_anomalies(lat|:, lon|:, time|0 ), precip_sst_reg_lag3 )

  ncdf = addfile("precip_regressed_sst.nc" ,"c")  ; open output netCDF file
  ncdf->precip_sst_reg = precip_sst_reg
  ncdf->precip_sst_reg_lag3 = precip_sst_reg_lag3
  ncdf->prob = prob
  ncdf->prob_lag3 = prob_lag3

;*******************************************************************
; create plot for total preciptation and wind_850hpa anomalies
;*******************************************************************
  res                       = True              ; plot mods desired
  res@gsnDraw               = False
  res@gsnFrame              = False
;  res@mpProjection          = "CylindricalEqualArea"; choose projection
  res@mpFillOn              = False             ; turn off map fill
  res@cnFillOn              = True              ; turn on color
  res@cnLinesOn             = False             ; turn off contour lines
  res@cnFillPalette         = "GMT_polar"          ; set color map

  res@cnLevelSelectionMode  = "ManualLevels"    ; set manual contour levels
  res@cnMinLevelValF        =  -1.5              ; set min contour level
  res@cnMaxLevelValF        =   1.5              ; set max contour level
  res@cnLevelSpacingF       =   0.5               ; set contour spacing
  res@gsnAddCyclic          = False             ; regional plot

;  res@mpMinLatF             = 20                ; min lat to mask
;  res@mpMaxLatF             = 80                ; max lat to mask
;  res@mpMinLonF             = -90               ; min lon to mask
;  res@mpMaxLonF             = 40                ; max lon to mask

;  res@mpMaxLatF                   = max(lat)           ; choose subregion           
;  res@mpMinLatF                   = min(lat)
   res@mpMaxLonF                   = 360.
   res@mpMinLonF                   = 0. 
   res@mpCenterLonF                = 180.     ; This is necessary to get the correct map

;  wks = gsn_open_wks("pdf","sst_anomalies")            ; send graphics to PNG file
;  gsn_define_colormap(wks, "GMT_polar")
;  gsn_reverse_colormap(wks)

;   plot = new(12, graphic)
;   do month = 0, 11
;	res@tiMainString = months(month) 
;	plot(month) = gsn_csm_contour_map(wks, sst_anomalies (:,:,month+468), res )          ; draw second plot
;    end do
; 	resP                     = True
;	resP@gsnMaximize         = True                ; maximize plots
;	resP@gsnPanelLabelBar    = True
; 	resP@lbAutoManage        = False
; 	resP@lbLabelFontHeightF  = 0.0075
;	gsn_panel(wks,plot,(/4,3/), resP)

;  res@gsnDraw               = False
;  res@cnMinLevelValF        =  -0.8              ; set min contour level
;  res@cnMaxLevelValF        =   0.8              ; set max contour level
;  res@cnLevelSpacingF       =   0.2               ; set contour spacing
;  res@gsnAddCyclic          =  False             ; regional plot

  opt                  = True
  opt@gsnShadeFillType = "pattern"
  opt@gsnShadeHigh     = 17 

  res2 = True                            ; res2 probability plots

  res2@gsnDraw             = False       ; Do not draw plot
  res2@gsnFrame            = False       ; Do not advance frome
  res2@cnInfoLabelOn       = False       ; turn off info label
  res2@cnLinesOn           = False       ; do not draw contour lines
  res2@cnLineLabelsOn      = False       ; do not draw contour labels
  res2@cnFillScaleF        = 0.6         ; add extra density

  plot = new(2, graphic)
   
  wks1 = gsn_open_wks("pdf","precip_sst_regression")            ; send graphics to PNG file
;  gsn_define_colormap(wks1, "GMT_polar")
;  gsn_reverse_colormap(wks1)

  plot(0)  = gsn_csm_contour_map(wks1, precip_sst_reg , res )          ; draw second plot
  plot11   = gsn_csm_contour_map(wks1, prob , res2 )          ; draw second plot
  plot11   = gsn_contour_shade(plot11, 0.05, 1, opt)  ; shade all areas less than the 0.05 contour level
  overlay (plot(0), plot11)

  wks2 = gsn_open_wks("pdf","precip_sst_regression_lag3")            ; send graphics to PNG file
  gsn_define_colormap(wks2, "GMT_polar")
  gsn_reverse_colormap(wks1)

  plot(1)   = gsn_csm_contour_map(wks1, precip_sst_reg_lag3 , res )          ; draw second plot
  plot22   = gsn_csm_contour_map(wks1, prob_lag3 , res2 )          ; draw second plot
  plot22   = gsn_contour_shade(plot22, 0.05, 0.95, opt)  ; shade all areas less than the 0.05 contour level
  overlay (plot(1), plot22)

  resP                     = True
  resP@gsnPanelLabelBar    = True
  resP@lbAutoManage        = False
  resP@lbLabelFontHeightF  = 0.0075
  gsn_panel(wks1,plot,(/2,1/),resP)

;*********************************************************************
; create plot for integrated moisture flux and omega_500hpa anomalies
;*********************************************************************

end;

