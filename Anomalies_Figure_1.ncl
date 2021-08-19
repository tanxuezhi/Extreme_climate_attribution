

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

	  dimice=dimsizes(dat)
	  ntim2=dimice(0)
	  nlat2=dimice(1)
 	  mlon2=dimice(2)

;====================================
; Calculate the JFM Clim.
;====================================
  temp   =  dat(latitude|:,longitude|:,time|:)  
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

;************************************************
; 850 hpa wind anomalies 
;************************************************
file_precip_0 = addfile("wind_850hpa.nc","r")
u = file_precip_0->u(time|0:503, expver|0, latitude|160:360,longitude|240:720)
v = file_precip_0->v(time|0:503, expver|0, latitude|160:360,longitude|240:720)
  temp   :=  u(latitude|:,longitude|:,time|:)  
	temp!0 = "lat"
	temp!1 = "lon"
	temp!2 = "time"
  tp_clm := clmMonLLT (short2flt(temp))             ;Calculate LTMM's
  u_anomalies  = calcMonAnomLLT(short2flt(temp), tp_clm) ;Calculate Anomalies from Means

  temp   :=  v(latitude|:,longitude|:,time|:)  
	temp!0 = "lat"
	temp!1 = "lon"
	temp!2 = "time"
  tp_clm := clmMonLLT (short2flt(temp))             ;Calculate LTMM's
  v_anomalies  = calcMonAnomLLT(short2flt(temp), tp_clm) ;Calculate Anomalies from Means

;  region_mean = dim_avg (tp_clm)
region_min := dim_min(u_anomalies)
region_mean := dim_avg(u_anomalies)
region_max  := dim_max(u_anomalies)

region_min := dim_min(v_anomalies)
region_mean := dim_avg(v_anomalies)
region_max  := dim_max(v_anomalies)


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
  cmap = read_colormap_file ("GMT_polar")
  res@cnFillPalette         =  cmap(::-1,:)     ; reverse color map         ; set color map

  res@cnLevelSelectionMode  = "ManualLevels"    ; set manual contour levels
  res@cnMinLevelValF        =  -300.              ; set min contour level
  res@cnMaxLevelValF        =   300.              ; set max contour level
  res@cnLevelSpacingF       =   50.               ; set contour spacing
  res@gsnAddCyclic          = False             ; regional plot

;  res@mpMinLatF             = 20                ; min lat to mask
;  res@mpMaxLatF             = 80                ; max lat to mask
;  res@mpMinLonF             = -90               ; min lon to mask
;  res@mpMaxLonF             = 40                ; max lon to mask

  res@mpMaxLatF                   = max(lat)           ; choose subregion           
  res@mpMinLatF                   = min(lat)
  res@mpMaxLonF                   = 150
  res@mpMinLonF                   = min(lon) 

  vcres = True
  vcres@gsnDraw   = False
  vcres@gsnFrame  = False
  vcres@vcGlyphStyle     = "CurlyVector"   ; curly vectors

;---Set up some vector resources.
  vcres@vcLevelSelectionMode     = "ManualLevels"
  vcres@vcMinLevelValF           = -1.5
  vcres@vcMaxLevelValF           =  1.5
  vcres@vcLevelSpacingF          =  0.3    ; 
;  vcres@vcLevelPalette           = "amwg_blueyellowred"   ; assign color map to vectors
;  vcres@lbLabelBarOn                = False

;---Vector lengths and color
;  vcres@vcFillArrowsOn           = True
;  vcres@vcLineArrowThicknessF    =  2.0

  vcres@vcMinFracLengthF         = 0.33
  vcres@vcMinMagnitudeF          = 10
  vcres@vcMinDistanceF           = 0.03
;  vcres@vcMonoFillArrowFillColor = False
;  vcres@vcMonoLineArrowColor     = False

  vcres@vcRefLengthF             = 0.045
  vcres@vcRefMagnitudeF          = 3.0
  vcres@vcRefAnnoOrthogonalPosF  = -0.12
  vcres@vcRefAnnoParallelPosF    =  0.997
  vcres@vcRefAnnoFontHeightF     = 0.015

  vcres@lbTitleString            = "wind flux"
  vcres@lbTitleOffsetF           = -0.25
  vcres@lbTitleFontHeightF       = 0.02
  vcres@lbLabelFontHeightF       = 0.015
  vcres@lbLabelAutoStride        = True

;---Make sure vectors are drawn in "Postdraw" phase.
  vcres@vcVectorDrawOrder        = "Postdraw"

  wks = gsn_open_wks("pdf","preciptation_wind_anomalies")            ; send graphics to PNG file
;  gsn_define_colormap(wks, "GMT_polar")
;  gsn_draw_colormap(wks) 

  plot = new(12, graphic)
    do month = 0, 11
	res@tiMainString = months(month) 
	plotA= gsn_csm_contour_map(wks, tp_anomalies (:,:,month+468), res )          ; draw second plot
	plotB = gsn_csm_vector(wks,u_anomalies(:,:, month+492), v_anomalies(:,:, month+492),vcres)         ; create a default plot
	overlay(plotA, plotB)
	plot(month) = plotA
;	maximize_output(wks1,True)
    end do
 	resP                     = True
	resP@gsnMaximize         = True                ; maximize plots
	resP@gsnPanelLabelBar    = True
 	resP@lbAutoManage        = False
 	resP@lbLabelFontHeightF  = 0.0075
	gsn_panel(wks,plot,(/4,3/), resP)
  gsn_reverse_colormap(wks)


;*********************************************************************
; create plot for integrated moisture flux and omega_500hpa anomalies
;*********************************************************************

   file_precip_0 = addfile("Asia_moisture_uq.nc","r")
   uq = file_precip_0->uq
   file_precip_0 = addfile("Asia_moisture_vq.nc","r")
   vq = file_precip_0->vq
;   file_precip_0 = addfile("Asia_moisture_dq.nc","r")
;   dq = file_precip_0->dq

  temp   :=  uq(lat|:,lon|:,time|0:503,expver|0)  
  uq_1  = temp
  tp_clm := clmMonLLT (temp)             ;Calculate LTMM's
  uq_anomalies  = calcMonAnomLLT(temp, tp_clm) ;Calculate Anomalies from Means

; region_min := dim_min(temp)
; region_mean := dim_avg(temp)
; region_max  := dim_max(temp)

  temp   :=  vq(lat|:,lon|:,time|0:503,expver|0)  
  vq_1  = temp
  tp_clm := clmMonLLT (temp)             ;Calculate LTMM's
  vq_anomalies  = calcMonAnomLLT(temp, tp_clm) ;Calculate Anomalies from Means

; region_min := dim_min(temp)
 region_mean := dim_avg(temp)
; region_max  := dim_max(temp)

print (region_mean(200,200))

 file_precip_0 = addfile("omega.nc","r")
 omega = file_precip_0->w (time|:, expver|0, level|30, latitude|160:360,longitude|240:720)

  temp   :=  omega(latitude|:,longitude|:,time|0:503)  
	temp!0 = "lat"
	temp!1 = "lon"
	temp!2 = "time"
  omega_1 = temp
  tp_clm := clmMonLLT (short2flt(temp))             ;Calculate LTMM's
  omega_anomalies  = calcMonAnomLLT(short2flt(temp), tp_clm) ;Calculate Anomalies from Means

; region_min := dim_min( omega_anomalies)
; region_mean := dim_avg( omega_anomalies)
; region_max  := dim_max( omega_anomalies)

 region_mean := dim_avg( omega_1)
 print (region_mean(100,200))

 region_anomalies_mean := dim_avg( omega_anomalies )
 print (region_anomalies_mean(100,200))

 region_anomalies_max := dim_max( omega_anomalies )
 print (region_anomalies_max(100,200))

 region_anomalies_min := dim_min( omega_anomalies )
 print (region_anomalies_min(100,200))

;---Set up some vector resources.
  vcres@vcLevelSelectionMode     = "AutomaticLevels"
  vcres@vcMinLevelValF           = -300.0
  vcres@vcMaxLevelValF           =  300.0
  vcres@vcLevelSpacingF          =  50.0    ; 
  vcres@vcRefLengthF             =  0.03
  vcres@vcRefMagnitudeF          =  200.0

  res@gsnDraw               = False
  res@cnMinLevelValF        =  -0.06              ; set min contour level
  res@cnMaxLevelValF        =   0.06              ; set max contour level
  res@cnLevelSpacingF       =   0.01               ; set contour spacing

   wks1 = gsn_open_wks("pdf","moisture_omega_anomalies")            ; send graphics to PNG file
  gsn_define_colormap(wks1, "GMT_polar")
;  gsn_reverse_colormap(wks1)

   plot = new(12, graphic)
    do month = 0, 11
	res@tiMainString = months(month) + "(mm/day)"
	plotA= gsn_csm_contour_map(wks1, omega_anomalies (:,:,month+492), res )          ; draw second plot
	plotB = gsn_csm_vector(wks1,uq_anomalies(:,:,month+492), vq_anomalies(:,:,month+492),vcres)         ; create a default plot
	overlay(plotA, plotB)
	plot(month) = plotA
;	maximize_output(wks1,True)
    end do
 	resP                     = True
	resP@gsnMaximize    = True                ; maximize plots
	resP@gsnPanelLabelBar    = True
 	resP@lbAutoManage        = False
 	resP@lbLabelFontHeightF  = 0.0075
	gsn_panel(wks1,plot,(/4,3/), resP)

  res@cnMinLevelValF        =  -2000.              ; set min contour level
  res@cnMaxLevelValF        =   0.             ; set max contour level
  res@cnLevelSpacingF       =   400.              ; set contour spacing
   
 wks2 = gsn_open_wks("pdf","moisture_omega")            ; send graphics to PNG file
  gsn_define_colormap(wks2, "GMT_polar")
;  gsn_reverse_colormap(wks)

    plot = new(12, graphic)
    do month = 0, 11
	res@tiMainString = months(month) + "(mm/day)"
	plotA= gsn_csm_contour_map(wks2, omega_1 (:,:,month+492), res )          ; draw second plot
	plotB = gsn_csm_vector(wks2, uq_1(:,:,month+492), vq_1(:,:,month+492),vcres)         ; create a default plot
	overlay(plotA, plotB)
	plot(month) = plotA
;	maximize_output(wks2,True)
    end do
 	resP                     = True
	resP@gsnMaximize    = True                ; maximize plots
	resP@gsnPanelLabelBar    = True
 	resP@lbAutoManage        = False
 	resP@lbLabelFontHeightF  = 0.0075
	gsn_panel(wks2,plot,(/4,3/), resP)


end;

