

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

;+++++++++++++++++++++++++++++++++++++++++++++++++++++++
;++++++++++++ advect varible of temperature ++++++++++++
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++

 file_precip_0 = addfile("GPH_T_U_V_500hpa.nc","r")
 t = file_precip_0->t (time|:, expver|0,  latitude|160:360,longitude|240:720)
 t := t(:,::-1,:)
 u = file_precip_0->u (time|:, expver|0,  latitude|160:360,longitude|240:720)
 u := u(:,::-1,:)
 v = file_precip_0->v (time|:, expver|0,  latitude|160:360,longitude|240:720)
 v := v(:,::-1,:)

	lat = file_precip_0->latitude(160:360)
	lon = file_precip_0->longitude(240:720)


   gridType  = 0   ; global gaussian grid
   opt_adv   = 0   ; return only the advected variable
 
   long_name = "advection of temperature"
   units     = "K/s"

   Tadv = advect_variable(short2flt(u), short2flt(v), short2flt(t),gridType,long_name,units,opt_adv)

  temp   :=  Tadv(latitude|:,longitude|:,time|0:503)  
	temp!0 = "lat"
	temp!1 = "lon"
	temp!2 = "time"
  tp_clm := clmMonLLT (short2flt(temp))             ;Calculate LTMM's
  Tadv_anomalies  = calcMonAnomLLT(short2flt(temp), tp_clm) ;Calculate Anomalies from Means

  temp   :=  u(latitude|:,longitude|:,time|0:503)  
	temp!0 = "lat"
	temp!1 = "lon"
	temp!2 = "time"
  tp_clm := clmMonLLT (short2flt(temp))             ;Calculate LTMM's
  u5_anomalies  = calcMonAnomLLT(short2flt(temp), tp_clm) ;Calculate Anomalies from Means

  temp   :=  v(latitude|:,longitude|:,time|0:503)  
	temp!0 = "lat"
	temp!1 = "lon"
	temp!2 = "time"
  tp_clm := clmMonLLT (short2flt(temp))             ;Calculate LTMM's
  v5_anomalies  = calcMonAnomLLT(short2flt(temp), tp_clm) ;Calculate Anomalies from Means


; region_min := dim_min( omega_anomalies)
; region_mean := dim_avg( omega_anomalies)
; region_max  := dim_max( omega_anomalies)

 region_mean := dim_avg(Tadv)
 print (region_mean(100,200))



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
  vcres@vcLevelSpacingF          =  0.5    ; 
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


;*********************************************************************
; create plot for tempature advction and wind_500hpa anomalies
;*********************************************************************

  res@gsnDraw               = False
  res@cnMinLevelValF        =  -2.e-5              ; set min contour level
  res@cnMaxLevelValF        =   2.e-5              ; set max contour level
  res@cnLevelSpacingF       =   5.e-6               ; set contour spacing

  vcres@vcRefLengthF             = 0.045
  vcres@vcRefMagnitudeF          = 4.0
  vcres@vcRefAnnoOrthogonalPosF  = -0.12
  vcres@vcRefAnnoParallelPosF    =  0.997
  vcres@vcRefAnnoFontHeightF     = 0.015

   wks3 = gsn_open_wks("pdf","Temperature_advection_wind_500hpa_anomalies")            ; send graphics to PNG file
  gsn_define_colormap(wks3, "GMT_polar")
;  gsn_reverse_colormap(wks3)


   plot = new(12, graphic)
    do month = 0, 11
	res@tiMainString = months(month) + "(mm/day)"
	plotA= gsn_csm_contour_map(wks3, Tadv_anomalies (:,:,month+492), res )          ; draw second plot
	plotB = gsn_csm_vector(wks3,u5_anomalies(:,:,month+492), v5_anomalies(:,:,month+492),vcres)         ; create a default plot
	overlay(plotA, plotB)
	plot(month) = plotA
;	maximize_output(wks3,True)
    end do
 	resP                     = True
	resP@gsnMaximize    = True                ; maximize plots
	resP@gsnPanelLabelBar    = True
 	resP@lbAutoManage        = False
 	resP@lbLabelFontHeightF  = 0.0075
	gsn_panel(wks3,plot,(/4,3/), resP)

end;

