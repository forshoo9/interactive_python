begin
  ;f    = addfile("uwnd.10m.mon.mean.nc","r")
  ;x    = f->uwnd(:,{20:},{40:240})


  ;--- Bootstrap sig. test
  nBoot       = 5000                ; user set
  stat        = 0                   ; (/0,0/) since they refer to the same dimension
  nDim        = 0                   ; (/0,0/) since they refer to the same dimension

  index = (/0,1,2,6,7,8,17,22,25,26,30,31,33/) ; random for test

  f    = addfile("uwnd.10m.mon.mean.nc","r")
  x    = f->uwnd(:,0,{20:90},{40:240})
  printVarSummary(x)

  clm = clmMonTLL (x) ; monthly climatology
  xAnom = calcMonAnomTLL (x, clm)
  xDJF = month_to_season(xAnom, "DJF")

  BootStrap   = bootstrap_stat(xDJF, stat, nBoot, nDim, False)
  diffBootLow = bootstrap_estimate(BootStrap[0], 0.025, False)
  diffBootHigh = bootstrap_estimate(BootStrap[0], 0.975, False)

  xAve = dim_avg_n_Wrap (xDJF(index,:,:), 0) ; average for the composite case
  udiff = where(xAve .le. diffBootLow .or. xAve .ge. diffBootHigh, 0.25, 999.)
  uwnd = where(xAve .le. diffBootLow .or. xAve .ge. diffBootHigh, xAve, xAve@_FillValue)
  copy_VarMeta(xAve, udiff)
  copy_VarMeta(xAve, uwnd)


  f    = addfile("vwnd.10m.mon.mean.nc","r")
  x    = f->vwnd(:,0,{20:90},{40:240})
  printVarSummary(x)

  clm = clmMonTLL (x) ; monthly climatology
  xAnom = calcMonAnomTLL (x, clm)
  xDJF = month_to_season(xAnom, "DJF")

  BootStrap   = bootstrap_stat(xDJF, stat, nBoot, nDim, False)
  diffBootLow = bootstrap_estimate(BootStrap[0], 0.025, False)
  diffBootHigh = bootstrap_estimate(BootStrap[0], 0.975, False)

  xAve = dim_avg_n_Wrap (xDJF(index,:,:), 0) ; average for the composite case
  vdiff = where(xAve .le. diffBootLow .or. xAve .ge. diffBootHigh, 0.25, 999.)
  vwnd = where(xAve .le. diffBootLow .or. xAve .ge. diffBootHigh, xAve, xAve@_FillValue)
  copy_VarMeta(xAve, vdiff)
  copy_VarMeta(xAve, vwnd)



;*******************************************************
; create plot
;*******************************************************
  wks    = gsn_open_wks("x11","conOncon")            ; send graphics to PNG file
  ;plots  = new(2,graphic)                           ; create graphical array
  colors = (/ "red", "green", "blue", "yellow", "cyan", "hotpink", \
  "red4", "skyblue", "navyblue", "lightyellow", "mediumorchid4", "orange", \
  "slateblue3", "palegreen", "magenta", "springgreen", "pink", \
  "forestgreen", "violet" /)

  res                      = True        ; plot mods desired
  res@gsnMaximize          = True        ; make as large as possible
  res@gsnDraw              = False       ; do not draw until the end
  res@gsnFrame             = False       ; do not automatically advance 'frame'
  res@gsnAddCyclic = False
  ;res@gsnPolar             = "NH"        ; choose which hemisphere
;  res@cnInfoLabelOn        = False            ; no info label
;  res@cnLineLabelsOn       = False            ; no line labels

  res2 = res
  vcres = res

  res@cnFillOn             = True 
  res@cnLinesOn = False
  res@cnFillPalette        = "BlueDarkRed18"
  res@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
  ;res@gsnContourNegLineDashPattern = 1 
  ;res@cnLineLabelPlacementMode = "constant"
  res@cnMinLevelValF       = -0.5        ; set min contour level
  res@cnMaxLevelValF       =  0.5        ; set max contour level
  res@cnLevelSpacingF      =  0.05        ; set contour spacing

  res@mpCenterLonF         = 180.
  res@mpMinLatF            = 20.
  res@mpMinLonF            = 40.
  res@mpMaxLonF            = 240
  res@mpShapeMode          = "FreeAspect"
  res@vpWidthF             = 0.5 ; Resize to a square
  res@vpHeightF            = 0.3
  res@gsnLeftString = ""
  res@gsnRightString = ""

; 1st plot generates the map + 1st contour

  ;res@cnLineThicknessF     = 3                ; thickness of contour lines
  res@tiMainString         = "Student-T" ; title 

  res2@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
  res2@cnMinLevelValF      = 0.00        ; set min contour level
  res2@cnMaxLevelValF      = 0.5        ; set max contour level
  res2@cnLevelSpacingF     = 0.01        ; set contour spacing
  res2@cnInfoLabelOn       = False
  res2@cnLinesOn           = False       ; do not draw contour lines
  res2@cnLineLabelsOn      = False       ; do not draw contour labels
  res2@gsnLeftString = ""
  res2@gsnRightString = ""

  plot = gsn_csm_contour_map(wks, uwnd, res)

  res@tiMainString         = "Bootstrapping" ; title 
  
;  opt                      = True
;  opt@gsnShadeMid          = 6           ; cross hatch pattern
;  opt@gsnShadeFillType     = "pattern"
;  opt@gsnShadeFillScaleF   = 0.5         ; add extra density
;  plot_sig1 = gsn_csm_contour(wks, udiff, res2)
;  plot_sig1 = gsn_contour_shade(plot_sig1, 0., 0.8, opt)
;  overlay(plots(0), plot_sig1)

  vectors1=gsn_csm_vector(wks, uwnd, vwnd, vcres)
  overlay(plot, vectors1)

;  plots(1) = gsn_csm_contour_map(wks, vwnd, res)
;
;  opt                      = True
;  opt@gsnShadeMid          = 6           ; cross hatch pattern
;  opt@gsnShadeFillType     = "pattern"
;  opt@gsnShadeFillScaleF   = 0.5         ; add extra density
;  plot_sig2 = gsn_csm_contour(wks, vdiff, res2)
;  plot_sig2 = gsn_contour_shade(plot_sig2, 0., 0.8, opt)
;  overlay(plots(1), plot_sig2)

;  res@tiMainString         = "Bootstrap 5% | nboot=5000" ; title 
;  plots(2) = gsn_csm_contour_map(wks, diffBootLow, res)
;  res@tiMainString         = "Bootstrap 95% | nboot=5000" ; title 
;  plots(3) = gsn_csm_contour_map(wks, diffBootHigh, res)
  ;plots(1) = gsn_csm_contour_map(wks, diffBootAvg, res)

;  opt@gsnShadeLow          = 6           ; cross hatch pattern
;  opt@gsnShadeFillType     = "pattern"
;  opt@gsnShadeFillScaleF   = 0.3         ; add extra density
;  plot_sig = gsn_csm_contour(wks, prob, res2)
;  plot_sig = gsn_contour_shade(plot_sig, 0.05, 0.8, opt)
;  overlay(plots(0), plot_sig)


;  resP = True
;  resP@gsnPanelMainString = "season=DJF | var=HGT500 | case=NA"
;  ;resP@gsnPanelLabelBar    = True                ; add common colorbar
;  resP@lbLabelFontHeightF  = 0.007               ; make labels smaller 
;  gsn_panel(wks,plots,(/1,1/),resP) 

  draw(plot)
  frame(wks)

;  plot_base = gsn_csm_contour_map_polar(wks,prob,res)

;  delete (res@gsnPolar)        ; avoid annoying warning messages
;;*******************************************
;; loop over other fields but only do contour
;; note the color index is changing. here we are selecting a new
;; color from the default color map.
;;******************************************
;
;  res@gsnLeftString    = ""  ; suppress automatic labels
;  res@gsnRightString   = ""
;  res@cnLineThicknessF = 1
;  
;  do i=0,18                    ; loop over the number of plots
;    res@cnLineColor      = colors(i)     ; change line color 
;    plots(i) = gsn_csm_contour(wks,x(12*i+1,:,:),res)
;    overlay(plot_base,plots(i))
;  end do

end
