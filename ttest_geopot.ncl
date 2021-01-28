undef("student_t_test")
function student_t_test(x:numeric, y:numeric, opt_mask:logical)
local dims, ntim, sigr, xtmp, ytmp, xAve, yAve, \
    xVar, yVar, xEqv, yEqv, rad, wgty, xN, yN, prob
begin

  dims = dimsizes(y)
  ntim = dims(0)

  ;--- Student ttest
  xtmp = x(lat|:,lon|:,time|:)
  ytmp = y(lat|:,lon|:,time|:)
  xAve = dim_avg_Wrap (xtmp)
  yAve = dim_avg_Wrap (ytmp)
  xVar = dim_variance_Wrap (xtmp)
  yVar = dim_variance_Wrap (ytmp)

  sigr = 0.05 ; critical sig lvl for r
  xEqv = equiv_sample_size (xtmp, sigr, 0)
  yEqv = equiv_sample_size (ytmp, sigr, 0)

  rad  = 4.0*atan(1.0)/180.0
  wgty = cos(x&lat*rad)
  xN   = wgt_areaave (xEqv, wgty, 1., 0)
  yN   = wgt_areaave (yEqv, wgty, 1., 0)

  prob = ttest(xAve,xVar,xN, yAve,yVar,yN, False, False) 
  copy_VarMeta(x(0,:,:), prob)
  ;printVarSummary(prob)

  return ([/xAve, prob/])
end
  
undef("bootstrap_test")
function bootstrap_test(x:numeric, y:numeric, opt_mask:logical)
local nBoot, stat, nDim, BootStrap, CI_95, CI_05, mean, prob
begin

  ;--- Bootstrap sig. test
  nBoot       = 5000
  stat        = 0 ; average
  nDim        = 0
  BootStrap   = bootstrap_stat(y, stat, nBoot, nDim, False)
  CI_05       = bootstrap_estimate(BootStrap[0], 0.025, False)
  CI_95       = bootstrap_estimate(BootStrap[0], 0.975, False)

  mean        = dim_avg_n_Wrap (x, 0)
  prob        := where(mean.le.CI_05 .or. mean.ge.CI_95, 0.025, 999.)
  if (opt_mask) then
      mean    = where(prob.eq.0.025, mean, mean@_FillValue)
      copy_VarMeta(x(0,:,:), mean)
  end if
  copy_VarMeta(x(0,:,:), prob)

  return([/mean, prob/])
end

begin

  f    = addfile("HGT500_MON_1958-1997.nc","r")
  x    = f->HGT(:,{20:},{40:240})
  printVarSummary(x)

  ; random sample for test
  index = (/0,2,3,5,7,8,10,11,12,19,20,21,22,23,27,28,29,38/) 

  clm = clmMonTLL (x) ; monthly climatology
  xAnom = calcMonAnomTLL (x, clm)

  xDJF = month_to_season(xAnom, "DJF")

  x := xDJF(index,:,:)
  y := xDJF

  output := student_t_test(x, y, False)
  mean1 := output[0]
  prob1 := output[1]

  output := bootstrap_test(x, y, False)
  mean2 := output[0]
  prob2 := output[1]


;*******************************************************
; create plot
;*******************************************************
  wks                      = gsn_open_wks("x11","sample")
  plot_base                = new(2,graphic)
  plot_sig                 = new(2,graphic)

  res                      = True
  res@gsnMaximize          = True
  res@gsnDraw              = False
  res@gsnFrame             = False
  res@gsnAddCyclic         = False

  res2 = res

  res@cnFillOn             = True 
  res@cnFillPalette        = "BlueDarkRed18"
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF       = -30.
  res@cnMaxLevelValF       =  30.
  res@cnLevelSpacingF      =  5.

  res@mpCenterLonF         = 180.
  res@mpMinLatF            = 20.
  res@mpMinLonF            = 40.
  res@mpMaxLonF            = 240
  res@mpShapeMode          = "FreeAspect"
  res@vpWidthF             = 0.5
  res@vpHeightF            = 0.3
  res@gsnLeftString        = ""
  res@gsnRightString       = ""

  res2@cnLevelSelectionMode = "ManualLevels"
  res2@cnMinLevelValF      = 0.00
  res2@cnMaxLevelValF      = 0.5
  res2@cnLevelSpacingF     = 0.01
  res2@cnInfoLabelOn       = False
  res2@cnLinesOn           = False
  res2@cnLineLabelsOn      = False
  res2@gsnLeftString       = ""
  res2@gsnRightString      = ""

  res@tiMainString         = "Student-T" ; title 
  plot_base(0) = gsn_csm_contour_map(wks, mean1, res)

  res@tiMainString         = "Bootstrapping" ; title 
  plot_base(1) = gsn_csm_contour_map(wks, mean2, res)

  opt                      = True
  opt@gsnShadeLow          = 6
  opt@gsnShadeFillType     = "pattern"
  opt@gsnShadeFillScaleF   = 0.5

  plot_sig(0) = gsn_csm_contour(wks, prob1, res2)
  plot_sig(0) = gsn_contour_shade(plot_sig(0), 0.05, 0.8, opt)
  overlay(plot_base(0), plot_sig(0))

  plot_sig(1) = gsn_csm_contour(wks, prob2, res2)
  plot_sig(1) = gsn_contour_shade(plot_sig(1), 0.05, 0.8, opt)
  overlay(plot_base(1), plot_sig(1))

  resP = True
  resP@gsnPanelMainString = "season=DJF | var=HGT500 | case=NA"
  gsn_panel(wks,plot_base,(/2,1/),resP) 

end
