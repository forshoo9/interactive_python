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

  if (opt_mask) then
      mean    = where(prob.le.0.05, xAve, xAve@_FillValue)
      copy_VarMeta(x(0,:,:), mean)
  end if

  return ([/mean, prob/])
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

  f    = addfile("uwnd.10m.mon.mean.nc","r")
  xu   = f->uwnd(:,0,{20:90},{40:240})

  f    = addfile("vwnd.10m.mon.mean.nc","r")
  xv   = f->vwnd(:,0,{20:90},{40:240})

  ; random sample for test
  index = (/0,1,2,6,7,8,17,22,25,26,30,31,33/) 

  clm = clmMonTLL (xu) ; monthly climatology
  Anom = calcMonAnomTLL (xu, clm)
  DJF = month_to_season(Anom, "DJF")

  x := DJF(index,:,:)
  y := DJF

  output := student_t_test(x, y, True)
  umean1 := output[0]
  uprob1 := output[1]

  output := bootstrap_test(x, y, True)
  umean2 := output[0]
  uprob2 := output[1]

  clm = clmMonTLL (xv) ; monthly climatology
  Anom = calcMonAnomTLL (xv, clm)
  DJF = month_to_season(Anom, "DJF")

  x := DJF(index,:,:)
  y := DJF

  output := student_t_test(x, y, True)
  vmean1 := output[0]
  vprob1 := output[1]

  output := bootstrap_test(x, y, True)
  vmean2 := output[0]
  vprob2 := output[1]


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
  res@mpCenterLonF         = 180.
  res@mpMinLatF            = 20.
  res@mpMinLonF            = 40.
  res@mpMaxLonF            = 240
  res@mpShapeMode          = "FreeAspect"
  res@vpWidthF             = 0.5
  res@vpHeightF            = 0.3
  res@gsnLeftString        = ""
  res@gsnRightString       = ""

  vcres = res

  res@tiMainString         = "Student-T" ; title 
  plot_base(0) = gsn_csm_vector_map(wks, umean1, vmean1, vcres)

  res@tiMainString         = "Bootstrapping" ; title 
  plot_base(1) = gsn_csm_vector_map(wks, umean2, vmean2, vcres)

;  opt                      = True
;  opt@gsnShadeLow          = 6
;  opt@gsnShadeFillType     = "pattern"
;  opt@gsnShadeFillScaleF   = 0.5
;
;  plot_sig(0) = gsn_csm_contour(wks, prob1, res2)
;  plot_sig(0) = gsn_contour_shade(plot_sig(0), 0.05, 0.8, opt)
;  overlay(plot_base(0), plot_sig(0))
;
;  plot_sig(1) = gsn_csm_contour(wks, prob2, res2)
;  plot_sig(1) = gsn_contour_shade(plot_sig(1), 0.05, 0.8, opt)
;  overlay(plot_base(1), plot_sig(1))

  resP = True
  resP@gsnPanelMainString = "season=DJF | var=HGT500 | case=NA"
  gsn_panel(wks,plot_base,(/2,1/),resP) 

end
