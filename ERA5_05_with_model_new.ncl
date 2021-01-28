;------------------------------------------------------------------------------------------
;  2019.09.25. Gaeun Kim
;------------------------------------------------------------------------------------------
;  PERPOSE:
;  ECMWF Sub-Seasonal hindcast Anomalies
;------------------------------------------------------------------------------------------
;  PRODUCTS:
;  Weekly AO and ENSO composite analysis map (Z500)
;  Time series plot of OBS AO and ENSO index
;  Box plot of KR winter mean temperature during AO and ENSO composite cases
;  Box plot of PCC
;------------------------------------------------------------------------------------------

load "/lfs/data25/gkim/myLIB/myLIB_1_calc_ttest.ncl"
load "/lfs/data25/gkim/AO/test/Procedure_plot1_TimeSeries_AO_and_ENSO.ncl"
;load "/data25/gkim/AO/test/Procedure_plot2_Composite_map.ncl"
load "/lfs/data25/gkim/AO/ForJournal/ERA5_04_composite_res_prob.ncl"


;------------------------------------------------------------------------------------------
; ENV Variables
;------------------------------------------------------------------------------------------
N2_DIR = "/lfs/data25/gkim/AO/data/ERA5/"
OI_DIR = "/lfs/data25/gkim/AO/data/OISST/"
EC_DIR = "/lfs/data25/gkim/EC_hind/"

dStrt = 19990101
dLast = 20111231

cStrt = 19990101
cLast = 20101231

wtype = "ps"

slat = 31.0
elat = 41.0
slon = 121.0
elon = 133.0
jpKor = (/37,36,36,35,35/)
ipKor = (/85,85,86,85,86/)

nlat = 121
nlon = 240

;------------------------------------------------------------------------------------------
; READ OBS DATA
;------------------------------------------------------------------------------------------
;; ERA5 T2M
   ifile  	= addfile(N2_DIR + "/t2m/t2m.nc","r") ; 121*240
   t2m   	= ifile->t2m
   
   xtime  	= ifile->time
   
   ctime  	= cd_calendar(xtime,-2)
   tinds  	= ind(ctime .ge. cStrt .and. ctime .le. cLast)
   
   t2mClm   	= clmDayTLL(t2m(tinds,:,:),yyyymmdd_to_yyyyddd(ctime(tinds)))
   t2mAnom  	= calcDayAnomTLL(t2m,yyyymmdd_to_yyyyddd(ctime),t2mClm)

;; ERA5 HGT 
   ifile        := addfile(N2_DIR + "/z/z.nc","r")
   hgt           = ifile->z
   hgt           = hgt/9.80665

   hgtClm        = clmDayTLL(hgt(tinds,:,:),yyyymmdd_to_yyyyddd(ctime(tinds)))
   hgtAnom       = calcDayAnomTLL(hgt,yyyymmdd_to_yyyyddd(ctime),hgtClm)

;; ERA5 SLP
   ifile	:= addfile(N2_DIR + "/msl/msl.nc","r")
   slp    	 = ifile->msl
   slp		 = slp/100.
    
   slpClm	 = clmDayTLL(slp(tinds,:,:),yyyymmdd_to_yyyyddd(ctime(tinds)))
   slpAnom 	 = calcDayAnomTLL(slp,yyyymmdd_to_yyyyddd(ctime),slpClm)
    

;; AO Index CPC
   ifile  	:= addfile("/lfs/data25/gkim/AO/data/daily_ao_1981-2018.nc","r")
   xtime  	:= ifile->time
   ctime  	:= cd_calendar(xtime,-2)
   tinds  	:= ind(ctime .ge. 19830101 .and. ctime .le. 20181231)
   
   CPC    	 = ifile->ao(tinds)
   ftime  	 = cd_calendar(CPC&time,-2)


;; KR temp
   ifile  	:= addfile("/lfs/data25/gkim/AO/data/KR_temp_45stn.nc","r")
   xtime  	:= ifile->time
   
   temp 	 = new((/dimsizes(xtime),1,1/),float)
   temp(:,0,0)   = ifile->t_std
   printVarSummary(temp)
   
   ctime  	:= cd_calendar(xtime,-2)
   tinds  	:= ind(ctime .ge. cStrt .and. ctime .le. cLast)
   
   tClm    	 = clmDayTLL(temp(tinds,:,:),yyyymmdd_to_yyyyddd(ctime(tinds)))
   tA_dump 	 = calcDayAnomTLL(temp,yyyymmdd_to_yyyyddd(ctime),tClm)
   tAnom   	 = tA_dump(:,0,0)
   printVarSummary(tAnom)


;; OISST
   ifile 	:= addfile(OI_DIR + "sst.day.mean.1983-2018.nc","r") ; 121*240
   
   xtime  	:= ifile->time
   ctime  	:= cd_calendar(xtime,-2)
   tinds  	:= ind(ctime .ge. cStrt .and. ctime .le. cLast)
   
   sst     	 = ifile->sst
   sstClm  	 = clmDayTLL(sst(tinds,:,:),yyyymmdd_to_yyyyddd(ctime(tinds)))
   sstAnom 	 = calcDayAnomTLL(sst,yyyymmdd_to_yyyyddd(ctime),sstClm)
   
   sstAnom_avg 	 = wgt_areaave_Wrap(sstAnom(:,{-5:5},{190:240}),1.0,1.0,1)
   tinds	:= ind(ctime .ge. dStrt .and. ctime .le. dLast)
   
   nrun		 = 150
   sstAnom_avg   = runave_n_Wrap(sstAnom_avg,nrun,1,0)
   sstStd	 = sstAnom_avg
   sstStd        = sstAnom_avg / stddev(sstAnom_avg(tinds))
;   print(cd_calendar(sstStd&time(tinds),-2) + " " + sstStd(tinds))

;-----------------------------------------------------------------------------------
; SET TIME VARIABLES 
;-----------------------------------------------------------------------------------
;; Time Original
   tinds   	       := ind(ctime .ge. dStrt .and. ctime .le. dLast)
   yyyymmdd		= ctime(tinds)
   yyyy     		= yyyymmdd / 10000
   mmdd     		= yyyymmdd - yyyy*10000
   mm       		= mmdd / 100
   dd       		= mmdd - mm*100
   
;; Ignore leap Year -- Date of Model Data Directories
   day29 		= ind( .not. (mm .eq. 2 .and. dd .eq. 29) )
   dir_time 		= yyyymmdd(day29)
   d_yyyy 		= dir_time / 10000
   d_mmdd 		= dir_time - d_yyyy * 10000
   d_mm   		= d_mmdd / 100
   d_dd   		= d_mmdd - d_mm * 100
   
   day1   		= day_of_year(d_yyyy,d_mm,d_dd)
   day1_ind		= ind(day1 .eq. 1)
   
   nfd 			= 42
   days_in_year 	= 365
   nweek 		= 6
   
   yrStrt 		= dStrt / 10000
   yrLast 		= dLast / 10000


;-----------------------------------------------------------------------------------
; SET TIME VARIABLES 
;-----------------------------------------------------------------------------------
   tKor_dum = new( (/nfd,5/), float)
   onedlat = new( (/dimsizes(ipKor)/), float)
   onedlon = onedlat


   cnt = 1
   do ii = 0, dimsizes(day1_ind) - 1
   
      year = yyyy(day1_ind(ii))
      yr_end_ind = day1_ind(ii) + days_in_year - 1
   
      do jj = day1_ind(ii), yr_end_ind, 7			;; date of producing reforecast
    
         if( (d_yyyy(jj) .eq. yrStrt .and. d_mm(jj) .ge. 11) \
   	.or. (d_yyyy(jj) .ne. yrStrt .and. ( d_mm(jj) .eq. 1 .or. d_mm(jj) .eq. 2 .or. d_mm(jj) .eq. 12 ) )) then 
   
             fdate = d_yyyy(jj)*10000 + d_mm(jj)*100 + d_dd(jj)
	     print(cnt + " | " + fdate + "")
	     cnt = cnt + 1
            ; t_match = ind(ftime .eq. fdate) 
            ; print(ftime(t_match:t_match+6))
   
   	  fname1 = EC_DIR + d_yyyy(jj) + "/" + fdate + "/S2S_hind_ECMWF_ano_" + fdate + "_Z500.nc"
   	  fname2 = EC_DIR + d_yyyy(jj) + "/" + fdate + "/S2S_hind_ECMWF_ano_" + fdate + "_T2M.nc"
   	  fname3 = EC_DIR + d_yyyy(jj) + "/" + fdate + "/S2S_hind_ECMWF_ano_" + fdate + "_MSLP.nc"
   	  fname4 = EC_DIR + d_yyyy(jj) + "/" + fdate + "/S2S_hind_ECMWF_ano_" + fdate + "_SST.nc"
	
   	  print(fname1)
   
   	  ifile  := addfile(fname2,"r")
   	  EC_t2m  = ifile->T2M(:,::-1,:)
   	  EC_time = cd_calendar(ifile->time,-2)
   	  printVarSummary(EC_t2m)
   	  print(EC_time)
   
   	  EC_lat = ifile->lat
   	  EC_lon = ifile->lon
   
          ifile := addfile(fname1, "r")
          EC_hgt = ifile->Z500(:,::-1,:)
          EC_hgt = EC_hgt / 9.80665

	  ifile := addfile(fname3,"r")
	  EC_mslp = ifile->MSLP(:,::-1,:)
	  EC_mslp = EC_mslp/100.

	  ifile := addfile(fname4,"r")
	  EC_sst = ifile->SST(:,{-5:5},{190:240})
	  EC_sst_avg = wgt_areaave_Wrap(EC_sst,1.0,1.0,1)
	  printVarSummary(EC_sst_avg)

   	  do i = 0, dimsizes(ipKor) - 1
   		onedlon(i)    = EC_lon(ipKor(i))
   		onedlat(i)    = EC_lat(jpKor(i))
		print(onedlon(i) + " " + onedlat(i) )
   		tKor_dum(:,i) = EC_t2m(:,{onedlat(i)},{onedlon(i)})
   	  end do

   	  tKor = dim_avg_n_Wrap(tKor_dum,1)
   	  printVarSummary(tKor)
   
          t_match = get1Dindex(ftime,EC_time)
   	  print(t_match)
   
          w0_t2m = new((/6,nlat,nlon/),float)
   	  w0_t2m!0 = "time"
	  w0_mslp = w0_t2m
          w0_hgt = w0_t2m
          EC_0_t2m = w0_t2m 
	  EC_0_mslp = w0_t2m
          EC_0_hgt = w0_t2m
   
   	  w0_time = new((/6/),float)
   	  w0_time!0 = "time"
   	  w0_CPC = w0_time
   	  w0_nino = w0_time
   	  w0_tKor = w0_time
   	  EC_0_tKor = w0_time
	  EC_0_nino = w0_time
   
   	  cnt = 0
   	  do ww = 0, nfd -1, 7
   	     print(ftime(t_match(ww):t_match(ww+6)))
   
   	     EC_0_t2m(cnt,:,:)  = dim_avg_n_Wrap(EC_t2m(ww:ww+6,:,:),0)
   	     EC_0_mslp(cnt,:,:) = dim_avg_n_Wrap(EC_mslp(ww:ww+6,:,:),0)
   	     EC_0_hgt(cnt,:,:) = dim_avg_n_Wrap(EC_hgt(ww:ww+6,:,:),0)
   	     EC_0_tKor(cnt)     = dim_avg_Wrap(tKor(ww:ww+6))
   	     EC_0_nino(cnt)     = dim_avg_Wrap(EC_sst_avg(ww:ww+6))
   
             w0_t2m(cnt,:,:)   = dim_avg_n_Wrap(t2mAnom(t_match(ww):t_match(ww+6),:,:),0)
             w0_hgt(cnt,:,:)   = dim_avg_n_Wrap(hgtAnom(t_match(ww):t_match(ww+6),:,:),0)
             w0_mslp(cnt,:,:)   = dim_avg_n_Wrap(slpAnom(t_match(ww):t_match(ww+6),:,:),0)
   
   	     w0_time(cnt) = CPC&time(t_match(ww))
   	     w0_CPC(cnt)  = dim_avg_Wrap(CPC(t_match(ww):t_match(ww+6)))
   	     w0_nino(cnt) = dim_avg_Wrap(sstStd(t_match(ww):t_match(ww+6)))
	     print(w0_time(cnt) + " " + w0_nino(cnt) ) 
   
   	     w0_tKor(cnt) = dim_avg_Wrap(tAnom(t_match(ww):t_match(ww+6)))
   	
   	     cnt = cnt + 1
   	  end do

   	  print(fdate)
          if ( fdate .eq. 19991105) then  ; 19991203
   	       w_time = w0_time
   	       w_CPC  = w0_CPC 
               w_nino = w0_nino
	       w_EC_nino = EC_0_nino
   
   	       w_t2m    = w0_t2m
   	       w_EC_t2m = EC_0_t2m

   	       w_hgt    = w0_hgt
   	       w_EC_hgt = EC_0_hgt

   	       w_mslp    = w0_mslp
   	       w_EC_mslp = EC_0_mslp
   
   	       w_tKor    = w0_tKor
   	       w_EC_tKor = EC_0_tKor
   	  else 
   	       w_time := array_append_record(w_time,w0_time,0)
   	       w_CPC  := array_append_record(w_CPC,w0_CPC,0)
   	       w_nino := array_append_record(w_nino,w0_nino,0)
   	       w_EC_nino := array_append_record(w_EC_nino,EC_0_nino,0)
   
   	       w_t2m    := array_append_record(w_t2m,w0_t2m,0)
   	       w_EC_t2m := array_append_record(w_EC_t2m,EC_0_t2m,0)
   
   	       w_hgt    := array_append_record(w_hgt,w0_hgt,0)
   	       w_EC_hgt := array_append_record(w_EC_hgt,EC_0_hgt,0)

   	       w_mslp    := array_append_record(w_mslp,w0_mslp,0)
   	       w_EC_mslp := array_append_record(w_EC_mslp,EC_0_mslp,0)

   	       w_tKor := array_append_record(w_tKor,w0_tKor,0)
   	       w_EC_tKor := array_append_record(w_EC_tKor,EC_0_tKor,0)
   	  end if

         delete([/w0_time,w0_CPC,w0_nino,w0_tKor,EC_0_hgt,w0_hgt,EC_0_t2m,w0_t2m,w0_mslp,EC_0_mslp,EC_0_nino/])
         end if 

      end do
   end do
;   print(cd_calendar(w_time,-2))

;   printVarSummary(w_CPC)
;   printVarSummary(w_hgt)
;   printVarSummary(w_nino)
;   printVarSummary(w_time)

   w1    = ispan(0,(dimsizes(w_time)-1),6) 			; includes NOV.
   week1 = ispan(24,(dimsizes(w_time)-1),6)			; NOT includes NOV.
;print(week1)


   rad  = 4.0*atan(1.)/180.
   clat = EC_lat
   clat = cos(rad*clat)
   printVarSummary(clat) 

;-------------------------------------------------------------------------------------------------
; Composite Cases
;-------------------------------------------------------------------------------------------------
   print(cd_calendar(w_time(week1),-2) + " | " + w_CPC(week1) + " | " + w_nino(week1))
   
   E_AO_p = ind(w_CPC(week1) .ge. 0.5 .and. w_nino(week1) .ge. 0.5)
   L_AO_n = ind(w_CPC(week1) .le. -0.5 .and. w_nino(week1) .le. -0.5)
   p_AO_p = ind(w_CPC(week1) .ge. 0.5 .and. abs(w_nino(week1)) .lt. 0.5)
   p_AO_n = ind(w_CPC(week1) .le. -0.5 .and. abs(w_nino(week1)) .lt. 0.5)
   
   print(dimsizes(E_AO_p))
   print(dimsizes(L_AO_n))
   print(dimsizes(p_AO_p))
   print(dimsizes(p_AO_n))
 
   print(cd_calendar(w_time(week1(E_AO_p)),-2) + " | " + w_CPC(week1(E_AO_p)) + " | " + w_nino(week1(E_AO_p)))
 
   myList = NewList("lifo")
   ListPush(myList,p_AO_n)
   ListPush(myList,p_AO_p)
   ListPush(myList,L_AO_n)
   ListPush(myList,E_AO_p)
   printVarSummary(myList)
   
;;;-- Call Procedure for Plot 2 - 
;;   Procedure_plot2(wtype,myList,w_time,w1,week1,w_CPC,w_nino,w_tKor,w_EC_tKor,w_hgt,w_mslp,w_EC_hgt,w_EC_mslp)

;-------------------------------------------------------------------------------------------------
; Plot Composite analysis Z500
;-------------------------------------------------------------------------------------------------
   wname = "ECMWF_ERA5_Fig8-11"
   wks = gsn_open_wks(wtype,wname)


   Title = (/"EP","LN","PA","NA"/)

   do li = 0, ListCount(myList) -1
 
      leads = 4
      iarr := new((/leads,dimsizes(myList[li])/),integer) ; index of episods 
      pcc  := new((/leads,dimsizes(myList[li])/),float)   ; pattern correlation of episods
      pcc2 := pcc

      jj = 0
      do j = 0, dimsizes(myList[li]) - 1

        do n = 0, leads -1 
           iarr(n,jj) = w1(myList[li](j)) -(n*5) + 24
           pcc(n,jj) = pattern_cor(w_t2m(iarr(n,jj),:,:),w_EC_t2m(iarr(n,jj),:,:),clat,0)
           pcc2(n,jj) = pattern_cor(w_hgt(iarr(n,jj),:,:),w_EC_hgt(iarr(n,jj),:,:),clat,0)
        end do

        jj = jj + 1
      end do ; done for jj

      ;;-- PCC BOX PLOT
        PCC_stat = new( (/4,30/), float)
        PCC2_stat = PCC_stat
        parr = new( (/4,5/), float)
        parr2 = parr

        do i = 0, 3
           PCC_stat(i,:) = stat_dispersion(pcc(i,:),True)
	   PCC2_stat(i,:) = stat_dispersion(pcc2(i,:),True)
           parr(i,:) = (/PCC_stat(i,2),PCC_stat(i,6),PCC_stat(i,8),PCC_stat(i,10),PCC_stat(i,14)/)
           parr2(i,:) = (/PCC2_stat(i,2),PCC2_stat(i,6),PCC2_stat(i,8),PCC2_stat(i,10),PCC2_stat(i,14)/)
        end do

           pccres                               = True
           pccres@tmXBLabels                    = (/"week1","week2","week3","week4"/)
           pccres@tmXBLabelFontHeightF          = 0.015
           pccres@tmYLLabelFontHeightF          = 0.015
           pccres@tiYAxisString                 = "Pattern Correlation"
           pccres@tmXBMajorLineColor            = "white"

           pccres@trYMinF                       = -.5
           pccres@trYMaxF                       = 1.

           opti                                 = True
           opti@boxWidth                        = .5

           pccplot = boxplot(wks,ispan(1,4,1),parr,opti,pccres,False)
           pccplot2 = boxplot(wks,ispan(1,4,1),parr2,opti,pccres,False)
           ;draw(pccplot)
           ;frame(wks)
           ;draw(pccplot2)
           ;frame(wks)

       ;--- PLOT
      
   res                          = True
   res@gsnFrame                 = False
   res@gsnDraw                  = False

   res@mpMinLatF                = 20.
   res@mpMinLonF                = 40.
   res@mpMaxLonF                = 240
   res@mpCenterLonF             = 180.

   res@cnFillOn                 = True
   res@cnFillPalette            = "BlueWhiteOrangeRed"
   res@cnLinesOn                = False

   res@lbLabelBarOn		= False

   res2 = res
   res3 = res
   res4 = res

   res@cnLevelSelectionMode	= "ManualLevels"
   res@cnMinLevelValF		= -40.
   res@cnMaxLevelValF		= 40.
   res@cnLevelSpacingF		= 5

   res3@cnLevelSelectionMode	= "ManualLevels"
   res3@cnMinLevelValF		= -4.
   res3@cnMaxLevelValF		= 4.
   res3@cnLevelSpacingF		= 0.5

   fplot = new(5,graphic) 
   fplot2 = fplot 

   tt_t2m_obs = calc_t(w_t2m(iarr(0,:),:,:),w_t2m(week1,:,:),clat)
   prb_t2m_obs = gsn_csm_contour(wks, gsn_add_cyclic_point(tt_t2m_obs), pbres)
   prb_t2m_obs = gsn_contour_shade(prb_t2m_obs, 0.05, 30, opt)
   plt_t2m_obs = gsn_csm_contour_map(wks,dim_avg_n_Wrap(w_t2m(iarr(0,:),:,:),0), res3)
   overlay(plt_t2m_obs,prb_t2m_obs)

   tt_hgt_obs = calc_t(w_hgt(iarr(0,:),:,:),w_hgt(week1,:,:),clat)
   prb_hgt_obs = gsn_csm_contour(wks, gsn_add_cyclic_point(tt_hgt_obs), pbres)
   prb_hgt_obs = gsn_contour_shade(prb_hgt_obs, 0.05, 30, opt)
   plt_hgt_obs = gsn_csm_contour_map(wks,dim_avg_n_Wrap(w_hgt(iarr(0,:),:,:),0), res)
   overlay(plt_hgt_obs,prb_hgt_obs)
   fplot(1) = plt_hgt_obs      
   fplot2(1) = plt_t2m_obs      

   plt_t2m_mdl = new(4,graphic)
   plt_hgt_mdl = plt_t2m_mdl
   txid = new(4,graphic)
   amid = txid

;--- PLOT ECMWF
   ;-- PCC inof
   txres                        = True
   txres@txPerimOn              = True
   txres@txFontHeightF          = 0.02
   txres@txBackgroundFillColor  = "white"

   amres                        = True
   amres@amParallelPosF         = 0.5
   amres@amOrthogonalPosF       = -0.5
   amres@amJust                 = "TopRight"


   ;-- PLOT RES
   res2@cnLevelSelectionMode	= "ManualLevels"
   res2@cnMinLevelValF		= -4.
   res2@cnMaxLevelValF		= 4.
   res2@cnLevelSpacingF		= .5

   res4@cnLevelSelectionMode	= "ManualLevels"
   res4@cnMinLevelValF		= -2.
   res4@cnMaxLevelValF		= 2.
   res4@cnLevelSpacingF		= .2

   do i = 0, 3

       tt = calc_t(w_EC_t2m(iarr(i,:),:,:),w_EC_t2m(w1,:,:),clat)
       prb_tt = gsn_csm_contour(wks, gsn_add_cyclic_point(tt), pbres)
       prb_tt = gsn_contour_shade(prb_tt, 0.05, 30, opt)

       plt_t2m_mdl(i)  = gsn_csm_contour_map(wks, dim_avg_n_Wrap(w_EC_t2m(iarr(i,:),:,:),0), res4)
       overlay(plt_t2m_mdl(i), prb_tt) 

       txid(i) = gsn_create_text(wks, sprintf("%3.2f",parr(i,2)),txres)
       amid(i) = gsn_add_annotation(plt_t2m_mdl(i),txid(i),amres)

       tt2 = calc_t(w_EC_hgt(iarr(i,:),:,:),w_EC_hgt(w1,:,:),clat)
       prb_tt2 = gsn_csm_contour(wks, gsn_add_cyclic_point(tt2), pbres)
       prb_tt2 = gsn_contour_shade(prb_tt2, 0.05, 30, opt)

       plt_hgt_mdl(i)  = gsn_csm_contour_map(wks, dim_avg_n_Wrap(w_EC_hgt(iarr(i,:),:,:),0), res2)
       overlay(plt_hgt_mdl(i), prb_tt2)

       txid(i) = gsn_create_text(wks, sprintf("%3.2f",parr2(i,2)),txres)
       amid(i) = gsn_add_annotation(plt_hgt_mdl(i),txid(i),amres)
   end do

   fplot(0) = plt_hgt_mdl(0)
   fplot(2) = plt_hgt_mdl(1)
   fplot(3) = plt_hgt_mdl(2)
   fplot(4) = plt_hgt_mdl(3)

   fplot2(0) = plt_t2m_mdl(0)
   fplot2(2) = plt_t2m_mdl(1)
   fplot2(3) = plt_t2m_mdl(2)
   fplot2(4) = plt_t2m_mdl(3)

   pres                         = True
   pres@gsnPanelLabelBar        = True
   pres@gsnPanelCenter          = False
   pres@gsnPanelRowSpec         = True
;   gsn_panel(wks, plt_t2m_mdl,(/4,1/),pres)
;   gsn_panel(wks, plt_hgt_mdl,(/4,1/),pres)
   gsn_panel(wks, fplot, (/2,1,1,1/),pres)
   gsn_panel(wks, fplot2, (/2,1,1,1/),pres)

 delete([/res,res2,res3,res4,pres,amres,txres/])

 end do ; done for li 



delete(wks)
oname = wname + ".ps"
pname = wname + ".png"

system("convert -alpha off -background white -density 500 -trim " + oname + " " + pname)
system("rm " + oname)
 

