pro oned_spec_std

;dir = '/Users/jhyoon/Tools/MosfireDRP-1.0/REDUX/mosFv2/2013jan4/Y/'
chr_rdtbl,'fname_std.list',0,flist
;chr_rdtbl,'fname_std_var.list',0,flist_var
N_file = N_elements(flist)

chr_rdtbl,'fname_mosR_star.list',0,flist_star
N_star = N_elements(flist_star)

@plot_setting
!p.charsize=1.2
;goto, skip
device,file='oned_spec_std.ps',/color,/cmyk,xs=15,ys=12

;;;
;;; stars in R123 masks
;;;
interval = 50
openw,1,'slit_loss_star2.dat'
ypix  = [70,31,48,47,54,40,71]
dy = [15,15,15,15,15,9]
for i=0,N_star-1 do begin
	print,i,flist_star[i]
	f_eps = flist_star[i]+'_eps.fits'
	;;;
	;;; psf fit and compute slit loss
	;;;
	;spec2d = mrdfits(f_eps,0,hdr)
	fits_read,f_eps,spec2d,hdr
	yy = indgen(dy[i]*2+1)
	pix_scale = 4.985e-5*3600 ; arcsec	
	x_all = []
	slit_loss_all = []
	for j=0+200,N_elements(spec2d[*,0])-200,interval do begin
		prof = spec2d[j,ypix[i]-dy[i]:ypix[i]+dy[i]]
		;;; fit psf
		prof_fit = gaussfit(yy,prof,coeff,nterms=3)
		prof_fit_lorentzian = mpfitpeak(yy,prof,A_lorentzian,/lorentzian)
		prof_fit_moffat = mpfitpeak(yy,prof,A_moffat,/moffat)

		;;; compute 1D slit loss gaussian and moffat psf
		slit_left = A_moffat[1]-0.5/pix_scale
		slit_left2 = A_moffat[1]-3*0.5/pix_scale
		slit_right = A_moffat[1]+0.5/pix_scale
		slit_right2 = A_moffat[1]+3*0.5/pix_scale
		vline,slit_left,linestyle=1
		vline,slit_right,linestyle=1
		sub_in_slit = where(yy ge slit_left and yy le slit_right)
		sub_out_slit = where((yy gt slit_left2 and yy lt slit_left) or $
									(yy gt slit_right and yy lt slit_right2))
		flux_in_slit = total(prof[sub_in_slit])
		flux_out_slit = total(prof[sub_out_slit])
		flux_in_slit_moffat = total(prof_fit_moffat[sub_in_slit])
		flux_out_slit_moffat = total(prof_fit_moffat[sub_out_slit])
		slit_loss = flux_out_slit/(flux_in_slit+flux_out_slit)*100
		slit_loss_moffat = flux_out_slit_moffat / $
								(flux_in_slit_moffat+flux_out_slit_moffat)*100
		x_all = [x_all,j]
		slit_loss_all = [slit_loss_all,slit_loss_moffat]

		;;; plot
		if j/100 eq j/100. then begin
			plot,yy,prof,psym=10,xr=[coeff[1]-10,coeff[1]+10],$
				tit=strmid(flist_star[i],61,25)+',x='+strtrim(j,2)
			oplot,yy,prof_fit,color=255
			oplot,yy,prof_fit_lorentzian,color=150,linestyle=1
			oplot,yy,prof_fit_moffat,color=70,linestyle=2
			legend,['center='+strtrim(coeff[1],2)+','+strtrim(A_moffat[1],2),$
					'FWHM='+strtrim(coeff[2]*2.3548*pix_scale,2)+'"',$
					'slit_loss='+strtrim(slit_loss,2)+'%',$
					'slit_loss(moffat)='+strtrim(slit_loss_moffat,2)+'%'],box=0
		endif
	endfor
	plot,x_all,slit_loss_all,psym=7,yr=[0,100],$
		xtit='x (pix)',ytit='slit loss (%)'
	median_slit_loss_all = median(slit_loss_all,/even)
	mean_slit_loss_all = mean(slit_loss_all)
	hline,median_slit_loss_all,linestyle=2
	hline,mean_slit_loss_all,linestyle=1
	legend,['median slit loss='+string(median_slit_loss_all,f='(f5.2)')+'%',$
			'mean slit loss='+string(mean_slit_loss_all,f='(f5.2)')+'%'],$
			box=0
	printf,1,strmid(flist_star[i],49,99),median_slit_loss_all
endfor
close,1

;;;
;;; standard stars
;;;
ypix = [1017,1062,1007,1073,1005,1006,1072,1072,1005,1072,$
		1017,1061,1017,1018,1062,1062,1018,1062,1018,1063,$
		1011,1065,1010,1066,1010,1065,1011,1065,1011,1065,1011,1065,1011,1065,1011,1065,1011,1065,1011,1065,1011,1065]
dy = [20,20,25,25,25,25,25,25,35,30,$
		20,20,20,20,20,20,20,20,15,20,$
		25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25]

openw,1,'slit_loss_std2.dat'
for i=0,N_file-1 do begin
	print,flist[i]
	f_eps = flist[i]+'.fits.gz'
	;f_out = strmid(f_eps,69,50)+'_1dspec.dat'
	f_out = strmid(flist[i],88,30)+'_'+strmid(flist[i],57,9)+'_1dspec.dat'
	;;;
	;;; Horne extraction of 1D spectra
	;;;
	horne_extract_all,f_eps,f_sig,f_out,$
					ypix=ypix[i],dy=dy[i]
	
	readcol,f_out,lambda,flux,flux_err
	
	plot,lambda,flux,psym=10,xr=[9750,11200],$
		xtit=textoidl('\lambda (\AA)'),ytit='flux',$
		tit=strmid(flist[i],69,50)
	oplot,lambda,flux_err,color=255,psym=10
	;;;
	;;; psf fit and compute slit loss
	;;;
	spec2d = mrdfits(f_eps,0,hdr)
	yy = indgen(dy[i]*2+1)
	pix_scale = 4.985e-5*3600 ; arcsec	
	x_all = []
	slit_loss_all = []
	for j=0+200,N_elements(spec2d[*,0])-200,interval do begin
		prof = spec2d[j,ypix[i]-dy[i]:ypix[i]+dy[i]]
		;;; fit psf
		prof_fit = gaussfit(yy,prof,coeff,nterms=3)
		prof_fit_lorentzian = mpfitpeak(yy,prof,A_lorentzian,/lorentzian)
		prof_fit_moffat = mpfitpeak(yy,prof,A_moffat,/moffat)

		;;; compute 1D slit loss gaussian and moffat psf
		slit_left = A_moffat[1]-0.6/pix_scale
		slit_left2 = A_moffat[1]-3*0.6/pix_scale
		slit_right = A_moffat[1]+0.6/pix_scale
		slit_right2 = A_moffat[1]+3*0.6/pix_scale
		vline,slit_left,linestyle=1
		vline,slit_right,linestyle=1
		sub_in_slit = where(yy ge slit_left and yy le slit_right)
		sub_out_slit = where((yy gt slit_left2 and yy lt slit_left) or $
									(yy gt slit_right and yy lt slit_right2))
		flux_in_slit = total(prof[sub_in_slit])
		flux_out_slit = total(prof[sub_out_slit])
		flux_in_slit_moffat = total(prof_fit_moffat[sub_in_slit])
		flux_out_slit_moffat = total(prof_fit_moffat[sub_out_slit])
		slit_loss = flux_out_slit/(flux_in_slit+flux_out_slit)*100
		slit_loss_moffat = flux_out_slit_moffat / $
								(flux_in_slit_moffat+flux_out_slit_moffat)*100
		x_all = [x_all,j]
		slit_loss_all = [slit_loss_all,slit_loss_moffat]
		;;; build 2D moffat psf
;		Npts = 200
;		Nsigma = 10.
;		x2d = (findgen(Npts+1)/(Npts/2.)-1.)*Nsigma*A_moffat[2] + A_moffat[1]
;		y2d = (findgen(Npts+1)/(Npts/2.)-1.)*Nsigma*A_moffat[2] + A_moffat[1]
;print,min(x2d),max(x2d),x2d[10]-x2d[9]
;		moffat_2d = dblarr(Npts+1,Npts+1)
;		for y_i=0,Npts do begin
;			u_2d = ((x2d - A_moffat[1])/A_moffat[2])^2. + $
;					((y2d[y_i] - A_moffat[1])/A_moffat[2])^2.
;			moffat_2d[*,y_i] = A_moffat[0]/(u_2d+1.)^A_moffat[3]
;		endfor
;		slit_left_2d = A_moffat[1]-0.6/pix_scale;A_moffat[2]*Nsigma
;		slit_right_2d = A_moffat[1]+0.6/pix_scale;A_moffat[2]*Nsigma
;		sub_in_slit_2d = where(x2d ge slit_left_2d and x2d le slit_right_2d)
;		sum_moffat_2d = total(moffat_2d)
;		sum_moffat_in_slit_2d = total(moffat_2d[sub_in_slit_2d,*])
;		slit_loss_2d = (1. - sum_moffat_in_slit_2d/sum_moffat_2d)*100.
	
		if j/100 eq j/100. then begin
			plot,yy,prof,psym=10,xr=[coeff[1]-10,coeff[1]+10],tit=strmid(flist[i],88,25)+',x='+strtrim(j,2)
			oplot,yy,prof_fit,color=255
			oplot,yy,prof_fit_lorentzian,color=150,linestyle=1
			oplot,yy,prof_fit_moffat,color=70,linestyle=2
			legend,['center='+strtrim(coeff[1],2)+','+strtrim(A_moffat[1],2),$
					'FWHM='+strtrim(coeff[2]*2.3548*pix_scale,2)+'"',$
					'slit_loss='+strtrim(slit_loss,2)+'%',$
					'slit_loss(moffat)='+strtrim(slit_loss_moffat,2)+'%'],box=0
		endif
;		plot,x2d,moffat_2d[*,Npts/2.],xr=[coeff[1]-10,coeff[1]+10]
;		oplot,x2d,moffat_2d[*,Npts/4.],color=255
;		oplot,x2d,moffat_2d[*,0],color=70
;		plot,x2d,moffat_2d[*,Npts/2.]/max(moffat_2d[*,Npts/2.]);,xr=[coeff[1]-10,coeff[1]+10]
;		oplot,x2d,moffat_2d[*,Npts/4.]/max(moffat_2d[*,Npts/4.]),color=255
;		oplot,x2d,moffat_2d[*,0]/max(moffat_2d[*,0]),color=70
	
;		levels = [10,100,500,(indgen(10)+1.)/10.*max(moffat_2d)]
;		contour,moffat_2d,x2d,y2d,/isotropic,levels=levels,$
;				xr=[coeff[1]-10,coeff[1]+10],yr=[coeff[1]-10,coeff[1]+10]
;		hline,y2d[Npts/4.],linestyle=1
;		legend,[strtrim(max(moffat_2d),2)]
	endfor
	plot,x_all,slit_loss_all,psym=7,yr=[0,100],$
		xtit='x (pix)',ytit='slit loss (%)'
	median_slit_loss_all = median(slit_loss_all,/even)
	mean_slit_loss_all = mean(slit_loss_all)
	hline,median_slit_loss_all,linestyle=2
	hline,mean_slit_loss_all,linestyle=1
	legend,['median slit loss='+string(median_slit_loss_all,f='(f5.2)')+'%',$
			'mean slit loss='+string(mean_slit_loss_all,f='(f5.2)')+'%'],$
			box=0
	printf,1,strmid(flist[i],49,99),median_slit_loss_all
endfor

device,/close
close,1
	
skip:
device,file='oned_spec_std_compare.ps',/color,xs=15,ys=20,$
		xoffset=0.5,yoffset=0.5
dir_2014dec29 = '/Users/jhyoon/Tools/MosfireDRP-1.0/REDUX/LONGSLIT-46x1.2/2014dec29/Y/'
chr_rdtbl,dir_2014dec29+'fname_std_iraf.list',0,flist_iraf
for i=0,N_elements(flist_iraf)-1 do begin
	print,i,flist_iraf[i]
	f_iraf = dir_2014dec29+flist_iraf[i]+'.fits'
	f_weight_iraf = f_iraf;dir_2014dec29+flist_iraf[i]+'_weight.fits'
	f_optimal = strmid(flist[i],88,30)+'_'+strmid(flist[i],57,9)+'_1dspec.dat'
	flux_iraf = mrdfits(f_iraf,0,hdr_iraf,/silent)
	w0_iraf = sxpar(hdr_iraf,'CRVAL1')
	dw_iraf = sxpar(hdr_iraf,'CD1_1')
	lambda_iraf = w0_iraf + indgen(N_elements(flux_iraf))*dw_iraf	
	flux_iraf_weight = mrdfits(f_weight_iraf,0,hdr_iraf_weight,/silent)
	w0_iraf_weight = sxpar(hdr_iraf_weight,'CRVAL1')
	dw_iraf_weight = sxpar(hdr_iraf_weight,'CD1_1')
	lambda_iraf_weight = w0_iraf_weight + $
						indgen(N_elements(flux_iraf_weight))*dw_iraf_weight
	readcol,f_optimal,lambda_optimal,flux_optimal,f='(f,f)',/silent

	multiplot,[1,3],ygap=0.02,/doxaxis
	if i lt 8 then $
		plot,lambda_optimal,flux_optimal,psym=10,tit=flist_iraf[i]
	if i ge 8 then $
		plot,lambda_optimal,flux_optimal,psym=10,tit=flist_iraf[i],$
			yr=[0,200000]
	oplot,lambda_iraf,flux_iraf,psym=10,color=255
	oplot,lambda_iraf_weight,flux_iraf_weight,psym=10,color=70
	legend,['optimal','apsum','apsum,var'],linestyle=0,$
			color=[0,255,70],box=0,/bottom

	multiplot,/doxaxis
	plot,lambda_optimal,flux_optimal,psym=10,$
		xr=[9970,10200]
	oplot,lambda_iraf,flux_iraf,psym=10,color=255
	oplot,lambda_iraf_weight,flux_iraf_weight,psym=10,color=70

	multiplot,/doxaxis
	plot,lambda_optimal,flux_optimal,psym=10,$
		xr=[10800,11030]
	oplot,lambda_iraf,flux_iraf,psym=10,color=255
	oplot,lambda_iraf_weight,flux_iraf_weight,psym=10,color=70
	multiplot,/reset
	erase

	;;; print out 1d spectra
	forprint,lambda_iraf,flux_iraf,$
				textout=flist_iraf[i]+'_2014dec29_iraf.dat',/nocomment
	forprint,lambda_iraf_weight,flux_iraf_weight,$
				textout=flist_iraf[i]+'_2014dec29_iraf_weight.dat',/nocomment
endfor

;;;
;;; std on 20141230
;;;
dir_2014dec30 = '/Users/jhyoon/Tools/MosfireDRP-1.0/REDUX/LONGSLIT-46x1.2/2014dec30/Y/'
chr_rdtbl,dir_2014dec30+'fname_std_iraf.list',0,flist_iraf
for i=0,N_elements(flist_iraf)-1 do begin
	print,'2014dec30',i,'-',flist_iraf[i]
	f_iraf = dir_2014dec30+flist_iraf[i]+'.fits'
	f_weight_iraf = f_iraf;dir_2014dec30+flist_iraf[i]+'_weight.fits'
	f_optimal = strmid(flist[i],88,30)+'_'+strmid(flist[i],57,9)+'_1dspec.dat'
	flux_iraf = mrdfits(f_iraf,0,hdr_iraf,/silent)
	w0_iraf = sxpar(hdr_iraf,'CRVAL1')
	dw_iraf = sxpar(hdr_iraf,'CD1_1')
	lambda_iraf = w0_iraf + indgen(N_elements(flux_iraf))*dw_iraf	
	flux_iraf_weight = mrdfits(f_weight_iraf,0,hdr_iraf_weight,/silent)
	w0_iraf_weight = sxpar(hdr_iraf_weight,'CRVAL1')
	dw_iraf_weight = sxpar(hdr_iraf_weight,'CD1_1')
	lambda_iraf_weight = w0_iraf_weight + $
						indgen(N_elements(flux_iraf_weight))*dw_iraf_weight
	readcol,f_optimal,lambda_optimal,flux_optimal,f='(f,f)',/silent

	multiplot,[1,3],ygap=0.02,/doxaxis
	if i lt 8 then $
		plot,lambda_optimal,flux_optimal,psym=10,tit=flist_iraf[i]
	if i ge 8 then $
		plot,lambda_optimal,flux_optimal,psym=10,tit=flist_iraf[i],$
			yr=[0,200000]
	oplot,lambda_iraf,flux_iraf,psym=10,color=255
	oplot,lambda_iraf_weight,flux_iraf_weight,psym=10,color=70
	legend,['optimal','apsum','apsum,var'],linestyle=0,$
			color=[0,255,70],box=0,/bottom

	multiplot,/doxaxis
	plot,lambda_optimal,flux_optimal,psym=10,$
		xr=[9970,10200]
	oplot,lambda_iraf,flux_iraf,psym=10,color=255
	oplot,lambda_iraf_weight,flux_iraf_weight,psym=10,color=70

	multiplot,/doxaxis
	plot,lambda_optimal,flux_optimal,psym=10,$
		xr=[10800,11030]
	oplot,lambda_iraf,flux_iraf,psym=10,color=255
	oplot,lambda_iraf_weight,flux_iraf_weight,psym=10,color=70
	multiplot,/reset
	erase

	;;; print out 1d spectra
	forprint,lambda_iraf,flux_iraf,$
				textout=flist_iraf[i]+'_2014dec30_iraf.dat',/nocomment
	forprint,lambda_iraf_weight,flux_iraf_weight,$
				textout=flist_iraf[i]+'_2014dec30_iraf_weight.dat',/nocomment
endfor

;;;
;;; std on 20130104
;;;
dir_2013jan04 = '/Users/jhyoon/Tools/MosfireDRP-1.0/REDUX/LONGSLIT-3x1/2013jan4/Y/'
chr_rdtbl,dir_2013jan04+'fname_std_iraf.list',0,flist_iraf
for i=0,N_elements(flist_iraf)-1 do begin
	print,'2013jan04',i,'-',flist_iraf[i]
	f_iraf = dir_2013jan04+flist_iraf[i]+'.fits'
	f_weight_iraf = f_iraf;dir_2013jan04+flist_iraf[i]+'_weight.fits'
	flux_iraf = mrdfits(f_iraf,0,hdr_iraf,/silent)
	w0_iraf = sxpar(hdr_iraf,'CRVAL1')
	dw_iraf = sxpar(hdr_iraf,'CD1_1')
	lambda_iraf = w0_iraf + indgen(N_elements(flux_iraf))*dw_iraf	
	flux_iraf_weight = mrdfits(f_weight_iraf,0,hdr_iraf_weight,/silent)
	w0_iraf_weight = sxpar(hdr_iraf_weight,'CRVAL1')
	dw_iraf_weight = sxpar(hdr_iraf_weight,'CD1_1')
	lambda_iraf_weight = w0_iraf_weight + $
						indgen(N_elements(flux_iraf_weight))*dw_iraf_weight

	multiplot,[1,3],ygap=0.02,/doxaxis
	if i lt 8 then $
		plot,lambda_optimal,flux_optimal,psym=10,tit=flist_iraf[i]
	if i ge 8 then $
		plot,lambda_optimal,flux_optimal,psym=10,tit=flist_iraf[i],$
			yr=[0,200000]
	oplot,lambda_iraf,flux_iraf,psym=10,color=255
	oplot,lambda_iraf_weight,flux_iraf_weight,psym=10,color=70
	legend,['optimal','apsum','apsum,var'],linestyle=0,$
			color=[0,255,70],box=0,/bottom

	multiplot,/doxaxis
	plot,lambda_optimal,flux_optimal,psym=10,$
		xr=[9970,10200]
	oplot,lambda_iraf,flux_iraf,psym=10,color=255
	oplot,lambda_iraf_weight,flux_iraf_weight,psym=10,color=70

	multiplot,/doxaxis
	plot,lambda_optimal,flux_optimal,psym=10,$
		xr=[10800,11030]
	oplot,lambda_iraf,flux_iraf,psym=10,color=255
	oplot,lambda_iraf_weight,flux_iraf_weight,psym=10,color=70
	multiplot,/reset
	erase

	;;; print out 1d spectra
	forprint,lambda_iraf,flux_iraf,$
				textout=flist_iraf[i]+'_2013jan04_iraf.dat',/nocomment
	forprint,lambda_iraf_weight,flux_iraf_weight,$
				textout=flist_iraf[i]+'_2013jan04_iraf_weight.dat',/nocomment
endfor



device,/close


stop
end
