pro oned_spec_std

;dir = '/Users/jhyoon/Tools/MosfireDRP-1.0/REDUX/mosFv2/2013jan4/Y/'
chr_rdtbl,'fname_std.list',0,flist
chr_rdtbl,'fname_std_var.list',0,flist_var
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
interval = 10

goto, skip

openw,1,'slit_loss_star2.dat'
ypix  = [71,32,20,49,60,47,55,40,71]
dy = [15,15,15,15,15,15,15,15,9]
;slit_center = 0.5-[1.05/3.275,1.0/3.275,1.4/3.275,1.4/3.275,1.2/3.275,1.4/3.275,2.95/6.65]

t_obs_star = []
fwhm_star = []
for i=0,N_star-1 do begin
	print,i,flist_star[i]
	f_eps = flist_star[i]+'_eps.fits'
	f_sig = flist_star[i]+'_sig.fits'
	f_out = strmid(flist_star[i],61,22)+'_'+strmid(flist_star[i],49,9)+'_1dspec.dat'
	;;;
	;;; psf fit and compute slit loss
	;;;
	;spec2d = mrdfits(f_eps,0,hdr)
	fits_read,f_eps,spec2d,hdr
	fits_read,f_sig,spec2d_sig,hdr_sig
	f_eps = f_eps[*,0,0]
	f_sig = f_sig[*,0,0]

	yy = indgen(dy[i]*2+1)
	pix_scale = sxpar(hdr,'PSCALE') ; arcsec	
	t_obs_str = sxpar(hdr,'UTC')
	t_obs = float(strmid(t_obs_str,0,2)) + $
			float(strmid(t_obs_str,3,2))/60. + $
			float(strmid(t_obs_str,6,5))/3600.	; UTC

	x_all = []
	slit_loss_all = []
	peak_all = []
	fwhm_all = []
	x_gaussfit = []
	mean_gaussfit = []
	sigma_gaussfit = []
	for j=0+10,N_elements(spec2d[*,0])-10,interval do begin
		;prof = spec2d[j,ypix[i]-dy[i]:ypix[i]+dy[i]]
		;prof_sig = spec2d_sig[j,ypix[i]-dy[i]:ypix[i]+dy[i]]
		prof = total(spec2d[j:j+interval-1,ypix[i]-dy[i]:ypix[i]+dy[i]],1)
		prof_sig = total(spec2d_sig[j:j+interval-1,ypix[i]-dy[i]:ypix[i]+dy[i]],1)
		;;; fit psf
		prof_fit = gaussfit(yy,prof,coeff,nterms=3)
		prof_fit_lorentzian = mpfitpeak(yy,prof,A_lorentzian,/lorentzian)
		prof_fit_moffat = mpfitpeak(yy,prof,A_moffat,/moffat)

		x_gaussfit = [x_gaussfit,j]
		mean_gaussfit = [mean_gaussfit,coeff[1]]
		sigma_gaussfit= [sigma_gaussfit,coeff[2]]
		;;; compute 1D slit loss gaussian and moffat psf
		;slit_left = dy[i]-slit_center[i]/pix_scale-0.5/pix_scale
		;slit_left2 = dy[i]-slit_center[i]/pix_scale-3*0.5/pix_scale
		;slit_right = dy[i]-slit_center[i]/pix_scale+0.5/pix_scale
		;slit_right2 = dy[i]-slit_center[i]/pix_scale+3*0.5/pix_scale
		slit_left = dy[i]-0.5/pix_scale
		slit_left2 = dy[i]-3*0.5/pix_scale
		slit_right = dy[i]+0.5/pix_scale
		slit_right2 = dy[i]+3*0.5/pix_scale

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
		sn_in_slit = total(prof[sub_in_slit])/sqrt(total(prof_sig[sub_in_slit]^2.))
		fwhm = coeff[2]*2.3548*pix_scale	; in arcsec
		fwhm_all = [fwhm_all,fwhm]
;		if abs(A_moffat[1]-dy[i]) lt 3 and $
;			A_moffat[0] gt 0 and $
;			sn_in_slit gt 5 and $
;			fwhm lt 1.5 and $
;			fwhm gt 0.4 then begin

			x_all = [x_all,j]
			slit_loss_all = [slit_loss_all,slit_loss_moffat]
			peak_all = [peak_all,A_moffat[1]]

			;;; plot
;			plot,yy,prof,psym=10,xr=[dy[i]-dy[i],dy[i]+dy[i]],$
;				tit=strmid(flist_star[i],61,25)+',x='+strtrim(j,2)
;			oplot,yy,prof_sig,color=220,psym=10
;			oplot,yy,prof_fit,color=255
;			oplot,yy,prof_fit_lorentzian,color=150,linestyle=1
;			oplot,yy,prof_fit_moffat,color=70,linestyle=2
;			legend,['center='+strtrim(coeff[1],2)+','+strtrim(A_moffat[1],2),$
;					'FWHM='+strtrim(fwhm,2)+'"',$
;					'slit_loss='+strtrim(slit_loss,2)+'%',$
;					'slit_loss(moffat)='+strtrim(slit_loss_moffat,2)+'%'],box=0
;		endif
;		vline,dy[i],linestyle=2
;		vline,slit_left,linestyle=1
;		vline,slit_right,linestyle=1
	endfor
	t_obs_star = [t_obs_star,t_obs]
	fwhm_star = [fwhm_star,median(fwhm_all)]

	median_slit_loss_all = median(slit_loss_all,/even)
	mean_slit_loss_all = mean(slit_loss_all)
	median_peak_all = median(peak_all,/even)
	mean_peak_all = mean(peak_all)
;	;;; plot mean slit loss
;	plot,x_all,slit_loss_all,psym=7,yr=[0,100],$
;		xtit='x (pix)',ytit='slit loss (%)'
;	hline,median_slit_loss_all,linestyle=2,color=100
;	hline,mean_slit_loss_all,linestyle=1,color=150
;	legend,['median slit loss='+ $
;			string(median_slit_loss_all,f='(f5.2)')+'%',$
;			'mean slit loss='+ $
;			string(mean_slit_loss_all,f='(f5.2)')+'%'],box=0
;	;;; plot mean center
;	plot,x_all,peak_all,psym=7,yr=[dy[i]-10,dy[i]+10],$
;		xtit='x (pix)',ytit='moffat center position'
;	fit = ladfit(x_all,peak_all)
;	yfit = fit[0]+fit[1]*x_all
;	oplot,x_all,yfit,color=255
;	hline,median_peak_all,linestyle=2,color=100
;	hline,mean_peak_all,linestyle=1,color=150
;	legend,['median center='+string(median_peak_all,f='(f5.2)'),$
;			'mean center='+string(mean_peak_all,f='(f5.2)')],$
;			box=0


	;;;
	;;; find and fit curvature of a 2D spectrum
	;;;
	ploterror,x_gaussfit,mean_gaussfit,sigma_gaussfit,psym=4,$
				xr=[0,1400],yr=[dy[i]-5,dy[i]+5],$
				tit=strmid(flist_star[i],61,25),xtit='x (pix)',$
				ytit='mu of gaussfit'
	dx = 100
	xx = findgen(1200)/dx
	count = 0
	for kkk=0,100 do begin
		N_mean = N_elements(mean_gaussfit)	
		;fit_legendre = svdfit(x_gaussfit,mean_gaussfit,3,$
		;							func='flegendre',yfit=yfit)
		;fit = poly_fit(x_gaussfit,mean_gaussfit,3,yfit=yfit,yerror=yerror)
		fit = ladfit(x_gaussfit,mean_gaussfit,absdev=yerror)
		yfit = fit[0]+fit[1]*x_gaussfit
		sub_inlier = where(abs(yfit-mean_gaussfit) lt 3*yerror,count)
		if N_mean eq count then break
		x_gaussfit = x_gaussfit[sub_inlier]
		mean_gaussfit = mean_gaussfit[sub_inlier]
		print, kkk,yerror
	endfor
	oplot,x_gaussfit,yfit,color=255

	;;;
	;;; coadd a spectrum over wavelength by correcting the curvature
	;;;
	prof_coadd = 0
	prof_var_coadd = 0	
	for j=0+10,N_elements(spec2d[*,0])-10 do begin
		prof = spec2d[j,ypix[i]-dy[i]:ypix[i]+dy[i]]
		prof_sig = spec2d_sig[j,ypix[i]-dy[i]:ypix[i]+dy[i]]
		;ycen = fit[0]+fit[1]*j+fit[2]*j^2.+fit[3]*j^3.
		;ycen = fit[0]+fit[1]*j+fit[2]*j^2.
		ycen = fit[0]+fit[1]*j
		prof_shifted = shift(prof,round(dy[i]-ycen))
		prof_sig_shifted = shift(prof_sig,round(dy[i]-ycen))
		prof_coadd += prof_shifted
		prof_var_coadd += (prof_sig_shifted)^2.
	endfor

	;;;
	;;; plot coadded profile
	;;;
	plot,yy,prof_coadd,psym=10,xtit='y',ytit='f',$
		tit=strmid(flist_star[i],61,25)+' coadded'
	prof_fit = gaussfit(yy,prof_coadd,coeff,nterms=3)
	prof_fit_moffat = mpfitpeak(yy,prof_coadd,A_moffat,/moffat)
	fwhm = coeff[2]*2.3548*pix_scale	; in arcsec
	oplot,yy,prof_fit,color=255
	oplot,yy,prof_fit_moffat,color=70,linestyle=2
	flux_in_slit_moffat = total(prof_fit_moffat[sub_in_slit])
	flux_out_slit_moffat = total(prof_fit_moffat[sub_out_slit])
	slit_loss_moffat = flux_out_slit_moffat / $
							(flux_in_slit_moffat+flux_out_slit_moffat)*100
	vline,dy[i],linestyle=2
	vline,slit_left,linestyle=1
	vline,slit_right,linestyle=1
	legend,['center='+strtrim(coeff[1],2)+','+strtrim(A_moffat[1],2),$
			'mean FWHM='+strtrim(sqrt(mean(fwhm_star^2.)),2)+'"',$
			'FWHM='+strtrim(fwhm,2)+'"',$
			'slit_loss='+strtrim(slit_loss,2)+'%',$
			'slit_loss(moffat)='+strtrim(slit_loss_moffat,2)+'%'],box=0

	printf,1,strmid(flist_star[i],49,99),median_slit_loss_all,$
			slit_loss_moffat,fwhm,f='(a35,f8.3,f8.3,f7.3)'

;	xx = findgen(N_elements(spec2d[*,0]))
;	ycen_all = ypix[i]-dy[i]+(fit[0]+fit[1]*xx)
;	horne_extract_all,f_eps,f_sig,f_out,ypix=round(ycen_all),dy=dy[i],/gauss,sigma_gauss=coeff[1]

endfor
close,1

dir = '/Users/jhyoon/Tools/MosfireDRP-1.0/REDUX/'
flist_1deps = ['mosR1v4/2014dec29/Y/mosR1v4_Y_star_1012740_eps_1dspec.fits',$
			'mosR1v4/2014dec29/Y/mosR1v4_Y_star_1013948_eps_1dspec.fits',$
			'mosR2v4/2014dec29/Y/mosR2v4_Y_star_1014109_eps_1dspec.fits',$
			'mosR2v4/2014dec29/Y/mosR2v4_Y_star_1036896_eps_1dspec.fits',$
			'mosR2v4/2014dec30/Y/mosR2v4_Y_star_1036896_eps_1dspec.fits',$
			'mosR2v4/2014dec30/Y/mosR2v4_Y_star_1014109_eps_1dspec.fits',$
			'mosR3v4/2014dec30/Y/mosR3v4_Y_star_771147_eps_1dspec.fits',$
			'mosR3v4/2014dec30/Y/mosR3v4_Y_star_776669_eps_1dspec.fits']
flist_1dvar = ['mosR1v4/2014dec29/Y/mosR1v4_Y_star_1012740_var_1dspec.fits',$
			'mosR1v4/2014dec29/Y/mosR1v4_Y_star_1013948_var_1dspec.fits',$
			'mosR2v4/2014dec29/Y/mosR2v4_Y_star_1014109_var_1dspec.fits',$
			'mosR2v4/2014dec29/Y/mosR2v4_Y_star_1036896_var_1dspec.fits',$
			'mosR2v4/2014dec30/Y/mosR2v4_Y_star_1036896_var_1dspec.fits',$
			'mosR2v4/2014dec30/Y/mosR2v4_Y_star_1014109_var_1dspec.fits',$
			'mosR3v4/2014dec30/Y/mosR3v4_Y_star_771147_var_1dspec.fits',$
			'mosR3v4/2014dec30/Y/mosR3v4_Y_star_776669_var_1dspec.fits']
flist_eps_out = ['mosR1v4_Y_star_1012740_2014dec29_1dspec.dat',$
			'mosR1v4_Y_star_1013948_2014dec29_1dspec.dat',$
			'mosR2v4_Y_star_1014109_2014dec29_1dspec.dat',$
			'mosR2v4_Y_star_1036896_2014dec29_1dspec.dat',$
			'mosR2v4_Y_star_1036896_2014dec30_1dspec.dat',$
			'mosR2v4_Y_star_1014109_2014dec30_1dspec.dat',$
			'mosR3v4_Y_star_771147_2014dec30_1dspec.dat',$
			'mosR3v4_Y_star_776669_2014dec30_1dspec.dat']
for i=0,N_elements(flist_1deps)-1 do begin
	fits_read,dir+flist_1deps[i],spec,hdr
	fits_read,dir+flist_1dvar[i],var,hdr2
	spec = spec[*,0,0]
	var = spec[*,0,0]
	w0 = sxpar(hdr,'CRVAL1')
	dw = 1.08597285067845;sxpar(hdr,'CD2_2')
	w = w0+abs(dw)*findgen(N_elements(spec))
	sub_crazy = where(finite(sqrt(var)) ne 1,count,complement=sub_normal)
	if count gt 0 then var[sub_crazy] = mean(var[sub_normal])
	forprint,w,spec,sqrt(var),textout=flist_eps_out[i],/nocomment
endfor	

;device,/close
;stop


plot,t_obs_star[0:2],fwhm_star[0:2],psym=-7,xr=[10,14],yr=[0.3,0.8],$
	xtit='t_obs (UTC)',ytit='FWHM (arcsec)'
oplot,t_obs_star[3:5],fwhm_star[3:5],psym=-1,color=70
oplot,[t_obs_star[6]],[fwhm_star[6]],psym=-5,color=150


skip: print,'skip'
;;;
;;; standard stars
;;;
ypix = [1016,1061,1006,1071,1005,1071,1004,1071,1005,1071,$

		1016,1060,1016,1062,1017,1061,1016,1060,1017,1062,$
		1011,1065,1011,1066,1010,1065,1011,1065,1010,1066,1010,1065,1010,1065,1010,1065,1010,1065,1011,1065,1010,1065]
dy = [20,20,25,25,25,25,35,30,25,25,$
		20,20,20,20,20,20,20,20,15,20,$
		25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25]

openw,1,'slit_loss_std2.dat'
for i=0,N_file-1 do begin
	print,i,'-',flist[i]
	f_eps = flist[i]+'.fits.gz'
	f_sig = flist_var[i]+'.fits.gz'
	;f_out = strmid(f_eps,69,50)+'_1dspec.dat'
	f_out = strmid(flist[i],88,30)+'_'+strmid(flist[i],57,9)+'_1dspec.dat'
	if i ge 20 then f_out = strmid(flist[i],84,34)+'_'+strmid(flist[i],54,8)+'_1dspec.dat'
	;;;
	;;; Horne extraction of 1D spectra
	;;;
;	horne_extract_all,f_eps,f_sig,f_out,$
;					ypix=ypix[i],dy=dy[i]

;	print,'####### f_out ',f_out
;	readcol,f_out,lambda,flux,flux_err,flux2,/silent
	
	;plot,lambda,flux,psym=10,xr=[9750,11200],$
	;	xtit=textoidl('\lambda (\AA)'),ytit='flux',$
	;	tit=strmid(flist[i],69,50)
	;oplot,lambda,flux_err,color=255,psym=10
	;;;
	;;; psf fit and compute slit loss
	;;;
	spec2d = mrdfits(f_eps,0,hdr)
	yy = indgen(dy[i]*2+1)
	pix_scale = sxpar(hdr,'PSCALE') ; arcsec	
	x_all = []
	slit_loss_all = []
	peak_all = []
	x_gaussfit = []
	mean_gaussfit = []
	sigma_gaussfit = []
	slit_width = 1.2
	if i ge 20 then slit_width = 1.0

	for j=0+10,N_elements(spec2d[*,0])-10,interval do begin
		prof = spec2d[j,ypix[i]-dy[i]:ypix[i]+dy[i]]
		;;; fit psf
		prof_fit = gaussfit(yy,prof,coeff,nterms=3)
		prof_fit_lorentzian = mpfitpeak(yy,prof,A_lorentzian,/lorentzian)
		prof_fit_moffat = mpfitpeak(yy,prof,A_moffat,/moffat)

		;;; compute 1D slit loss gaussian and moffat psf
		x_gaussfit = [x_gaussfit,j]
		mean_gaussfit = [mean_gaussfit,coeff[1]]
		sigma_gaussfit= [sigma_gaussfit,coeff[2]]
		slit_center = dy[i];coeff[1]
		slit_left = slit_center-slit_width/2./pix_scale
		slit_left2 = slit_center-5*slit_width/2./pix_scale
		slit_right = slit_center+slit_width/2./pix_scale
		slit_right2 = slit_center+5*slit_width/2./pix_scale
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
		peak_all = [peak_all,A_moffat[1]]
	
		f_tit = strmid(flist[i],88,25)+',x='+strtrim(j,2)
		if j/200 eq j/200. then begin
			plot,yy,prof,psym=10,xr=[coeff[1]-10,coeff[1]+10],tit=f_tit
			oplot,yy,prof_fit,color=255
			oplot,yy,prof_fit_lorentzian,color=150,linestyle=1
			oplot,yy,prof_fit_moffat,color=70,linestyle=2
			vline,slit_center,linestyle=2
			vline,slit_left,linestyle=1
			vline,slit_right,linestyle=1
			legend,['center='+strtrim(coeff[1],2)+','+strtrim(A_moffat[1],2),$
					'FWHM='+strtrim(coeff[2]*2.3548*pix_scale,2)+'"',$
					'slit_loss='+strtrim(slit_loss,2)+'%',$
					'slit_loss(moffat)='+strtrim(slit_loss_moffat,2)+'%'],box=0
		endif
	endfor

	;;; plot mean slit loss
	plot,x_all,slit_loss_all,psym=7,yr=[0,100],$
		xtit='x (pix)',ytit='slit loss (%)'
	median_slit_loss_all = median(slit_loss_all,/even)
	mean_slit_loss_all = mean(slit_loss_all)
	std_slit_loss_all = stddev(slit_loss_all)
	hline,median_slit_loss_all,linestyle=2
	hline,mean_slit_loss_all,linestyle=1
	legend,['median slit loss='+string(median_slit_loss_all,f='(f5.2)')+'%',$
			'mean slit loss='+string(mean_slit_loss_all,f='(f5.2)')+'%',$
			'std slit loss='+string(std_slit_loss_all,f='(f5.2)')+'%'],$
			box=0


	;;; plot mean center
	plot,x_all,peak_all,psym=7,yr=[dy[i]-5,dy[i]+5],$
		xtit='x (pix)',ytit='moffat center position'
	fit = ladfit(x_all,peak_all)
	yfit = fit[0]+fit[1]*x_all
	oplot,x_all,yfit,color=255
	median_peak_all = median(peak_all,/even)
	mean_peak_all = mean(peak_all)
	hline,median_peak_all,linestyle=2,color=100
	hline,mean_peak_all,linestyle=1,color=150
	legend,['median center='+string(median_peak_all,f='(f5.2)'),$
			'mean center='+string(mean_peak_all,f='(f5.2)')],$
			box=0


	;;;
	;;; find and fit curvature of a 2D spectrum
	;;;
	f_tit = strmid(flist[i],88,25)
	ploterror,x_gaussfit,mean_gaussfit,sigma_gaussfit,psym=4,$
				xr=[0,1400],yr=[dy[i]-5,dy[i]+5],$
				tit=f_tit,xtit='x (pix)',$
				ytit='mu of gaussfit'
	dx = 100
	xx = findgen(1200)/dx
	count = 0
	for kkk=0,100 do begin
		N_mean = N_elements(mean_gaussfit)	
		;fit_legendre = svdfit(x_gaussfit,mean_gaussfit,3,$
		;							func='flegendre',yfit=yfit)
		;fit = poly_fit(x_gaussfit,mean_gaussfit,3,yfit=yfit,yerror=yerror)
		fit = ladfit(x_gaussfit,mean_gaussfit,absdev=yerror)
		yfit = fit[0]+fit[1]*x_gaussfit
		sub_inlier = where(abs(yfit-mean_gaussfit) lt 3*yerror,count)
		if N_mean eq count then break
		x_gaussfit = x_gaussfit[sub_inlier]
		mean_gaussfit = mean_gaussfit[sub_inlier]
		print, kkk,yerror
	endfor
	oplot,x_gaussfit,yfit,color=255

	;;;
	;;; coadd a spectrum over wavelength by correcting the curvature
	;;;
	prof_coadd = 0
	prof_var_coadd = 0	
	for j=0+10,N_elements(spec2d[*,0])-10 do begin
		prof = spec2d[j,ypix[i]-dy[i]:ypix[i]+dy[i]]
		;prof_sig = spec2d_sig[j,ypix[i]-dy[i]:ypix[i]+dy[i]]
		;ycen = fit[0]+fit[1]*j+fit[2]*j^2.+fit[3]*j^3.
		;ycen = fit[0]+fit[1]*j+fit[2]*j^2.
		ycen = fit[0]+fit[1]*j
		prof_shifted = shift(prof,round(dy[i]-ycen))
		;prof_sig_shifted = shift(prof_sig,round(dy[i]-ycen))
		prof_coadd += prof_shifted
		;prof_var_coadd += (prof_sig_shifted)^2.
	endfor

	;;;
	;;; plot coadded profile
	;;;
	plot,yy,prof_coadd,psym=10,xtit='y',ytit='f',$
		tit=f_tit+' coadded',xr=[dy[i]-15,dy[i]+15]
	prof_fit = gaussfit(yy,prof_coadd,coeff,nterms=3)
	prof_fit_moffat = mpfitpeak(yy,prof_coadd,A_moffat,/moffat)
	fwhm = coeff[2]*2.3548*pix_scale	; in arcsec
	oplot,yy,prof_fit,color=255
	oplot,yy,prof_fit_moffat,color=70,linestyle=2
	flux_in_slit_moffat = total(prof_fit_moffat[sub_in_slit])
	flux_out_slit_moffat = total(prof_fit_moffat[sub_out_slit])
	slit_loss_moffat = flux_out_slit_moffat / $
							(flux_in_slit_moffat+flux_out_slit_moffat)*100
	vline,dy[i],linestyle=2
	vline,slit_left,linestyle=1
	vline,slit_right,linestyle=1
	legend,['center='+strtrim(coeff[1],2)+','+strtrim(A_moffat[1],2),$
			'mean FWHM='+strtrim(sqrt(mean(fwhm^2.)),2)+'"',$
			'FWHM='+strtrim(fwhm,2)+'"',$
			'slit_loss='+strtrim(slit_loss,2)+'%',$
			'slit_loss(moffat)='+strtrim(slit_loss_moffat,2)+'%'],box=0

	slit_loss_avg = 0
	if i mod 2 eq 0 then slit_loss_pre = slit_loss_moffat
	if i mod 2 eq 1 then slit_loss_avg = (slit_loss_pre+slit_loss_moffat)/2.

	f_tit = strmid(flist[i],57,9)+'__'+strmid(flist[i],88,29)
	if i ge 20 then f_tit = strmid(flist[i],54,8)+'__'+strmid(flist[i],84,41)
	printf,1,f_tit,median_slit_loss_all,$
			slit_loss_moffat,slit_loss_avg,fwhm,f='(a41,3(f8.3),f7.3)'

endfor


device,/close
close,1
	

;skip:
device,file='oned_spec_std_compare.ps',/color,xs=15,ys=20,$
		xoffset=0.5,yoffset=0.5
dir_2014dec29 = '/Users/jhyoon/Tools/MosfireDRP-1.0/REDUX/LONGSLIT-46x1.2/2014dec29/Y/'
chr_rdtbl,dir_2014dec29+'fname_std_iraf.list',0,flist_iraf
for i=0,N_elements(flist_iraf)-1 do begin
	print,i,flist_iraf[i]
	f_iraf = dir_2014dec29+flist_iraf[i]+'.fits'
	f_optimal = strmid(flist[i],88,30)+'_'+strmid(flist[i],57,9)+'_1dspec.dat'
	print,'------f_optimal ',f_optimal

	flux_iraf = mrdfits(f_iraf,0,hdr_iraf,/silent)
	flux_iraf = flux_iraf[*,0,0]
	w0_iraf = sxpar(hdr_iraf,'CRVAL1')
	dw_iraf = sxpar(hdr_iraf,'CD1_1')
	lambda_iraf = w0_iraf + indgen(N_elements(flux_iraf))*dw_iraf	
	readcol,f_optimal,lambda_optimal,flux_optimal,f='(f,f)',/silent

	multiplot,[1,3],ygap=0.02,/doxaxis
	if i lt 8 then $
		plot,lambda_optimal,flux_optimal,psym=10,tit=flist_iraf[i]
	if i ge 8 then $
		plot,lambda_optimal,flux_optimal,psym=10,tit=flist_iraf[i],$
			yr=[0,200000]
	oplot,lambda_iraf,flux_iraf,psym=10,color=255
	legend,['optimal','apsum','apsum,var'],linestyle=0,$
			color=[0,255,70],box=0,/bottom

	multiplot,/doxaxis
	plot,lambda_optimal,flux_optimal,psym=10,$
		xr=[9970,10200]
	oplot,lambda_iraf,flux_iraf,psym=10,color=255

	multiplot,/doxaxis
	plot,lambda_optimal,flux_optimal,psym=10,$
		xr=[10800,11030]
	oplot,lambda_iraf,flux_iraf,psym=10,color=255
	multiplot,/reset
	erase

	;;; print out 1d spectra
	forprint,lambda_iraf,flux_iraf,$
				textout=flist_iraf[i]+'_2014dec29_iraf.dat',/nocomment,/silent
endfor

;;;
;;; std on 20141230
;;;
dir_2014dec30 = '/Users/jhyoon/Tools/MosfireDRP-1.0/REDUX/LONGSLIT-46x1.2/2014dec30/Y/'
chr_rdtbl,dir_2014dec30+'fname_std_iraf.list',0,flist_iraf
for i=0,N_elements(flist_iraf)-1 do begin
	print,'2014dec30',i,'-',flist_iraf[i]
	f_iraf = dir_2014dec30+flist_iraf[i]+'.fits'
	f_optimal = strmid(flist[i+10],88,30)+'_'+strmid(flist[i+10],57,9)+'_1dspec.dat'
	print,'------f_optimal ',f_optimal
	flux_iraf = mrdfits(f_iraf,0,hdr_iraf,/silent)
	flux_iraf = flux_iraf[*,0,0]
	w0_iraf = sxpar(hdr_iraf,'CRVAL1')
	dw_iraf = sxpar(hdr_iraf,'CD1_1')
	lambda_iraf = w0_iraf + indgen(N_elements(flux_iraf))*dw_iraf	
	readcol,f_optimal,lambda_optimal,flux_optimal,f='(f,f)',/silent

	multiplot,[1,3],ygap=0.02,/doxaxis
	if i lt 8 then $
		plot,lambda_optimal,flux_optimal,psym=10,tit=flist_iraf[i]
	if i ge 8 then $
		plot,lambda_optimal,flux_optimal,psym=10,tit=flist_iraf[i],$
			yr=[0,200000]
	oplot,lambda_iraf,flux_iraf,psym=10,color=255
	legend,['optimal','apsum','apsum,var'],linestyle=0,$
			color=[0,255,70],box=0,/bottom

	multiplot,/doxaxis
	plot,lambda_optimal,flux_optimal,psym=10,$
		xr=[9970,10200]
	oplot,lambda_iraf,flux_iraf,psym=10,color=255

	multiplot,/doxaxis
	plot,lambda_optimal,flux_optimal,psym=10,$
		xr=[10800,11030]
	oplot,lambda_iraf,flux_iraf,psym=10,color=255
	multiplot,/reset
	erase

	;;; print out 1d spectra
	forprint,lambda_iraf,flux_iraf,$
				textout=flist_iraf[i]+'_2014dec30_iraf.dat',/nocomment,/silent
endfor

;;;
;;; std on 20130104
;;;
dir_2013jan04 = '/Users/jhyoon/Tools/MosfireDRP-1.0/REDUX/LONGSLIT-3x1/2013jan4/Y/'
chr_rdtbl,dir_2013jan04+'fname_std_iraf.list',0,flist_iraf
for i=0,N_elements(flist_iraf)-1 do begin
	print,'2013jan04',i,'-',flist_iraf[i]
	f_iraf = dir_2013jan04+flist_iraf[i]+'.fits'
	flux_iraf = mrdfits(f_iraf,0,hdr_iraf,/silent)
	flux_iraf = flux_iraf[*,0,0]
	w0_iraf = sxpar(hdr_iraf,'CRVAL1')
	dw_iraf = sxpar(hdr_iraf,'CD1_1')
	lambda_iraf = w0_iraf + indgen(N_elements(flux_iraf))*dw_iraf	

	multiplot,[1,3],ygap=0.02,/doxaxis
	if i lt 8 then $
		plot,lambda_optimal,flux_optimal,psym=10,tit=flist_iraf[i]
	if i ge 8 then $
		plot,lambda_optimal,flux_optimal,psym=10,tit=flist_iraf[i],$
			yr=[0,200000]
	oplot,lambda_iraf,flux_iraf,psym=10,color=255
	legend,['optimal','apsum','apsum,var'],linestyle=0,$
			color=[0,255,70],box=0,/bottom

	multiplot,/doxaxis
	plot,lambda_optimal,flux_optimal,psym=10,$
		xr=[9970,10200]
	oplot,lambda_iraf,flux_iraf,psym=10,color=255

	multiplot,/doxaxis
	plot,lambda_optimal,flux_optimal,psym=10,$
		xr=[10800,11030]
	oplot,lambda_iraf,flux_iraf,psym=10,color=255
	multiplot,/reset
	erase

	;;; print out 1d spectra
	forprint,lambda_iraf,flux_iraf,$
				textout=flist_iraf[i]+'_2013jan04_iraf.dat',/nocomment,/silent
endfor



device,/close


stop
end
