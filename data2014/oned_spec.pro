pro oned_spec

;dir = '/Users/jhyoon/Tools/MosfireDRP-1.0/REDUX/mosFv2/2013jan4/Y/'
;chr_rdtbl,'fname.list',0,flist
;N_file = (size(flist))[2]

;flist = ['/Users/jhyoon/Tools/MosfireDRP-1.0/REDUX/mosR1v4/2014dec29/Y/mosR1v4_Y_R1_1009808',$
;			'/Users/jhyoon/Tools/MosfireDRP-1.0/REDUX/mosR1v4/2014dec29/Y/mosR1v4_Y_R1_1015084',$
;			'/Users/jhyoon/Tools/MosfireDRP-1.0/REDUX/mosR1v4/2014dec29/Y/mosR1v4_Y_R2_1015516',$
;			'/Users/jhyoon/Tools/MosfireDRP-1.0/REDUX/mosR2v4/2014dec29/Y/mosR2v4_Y_R2_1015516',$
;			'/Users/jhyoon/Tools/MosfireDRP-1.0/REDUX/mosR2v4/2014dec30/Y/mosR2v4_Y_R2_1015516',$
;			'/Users/jhyoon/Tools/MosfireDRP-1.0/REDUX/mosR3v4/2014dec30/Y/mosR3v4_Y_R2_1015516',$
;			'/Users/jhyoon/Tools/MosfireDRP-1.0/REDUX/mosR2v4/2014dec29/Y/mosR2v4_Y_R2_1010583',$
;			'/Users/jhyoon/Tools/MosfireDRP-1.0/REDUX/mosR2v4/2014dec30/Y/mosR2v4_Y_R2_1010583',$
;			'/Users/jhyoon/Tools/MosfireDRP-1.0/REDUX/mosR2v4/2014dec29/Y/mosR2v4_Y_R2_1010994',$
;			'/Users/jhyoon/Tools/MosfireDRP-1.0/REDUX/mosR2v4/2014dec30/Y/mosR2v4_Y_R2_1010994',$
;			'/Users/jhyoon/Tools/MosfireDRP-1.0/REDUX/mosR2v4/2014dec29/Y/mosR2v4_Y_R2_1013478',$
;			'/Users/jhyoon/Tools/MosfireDRP-1.0/REDUX/mosR2v4/2014dec30/Y/mosR2v4_Y_R2_1013478',$
;			'/Users/jhyoon/Tools/MosfireDRP-1.0/REDUX/mosR3v4/2014dec30/Y/mosR3v4_Y_R3_772773',$
;			'/Users/jhyoon/Tools/MosfireDRP-1.0/REDUX/mosR3v4/2014dec30/Y/mosR3v4_Y_R3_775548']

;plot_list = [1,1,0,0,0,1,0,1,0,1,0,1,1,1]
dir = '/Users/jhyoon/Tools/MosfireDRP-1.0/REDUX/'
readcol,'fname_mosR123.list',flist,redshift,ypix,gr,f='(a,f,i,a)',/silent
flist = dir+flist
N_file = N_elements(flist)
dy = replicate(15.,N_file)

;;; normal aperture width is +/-15 pix. These are different due to contaminants from close companion, etc.
dy[10] = 14
ypix[34] = 41
dy[34] = 11
dy[38] = 12
ypix[40] = 39.5
dy[40] = 11.5
dy[49] = 8
dy[53] = 12
ypix[55] = 37.5 
dy[55] = 11.5

;xpix1 = [785,1163,1169,1169,1169,1170,920,921,896,897,902,902,785,1112]
;xpix2 = [791,1172,1174,1174,1174,1175,931,931,903,901,908,908,792,1120]
;ypix = [32,38,49,47,44,102,64,62,97,95,64,63,57,53]
;xmin = [10000,10400,10420,10420,10420,10420,10150,10150,10130,10130,10140,10140,10000,10360]
;xmax = [11200,11600,11620,11620,11620,11620,11350,11350,11330,11330,11340,11340,11200,11560]

chr_rdtbl,'/Users/jhyoon/Dropbox/Research/mosfire/t_exp4obs/t_exp_Henry2013.dat',1,henry2013
id_henry2013 = henry2013[0,*]
redshift_henry2013 = double(henry2013[3,*]) 

;readcol,'mosR_Ha_redshift.dat',redshift
;openw,1,'mosR_Ha_redshift2.dat'

@plot_setting
device,file='oned_spec.ps',/color,/cmyk
for i=0,N_file-1 do begin
	fname = flist[i]
	print,i,'-----',fname
	if i eq 43 then goto,skip
	f_eps = fname+'_eps.fits'
	f_sig = fname+'_sig.fits'
	f_out = strmid(fname,61,20)+'_'+strmid(fname,49,9)+'_1dspec.dat'
	if strmid(fname,71,4) eq 'star' then begin
		f_out = strmid(fname,61,22)+'_'+strmid(fname,49,9)+'_1dspec.dat'
	endif

	sub_match = where(strmid(fname,74,7) eq id_henry2013)
	if sub_match[0] ne -1 then $	
		Ha_obs_henry2013 = 6563.*(1.+redshift_henry2013[sub_match])[0]

;	NII_obs = 6584.*(1.+redshift_henry2013[sub_match])[0]
	Ha_obs = 6563.*(1.+abs(redshift[i]))[0]
	NII_obs = 6584.*(1.+abs(redshift[i]))[0]
	lamguess_ha = Ha_obs	; guess observed Ha wavelength

	print,ypix[i],dy[i]
	if i eq 10 then flux2 = flux
	horne_extract,f_eps,f_sig,f_out,ypix=ypix[i],dy_input=dy[i],$
;					ypix=ypix[i],xpix=[xpix1[i],xpix2[i]],$
					lamguess_ha=lamguess_ha,/gauss
	vline,Ha_obs,linestyle=2
	vline,NII_obs,linestyle=1
	
	readcol,f_out,lambda,flux,flux_err,/silent
	sub_Ha = where(abs(Ha_obs-lambda) le 7)
	max_flux = max(flux[sub_Ha],max_sub)
	print,lambda[sub_Ha[max_sub]]
;	printf,1,(lambda[sub_Ha[max_sub]]-6563)/6563
	if i eq 10 then flux = flux2
	
	plot,lambda,flux,psym=10,xr=[9750,11200],$
		xtit=textoidl('\lambda (\AA)'),ytit='flux',tit=strmid(fname,74,10)
	oplot,lambda,flux_err,color=255,psym=10
	vline,Ha_obs,linestyle=2,color=70
	vline,NII_obs,linestyle=1,color=150
	
	plot,lambda,flux,psym=10,xr=[Ha_obs-50,NII_obs+50],$
		xtit=textoidl('\lambda (\AA)'),ytit='flux',tit=strmid(fname,74,10)
	oplot,lambda,flux_err,psym=10,color=255
	vline,Ha_obs,linestyle=2,color=70
	vline,NII_obs,linestyle=1,color=150
	vline,Ha_obs_henry2013,color=230,linestyle=3
	vline,lambda[sub_Ha[max_sub]],color=250,linestyle=3
	dw = 5
	sub_Ha = where(abs(Ha_obs-lambda) le dw)
	sub_NII = where(abs(NII_obs-lambda) le dw)
	sub_cont = where( ((lambda ge Ha_obs-80) and $
							(lambda le Ha_obs-10)) or $
							((lambda ge NII_obs+10) and $
							(lambda le NII_obs+80)) )

	mean_cont = mean(flux[sub_cont])
	Ha_flux = total(flux[sub_Ha]-mean_cont)
	Ha_flux_err = sqrt(total(flux_err[sub_Ha]^2.))
	NII_flux = total(flux[sub_NII]-mean_cont)
	NII_flux_err = sqrt(total(flux_err[sub_NII]^2.))
	vline,Ha_obs-dw,linestyle=2
	vline,Ha_obs+dw,linestyle=2
	vline,NII_obs-dw,linestyle=1
	vline,NII_obs+dw,linestyle=1
	xpos = NII_obs+7
	xpos1 = Ha_obs-70
	ypos = max(flux[sub_Ha])*0.8
	xyouts,xpos1,ypos,textoidl('H\alpha=')+$
			string(Ha_flux,f='(f6.2)')+textoidl('\pm')+$
			string(Ha_flux_err,f='(f6.2)'),chars=1
	xyouts,xpos1,ypos*0.8,'[NII]='+$
			string(NII_flux,f='(f6.2)')+textoidl('\pm')+$
			string(NII_flux_err,f='(f6.2)'),chars=1
	xyouts,xpos,ypos,'log([NII]/'+textoidl('H\alpha)=')+$
			string(alog10(NII_flux/Ha_flux),f='(f5.1)'),chars=1
	xyouts,xpos,ypos*0.8,'[NII]/[NII]_err='+$
			string(NII_flux/NII_flux_err,f='(f5.1)'),chars=1
	xyouts,xpos,ypos*0.6,'log(3*[NII]_err/'+textoidl('H\alpha)=')+$
			string(alog10(3.*NII_flux_err/Ha_flux),f='(f5.1)'),chars=1
	skip:
endfor

device,/close
;close,1
		
stop
end
