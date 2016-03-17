pro set_continuum,i,lambda_optimal,lambda_iraf,$
						sub_continuum_optimal,sub_continuum_iraf
if i eq 0 then begin
	sub_continuum_optimal = where( $
				((lambda_optimal gt 10475) and (lambda_optimal lt 10505)) or $
				((lambda_optimal gt 10530) and (lambda_optimal lt 10570)) or $
				((lambda_optimal gt 10615) and (lambda_optimal lt 10635)) or $
				((lambda_optimal gt 10660) and (lambda_optimal lt 10680)))
	sub_continuum_iraf = where( $
				((lambda_iraf gt 10475) and (lambda_iraf lt 10505)) or $
				((lambda_iraf gt 10530) and (lambda_iraf lt 10570)) or $
				((lambda_iraf gt 10615) and (lambda_iraf lt 10635)) or $
				((lambda_iraf gt 10660) and (lambda_iraf lt 10680)))
endif else if i eq 1 then begin
	sub_continuum_optimal = where( $
				((lambda_optimal gt 10865) and (lambda_optimal lt 10895)) or $
				((lambda_optimal gt 10905) and (lambda_optimal lt 10922)) or $
				((lambda_optimal gt 10955) and (lambda_optimal lt 10970)) or $
				((lambda_optimal gt 11097) and (lambda_optimal lt 11140)) or $
				((lambda_optimal gt 11160) and (lambda_optimal lt 11195)))
	sub_continuum_iraf = where( $
				((lambda_iraf gt 10865) and (lambda_iraf lt 10895)) or $
				((lambda_iraf gt 10905) and (lambda_iraf lt 10922)) or $
				((lambda_iraf gt 10955) and (lambda_iraf lt 10970)) or $
				((lambda_iraf gt 11097) and (lambda_iraf lt 11140)) or $
				((lambda_iraf gt 11160) and (lambda_iraf lt 11195)))
endif else if (i eq 2 or i eq 4 or i eq 5) then begin
	sub_continuum_optimal = where( $
				((lambda_optimal gt 10905) and (lambda_optimal lt 10921)) or $
				((lambda_optimal gt 10933) and (lambda_optimal lt 10945)) or $
				((lambda_optimal gt 10955) and (lambda_optimal lt 10970)) or $
				((lambda_optimal gt 11095) and (lambda_optimal lt 11138)) or $
				((lambda_optimal gt 11160) and (lambda_optimal lt 11195)))
	sub_continuum_iraf = where( $
				((lambda_iraf gt 10905) and (lambda_iraf lt 10921)) or $
				((lambda_iraf gt 10933) and (lambda_iraf lt 10945)) or $
				((lambda_iraf gt 10955) and (lambda_iraf lt 10970)) or $
				((lambda_iraf gt 11095) and (lambda_iraf lt 11138)) or $
				((lambda_iraf gt 11160) and (lambda_iraf lt 11195)))
endif else if i eq 3 then begin
	sub_continuum_optimal = where( $
				((lambda_optimal gt 10595) and (lambda_optimal lt 10637)) or $
				((lambda_optimal gt 10660) and (lambda_optimal lt 10680)) or $
				((lambda_optimal gt 10808) and (lambda_optimal lt 10830)) or $
				((lambda_optimal gt 10905) and (lambda_optimal lt 10922)))
	sub_continuum_iraf = where( $
				((lambda_iraf gt 10595) and (lambda_iraf lt 10637)) or $
				((lambda_iraf gt 10660) and (lambda_iraf lt 10680)) or $
				((lambda_iraf gt 10808) and (lambda_iraf lt 10830)) or $
				((lambda_iraf gt 10905) and (lambda_iraf lt 10922)))
endif else if i eq 6 then begin
	sub_continuum_optimal = where( $
				((lambda_optimal gt 10780) and (lambda_optimal lt 10825)) or $
				((lambda_optimal gt 10863) and (lambda_optimal lt 10870)) or $
				((lambda_optimal gt 10905) and (lambda_optimal lt 10920)) or $
				((lambda_optimal gt 10980) and (lambda_optimal lt 10992)) or $
				((lambda_optimal gt 11035) and (lambda_optimal lt 11050)) or $
				((lambda_optimal gt 11095) and (lambda_optimal lt 11030)))
	sub_continuum_iraf = where( $
				((lambda_iraf gt 10780) and (lambda_iraf lt 10825)) or $
				((lambda_iraf gt 10863) and (lambda_iraf lt 10870)) or $
				((lambda_iraf gt 10905) and (lambda_iraf lt 10920)) or $
				((lambda_iraf gt 10980) and (lambda_iraf lt 10992)) or $
				((lambda_iraf gt 11035) and (lambda_iraf lt 11050)) or $
				((lambda_iraf gt 11095) and (lambda_iraf lt 11030)))
endif

return
end

pro set_emission,i,lambda_optimal,lambda_iraf,$
					sub_Ha_optimal,sub_Ha_iraf,$
					sub_NII_optimal,sub_NII_iraf

if i eq 0 then begin
	Ha1 = 10601
	Ha2 = 10611
	NII1 = 10638
	NII2 = 10642
	sub_Ha_optimal = where( $
				(lambda_optimal gt Ha1) and (lambda_optimal lt Ha2))
	sub_Ha_iraf = where( $
				(lambda_iraf gt Ha1) and (lambda_iraf lt Ha2))
	sub_NII_optimal = where( $
				(lambda_optimal gt NII1) and (lambda_optimal lt NII2))
	sub_NII_iraf = where( $
				(lambda_iraf gt NII1) and (lambda_iraf lt NII2))
endif else if i eq 1 then begin
	Ha1 = 11008
	Ha2 = 11025
	NII1 = Ha1 + 34
	NII2 = Ha2 + 34
	sub_Ha_optimal = where( $
				(lambda_optimal gt Ha1) and (lambda_optimal lt Ha2))
	sub_Ha_iraf = where( $
				(lambda_iraf gt Ha1) and (lambda_iraf lt Ha2))
	sub_NII_optimal = where( $
				(lambda_optimal gt NII1) and (lambda_optimal lt NII2))
	sub_NII_iraf = where( $
				(lambda_iraf gt NII1) and (lambda_iraf lt NII2))
endif else if (i eq 2 or i eq 4 or i eq 5) then begin
	Ha1 = 11016
	Ha2 = 11027
	NII1 = 11056
	NII2 = 11060
	sub_Ha_optimal = where( $
				(lambda_optimal gt Ha1) and (lambda_optimal lt Ha2))
	sub_Ha_iraf = where( $
				(lambda_iraf gt Ha1) and (lambda_iraf lt Ha2))
	sub_NII_optimal = where( $
				(lambda_optimal gt NII1) and (lambda_optimal lt NII2))
	sub_NII_iraf = where( $
				(lambda_iraf gt NII1) and (lambda_iraf lt NII2))
endif else if i eq 3 then begin
	Ha1 = 10747
	Ha2 = 10763
	NII1 = Ha1 + 34
	NII2 = Ha2 + 34
	sub_Ha_optimal = where( $
				(lambda_optimal gt Ha1) and (lambda_optimal lt Ha2))
	sub_Ha_iraf = where( $
				(lambda_iraf gt Ha1) and (lambda_iraf lt Ha2))
	sub_NII_optimal = where( $
				(lambda_optimal gt NII1) and (lambda_optimal lt NII2))
	sub_NII_iraf = where( $
				(lambda_iraf gt NII1) and (lambda_iraf lt NII2))
endif else if i eq 6 then begin
	Ha1 = 10956
	Ha2 = 10968
	NII1 = 10993
	NII2 = 10998
	sub_Ha_optimal = where( $
				(lambda_optimal gt Ha1) and (lambda_optimal lt Ha2))
	sub_Ha_iraf = where( $
				(lambda_iraf gt Ha1) and (lambda_iraf lt Ha2))
	sub_NII_optimal = where( $
				(lambda_optimal gt NII1) and (lambda_optimal lt NII2))
	sub_NII_iraf = where( $
				(lambda_iraf gt NII1) and (lambda_iraf lt NII2))
endif

return
end

pro oned_spec_compare

;dir = '/Users/jhyoon/Tools/MosfireDRP-1.0/REDUX/mosFv2/2013jan4/Y/'
chr_rdtbl,'fname_1dspec_iraf.list',0,flist
flist_iraf = flist[0:6]
flist_optimal = flist[7:*]
N_file = N_elements(flist_iraf)

xmin = [10560,10960,10960,10700,10960,10960,10900]
xmax = [10660,11070,11080,10820,11080,11080,11020]
@plot_setting
!p.charsize=1.2
device,file='oned_spec_compare.ps',/color,xs=15,ys=19,$
		xoffset=0.5,yoffset=0.5
for i=0,N_file-1 do begin
	print,flist_iraf[i]

	f_iraf = flist_iraf[i]+'.fits'
	f_iraf_weight= flist_iraf[i]+'_weight.fits'
	f_optimal = flist_optimal[i]

	flux_iraf = mrdfits(f_iraf,0,hdr_iraf,/silent)
	w0_iraf = sxpar(hdr_iraf,'CRVAL1')
	dw_iraf = sxpar(hdr_iraf,'CD1_1')
	lambda_iraf = w0_iraf + indgen(N_elements(flux_iraf))*dw_iraf	
	flux_iraf_weight = mrdfits(f_iraf_weight,0,hdr_iraf_weight,/silent)
	w0_iraf_weight = sxpar(hdr_iraf_weight,'CRVAL1')
	dw_iraf_weight = sxpar(hdr_iraf_weight,'CD1_1')
	lambda_iraf_weight = w0_iraf_weight + $
						indgen(N_elements(flux_iraf_weight))*dw_iraf_weight
	readcol,f_optimal,lambda_optimal,flux_optimal,f='(f,f)',/silent

	set_continuum,i,lambda_optimal,lambda_iraf,$
					sub_continuum_optimal,sub_continuum_iraf
	set_emission,i,lambda_optimal,lambda_iraf,$
					sub_Ha_optimal,sub_Ha_iraf,$
					sub_NII_optimal,sub_NII_iraf

	;;; estimate rms of continuum
	rms_optimal = sqrt(mean(flux_optimal[sub_continuum_optimal]^2.))
	rms_iraf = sqrt(mean(flux_iraf[sub_continuum_iraf]^2.))
	rms_iraf_weight = sqrt(mean(flux_iraf_weight[sub_continuum_iraf]^2.))

	;;; estimate Ha flux
	Ha_optimal = total(flux_optimal[sub_Ha_optimal])
	Ha_iraf = total(flux_iraf[sub_Ha_iraf])
	Ha_iraf_weight = total(flux_iraf_weight[sub_Ha_iraf])
	xp1_Ha = min(lambda_iraf[sub_Ha_iraf])
	xp2_Ha = max(lambda_iraf[sub_Ha_iraf])
	NII_optimal = total(flux_optimal[sub_NII_optimal])
	NII_iraf = total(flux_iraf[sub_NII_iraf])
	NII_iraf_weight = total(flux_iraf_weight[sub_NII_iraf])
	xp1_NII = min(lambda_iraf[sub_NII_iraf])
	xp2_NII = max(lambda_iraf[sub_NII_iraf])

	multiplot,[1,2],ygap=0.02,/doxaxis
	plot,lambda_optimal,flux_optimal,psym=10,tit=flist_optimal[i]
	oplot,lambda_iraf,flux_iraf,psym=10,color=255
	oplot,lambda_iraf_weight,flux_iraf_weight,psym=10,color=70
	legend,['optimal, rms='+string(rms_optimal,f='(f6.3)'),$
			'apsum, rms= '+string(rms_iraf,f='(f6.3)'),$
			'apsum(var) rms= '+string(rms_iraf_weight,f='(f6.3)')],$
			linestyle=0,color=[0,255,70],box=0

	multiplot,/doxaxis
	plot,lambda_optimal,flux_optimal,psym=10,$
		xr=[xmin[i],xmax[i]]
	oplot,lambda_iraf,flux_iraf,psym=10,color=255
	oplot,lambda_iraf_weight,flux_iraf_weight,psym=10,color=70
;	oplot,lambda_optimal[sub_continuum_optimal],$
;			flux_optimal[sub_continuum_optimal],psym=7
	vline,xp1_Ha,linestyle=1
	vline,xp2_Ha,linestyle=1
	vline,xp1_NII,linestyle=1
	vline,xp2_NII,linestyle=1
	yp = max(flux_iraf[sub_Ha_iraf])*0.85
	xyouts,xp2_Ha,yp,'Ha='+string(Ha_optimal,f='(f5.2)')+$
			' S/N='+string(Ha_optimal/rms_optimal,f='(f5.1)')
	xyouts,xp2_Ha,yp*0.85,'Ha='+string(Ha_iraf,f='(f5.2)')+$
			' S/N='+string(Ha_iraf/rms_iraf,f='(f5.1)'),color=255
	xyouts,xp2_Ha,yp*0.70,'Ha='+string(Ha_iraf_weight,f='(f5.2)')+$
			' S/N='+string(Ha_iraf_weight/rms_iraf_weight,f='(f5.1)'),color=70

	multiplot,/reset
	erase
endfor
device,/close

stop
end
