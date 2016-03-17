pro sky_level

dir29 = '/Users/jhyoon/Tools/MosfireDRP-1.0/Data_20141229/'
dir30 = '/Users/jhyoon/Tools/MosfireDRP-1.0/Data_20141230/'

@plot_setting
device,file='sky_level.ps',/color;,xsize=15,ysize=20
!p.multi=[0,2,3]

for i=0,2 do begin
	if i eq 0 then begin 
		fits_read,dir29+'m141229_0273.fits',img1,hdr1
		fits_read,dir30+'m141230_0260.fits',img2,hdr2
	endif else if i eq 1 then begin 
		fits_read,dir29+'m141229_0293.fits',img1,hdr1
		fits_read,dir30+'m141230_0280.fits',img2,hdr2
	endif else if i eq 2 then begin 
		fits_read,dir29+'m141229_0313.fits',img1,hdr1
		fits_read,dir30+'m141230_0300.fits',img2,hdr2
	endif
	time1 = sxpar(hdr1,'UTC')
	time2 = sxpar(hdr2,'UTC')

	img1_sub1 = img1[1240:1256,667:830]
	img1_sub2 =	img1[1043:1070,843:920]
	img1_sub3 = img1[1004:1030,1022:1078]
	img1_sub4 = img1[911:931,1287:1361]
	img1_sub5 = img1[989:1019,1460:1625]
	img1_sky = [reform(reform(img1_sub1,1,N_elements(img1_sub1))), $
				reform(reform(img1_sub2,1,N_elements(img1_sub2))), $
				reform(reform(img1_sub3,1,N_elements(img1_sub3))), $
				reform(reform(img1_sub4,1,N_elements(img1_sub4))), $
				reform(reform(img1_sub5,1,N_elements(img1_sub5)))]
	img2_sub1 = img2[1240:1256,667:830]
	img2_sub2 =	img2[1043:1070,843:920]
	img2_sub3 = img2[1004:1030,1022:1078]
	img2_sub4 = img2[911:931,1287:1361]
	img2_sub5 = img2[989:1019,1460:1625]
	img2_sky = [reform(reform(img2_sub1,1,N_elements(img2_sub1))), $
				reform(reform(img2_sub2,1,N_elements(img2_sub2))), $
				reform(reform(img2_sub3,1,N_elements(img2_sub3))), $
				reform(reform(img2_sub4,1,N_elements(img2_sub4))), $
				reform(reform(img2_sub5,1,N_elements(img2_sub5)))]

	t_moonset1 = '11:10:00'	
	t_moonset2 = '12:07:00'	
	xtit = ''
	if i ge 2 then xtit='count'
	multiplot
	plothist,img1_sky,bin=1,xr=[-10,99],yr=[0,1300],xtit=xtit
	median_sky1 = median(img1_sky)
	mean_sky1 = mean(img1_sky)
	vline,median_sky1,linestyle=2
	vline,mean_sky1,linestyle=1
	legend,['night1','t_obs(UTC)='+time1,$
			't_moonset(UTC)='+t_moonset1,'median='+strtrim(median_sky1,2),$
			'mean='+strtrim(mean_sky1,2)],/right,box=0,chars=0.8
	multiplot
	plothist,img2_sky,bin=1,xr=[-10,99],yr=[0,1300],xtit=xtit
	median_sky2 = median(img2_sky)
	mean_sky2 = mean(img2_sky)
	vline,median_sky2,linestyle=2
	vline,mean_sky2,linestyle=1
	legend,['night2','t_obs(UTC)='+time2,$
			't_moonset(UTC)='+t_moonset2,'median='+strtrim(median_sky2,2),$
			'mean='+strtrim(mean_sky2,2)],/right,box=0,chars=0.8
endfor
multiplot,/reset
erase
device,/close

stop
end
