pro image_inspect,star=star

dir = '/Users/jhyoon/Tools/MosfireDRP-1.0/REDUX/'
readcol,'fname_mosR123_strict.list',flist,z_spec,z_phot,ypix,gr,f='(a,f,f,i,a)'
N_file = N_elements(flist)
dir_ACS = '/Users/jhyoon/Dropbox/Research/mosfire/t_exp4obs/ACS_cutout/'
dir_imacs = '/Users/jhyoon/Dropbox/Research/mosfire/t_exp4obs/imacs_img/'

dir_img_maskF = '/Users/jhyoon/Dropbox/Research/mosfire/flux_calib/stamps/targets/'
dir_eps_maskF = '/Users/jhyoon/Tools/MosfireDRP-1.0/REDUX/mosFv2/2013jan4/Y/'
flist_maskF = ['F15','F17','F31','F32','maskF_978834']

if not keyword_set(star) then begin
	;for i=0,N_file-1 do begin
	for i=30,N_file-1 do begin
		sep = strpos(flist[i],'_eps')
		sep2 = strpos(flist[i],'_Y_')
		spawn,'ds9 -tile grid layout 1 3 -scale mode zscale -zoom 0.8 '+ $
				dir_imacs+strmid(flist[i],sep2+3,20)+'.fits -invert '+$
				' -regions '+dir_ACS+'mosR1v4_SlitRegions.reg '+$
				' -regions '+dir_ACS+'mosR2v4_SlitRegions.reg '+$
				' -regions '+dir_ACS+'mosR3v4_SlitRegions.reg '+$
				dir_ACS+$
				strmid(flist[i],30,40)+'_acs_I_mosaic_30mas_sci.fits '+$
				dir+flist[i]+'_eps.fits '
	endfor
	for i=0,N_elements(flist_maskF)-1 do begin
		sep = strpos(flist_maskF[i],'_eps')
		spawn,'ds9 -tile grid layout 1 3 -scale mode zscale -zoom 0.8 '+ $
				dir_imacs+'F_'+flist_maskF[i]+'.fits -invert '+$
				' -regions '+dir_img_maskF+'mosFv2_SlitRegions.reg '+$
				dir_img_maskF+$
				flist_maskF[i]+'_acs_I_mosaic_30mas_sci.fits '+$
				dir_eps_maskF+'mosFv2_Y_'+flist_maskF[i]+'_eps.fits '
	endfor	
endif



if keyword_set(star) then begin
	spawn,'ds9 -tile grid layout 2 3 -scale mode zscale '+ $
			dir_ACS+$
			strmid(flist[28],30,40)+'_acs_I_mosaic_30mas_sci.fits -invert '+$
			' -regions '+dir_ACS+'mosR1v4_SlitRegions.reg '+$
			' -regions '+dir_ACS+'mosR2v4_SlitRegions.reg '+$
			' -regions '+dir_ACS+'mosR3v4_SlitRegions.reg '+$
			dir_ACS+$
			strmid(flist[29],30,40)+'_acs_I_mosaic_30mas_sci.fits -invert '+$
			' -regions '+dir_ACS+'mosR1v4_SlitRegions.reg '+$
			' -regions '+dir_ACS+'mosR2v4_SlitRegions.reg '+$
			' -regions '+dir_ACS+'mosR3v4_SlitRegions.reg '+$
			dir_ACS+$
			strmid(flist[57],30,40)+'_acs_I_mosaic_30mas_sci.fits -invert '+$
			' -regions '+dir_ACS+'mosR1v4_SlitRegions.reg '+$
			' -regions '+dir_ACS+'mosR2v4_SlitRegions.reg '+$
			' -regions '+dir_ACS+'mosR3v4_SlitRegions.reg '+$
			dir_ACS+$
			strmid(flist[58],30,40)+'_acs_I_mosaic_30mas_sci.fits -invert '+$
			' -regions '+dir_ACS+'mosR1v4_SlitRegions.reg '+$
			' -regions '+dir_ACS+'mosR2v4_SlitRegions.reg '+$
			' -regions '+dir_ACS+'mosR3v4_SlitRegions.reg '+$
			dir_ACS+$
			strmid(flist[72],30,40)+'_acs_I_mosaic_30mas_sci.fits -invert '+$
			' -regions '+dir_ACS+'mosR1v4_SlitRegions.reg '+$
			' -regions '+dir_ACS+'mosR2v4_SlitRegions.reg '+$
			' -regions '+dir_ACS+'mosR3v4_SlitRegions.reg '+$
			dir_ACS+$
			strmid(flist[73],30,40)+'_acs_I_mosaic_30mas_sci.fits -invert '+$
			' -regions '+dir_ACS+'mosR1v4_SlitRegions.reg '+$
			' -regions '+dir_ACS+'mosR2v4_SlitRegions.reg '+$
			' -regions '+dir_ACS+'mosR3v4_SlitRegions.reg &'
endif



stop
end
