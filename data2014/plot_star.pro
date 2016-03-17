pro plot_star

flist = ['mosR1v4_Y_star_1012740_2014dec29_1dspec.dat',$
		'mosR1v4_Y_star_1013948_2014dec29_1dspec.dat',$
		'mosR2v4_Y_star_1014109_2014dec29_1dspec.dat',$
		'mosR2v4_Y_star_1014109_2014dec30_1dspec.dat',$
		'mosR2v4_Y_star_1036896_2014dec30_1dspec.dat',$
		'mosR3v4_Y_star_771147_2014dec30_1dspec.dat',$
		'mosR3v4_Y_star_776669_2014dec30_1dspec.dat']

@plot_setting
device,file='plot_star.ps',/color
for i=0,N_elements(flist)-1 do begin
	readcol,flist[i],w,f,ferr
	ymax = 1
	if i eq 6 then ymax = 3.5
	plot,w,f,psym=10,xr=[9750,11000],yr=[-0.4,ymax],tit=flist[i]
	oplot,w,ferr,psym=10,color=255
endfor
device,/close

stop
end
