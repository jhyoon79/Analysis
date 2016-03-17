pro test

file1 = '/Users/jhyoon/Tools/MosfireDRP-1.0/REDUX/mosR2v4/2014dec29/Y/mosR2v4_Y_star_1014109_eps.fits'
file2 = '/Users/jhyoon/Tools/MosfireDRP-1.0/REDUX/mosR2v4/2014dec30/Y/mosR2v4_Y_star_1014109_eps.fits'

@plot_setting
device,file='test.ps',/color

fits_read,file1,x1,hdr1
fits_read,file2,x2,hdr2


x1_sub = total(x1[300:600,40:60],2)
x2_sub = total(x2[300:600,40:60],2)
ymin = min([x1_sub,x2_sub])
ymax = max([x1_sub,x2_sub])
plot,x1_sub,psym=10,yr=[ymin,ymax]
hline,mean(x1_sub)
oplot,x2_sub,psym=10,color=255
hline,mean(x2_sub),color=255

x1_sub = total(x1[600:800,40:60],2)
x2_sub = total(x2[600:800,40:60],2)
ymin = min([x1_sub,x2_sub])
ymax = max([x1_sub,x2_sub])
plot,x1_sub,psym=10,yr=[ymin,ymax]
hline,mean(x1_sub)
oplot,x2_sub,psym=10,color=255
hline,mean(x2_sub),color=255

x1_sub = total(x1[800:980,40:60],2)
x2_sub = total(x2[800:980,40:60],2)
ymin = min([x1_sub,x2_sub])
ymax = max([x1_sub,x2_sub])
plot,x1_sub,psym=10,yr=[ymin,ymax]
hline,mean(x1_sub)
oplot,x2_sub,psym=10,color=255
hline,mean(x2_sub),color=255

xr1 = [455,626,725,780]
xr2 = [475,641,754,850]
for i=0,N_elements(xr1)-1 do begin
	x1_sub = total(x1[xr1[i]:xr2[i],40:60],2)
	x2_sub = total(x2[xr1[i]:xr2[i],40:60],2)
	ymin = min([x1_sub,x2_sub])
	ymax = max([x1_sub,x2_sub])
	plot,x1_sub,psym=10,yr=[ymin,ymax]
	hline,mean(x1_sub)
	oplot,x2_sub,psym=10,color=255
	hline,mean(x2_sub),color=255
endfor
erase

file1 = '/Users/jhyoon/Tools/MosfireDRP-1.0/REDUX/mosR2v4/2014dec29/Y/mosR2v4_Y_star_1014109_eps_1dspec.fits'
file2 = '/Users/jhyoon/Tools/MosfireDRP-1.0/REDUX/mosR2v4/2014dec30/Y/mosR2v4_Y_star_1014109_eps_1dspec.fits'
f1 = 'mosR2v4_Y_star_1014109_2014dec29_1dspec.dat'
f2 = 'mosR2v4_Y_star_1014109_2014dec30_1dspec.dat'
f3 = 'mosR2v4_Y_star_1014109_2014dec29_1dspec2.dat'
f4 = 'mosR2v4_Y_star_1014109_2014dec30_1dspec2.dat'
fits_read,file1,x1,hdr1
fits_read,file2,x2,hdr2
readcol,f1,w_optimal1,f_optimal1
readcol,f2,w_optimal2,f_optimal2
readcol,f3,w_optimal3,f_optimal3
readcol,f4,w_optimal4,f_optimal4

for i=0,1100,200 do begin
	ymin = min([x1[i:i+200],x2[i:i+200],f_optimal1[i:i+200],f_optimal2[i:i+200]])
	ymax = max([x1[i:i+200],x2[i:i+200],f_optimal1[i:i+200],f_optimal2[i:i+200]])
	multiplot,[1,4]
	plot,x1,psym=10,tit='IRAF 1D x='+strtrim(i,2),$
		xr=[i,i+200],yr=[ymin,ymax]
	oplot,x2,psym=10,color=255
	multiplot
	plot,x1,psym=10,$
		xr=[i,i+200],yr=[ymin,ymax]
	oplot,f_optimal1,psym=10,color=100
	oplot,f_optimal3,psym=10,color=150
	multiplot
	plot,x2,psym=10,color=255,$
		xr=[i,i+200],yr=[ymin,ymax]
	oplot,f_optimal2,psym=10,color=70
	oplot,f_optimal4,psym=10,color=200
	multiplot
	plot,f_optimal1,psym=10,color=100,$
		xr=[i,i+200],yr=[ymin,ymax]
	oplot,f_optimal2,psym=10,color=70
	legend,['mean all n1='+strtrim(mean(x1),2),$
			'mean all n2='+strtrim(mean(x2),2),$
			'mean n1='+strtrim(mean(x1[i:i+200]),2),$
			'mean n2='+strtrim(mean(x2[i:i+200]),2),$
			'mean opt n1='+strtrim(mean(f_optimal1[i:i+200]),2),$
			'mean opt n2='+strtrim(mean(f_optimal2[i:i+200]),2),$
			'mean opt2 n1='+strtrim(mean(f_optimal3[i:i+200]),2),$
			'mean opt2 n2='+strtrim(mean(f_optimal4[i:i+200]),2)],$
			color=[0,255,0,255,100,70],chars=0.6
	multiplot,/reset
	erase
endfor

flist = ['mosR2v4_Y_star_1014109_2014dec29_1dspec.dat',$
			'mosR2v4_Y_star_1014109_2014dec30_1dspec.dat']
readcol,flist[0],w1,f1,e1;,f2,e2
readcol,flist[1],w2,f2,e2;,f2,e2

plot,w1,f1,psym=10,xr=[10000,10300]
oplot,w2,f2,psym=10,color=255
plot,w1,f1,psym=10,xr=[10300,10600]
oplot,w2,f2,psym=10,color=255
plot,w1,f1,psym=10,xr=[10600,10900]
oplot,w2,f2,psym=10,color=255

print,mean(f1),mean(f2)
print,mean(f1)/mean(f2)

readcol,'../flux_calib/sensfunc_2014dec29_mean.dat',w29,f29
readcol,'../flux_calib/sensfunc_2014dec30_mean.dat',w30,f30

plot,w29,f29,psym=10,xr=[10000,10500]
oplot,w30,f30,psym=10,color=255
plot,w29,f29,psym=10,xr=[10500,11000]
oplot,w30,f30,psym=10,color=255
print,mean(f29)/mean(f30)


device,/close

stop
end
