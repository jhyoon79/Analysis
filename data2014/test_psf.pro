pro test_psf

x = indgen(200)/10.+15
y = indgen(200)/10.+15
erase
y = 18
u_2d = ((x-18.)/3.1)^2.+((y-18.)/3.1)^2.
f_moffat_2d = 0 + 9755/(u_2d+1.)^2.6
multiplot,[2,1]
plot,x,f_moffat_2d/max(f_moffat_2d),title='Moffat'
y = 17
u_2d = ((x-18.)/3.1)^2.+((y-18.)/3.1)^2.
f_moffat_2d = 9755/(u_2d+1.)^2.6
oplot,x,f_moffat_2d/max(f_moffat_2d),color=255,linestyle=2

x = indgen(100)/10.
y = indgen(100)/10.
y = 5
u = ((x-5.)/2.)^2.+((y-5.)/2.)^2.
f_gauss = exp(-1.*u)
multiplot
plot,x,f_gauss/max(f_gauss),xr=[0,10],tit='Gaussian'
y = 4
u = ((x-5.)/2.)^2.+((y-5.)/2.)^2.
f_gauss = exp(-1.*u)
oplot,x,f_gauss/max(f_gauss),color=255,linestyle=2

stop
end
