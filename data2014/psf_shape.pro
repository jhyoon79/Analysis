pro psf_shape

dir = '/Users/jhyoon/Tools/MosfireDRP-1.0/Data_20141229/'

x = mrdfits(dir+'HIP43018.fits',0,hdr)
median_x = median(x)
x_sub = x[700:850,937] - median_x
x_max = max(x_sub)
x2_sub = x[700:850,951] - median_x
x2_max = max(x2_sub)
x3_sub = x[700:850,968] - median_x
x3_max = max(x3_sub)

yy = indgen(N_elements(x_sub))
@plot_setting
device,file='psf_shape.ps',/color
plot,x_sub/x_max,tit='HIP43018'
prof_fit = gaussfit(yy,x_sub/x_max,coeff,nterms=3)
prof_fit_moffat = mpfitpeak(yy,x_sub/x_max,A_moffat,/moffat)
oplot,yy,prof_fit,color=150,linestyle=1
oplot,yy,prof_fit_moffat,color=100,linestyle=2
oplot,x2_sub/x2_max,color=255,linestyle=1
oplot,x3_sub/x3_max,color=70,linestyle=2

x = mrdfits(dir+'HD13027.fits',0,hdr)
median_x = median(x)
x_sub = x[80:180,343] - median_x
x_max = max(x_sub)
x2_sub = x[80:180,350] - median_x
x2_max = max(x2_sub)
x3_sub = x[80:180,360] - median_x
x3_max = max(x3_sub)
plot,x_sub/x_max,tit='HD13027'
oplot,x2_sub/x2_max,color=255,linestyle=1
oplot,x3_sub/x3_max,color=70,linestyle=2
device,/close


stop
end

