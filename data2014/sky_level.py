import numpy as np
import pyfits
import matplotlib.pylab as plt

dir29 = '/Users/jhyoon/Tools/MosfireDRP-1.0/Data_20141229/'
dir30 = '/Users/jhyoon/Tools/MosfireDRP-1.0/Data_20141230/'

img1,hdr1 = pyfits.getdata(dir29+'m141229_0273.fits',header=True)
img2,hdr2 = pyfits.getdata(dir30+'m141230_0260.fits',header=True)

img1_sky = list(np.reshape(img1[667:830,1240:1256],(1,(830-667)*(1256-1240)))) + \
			list(np.reshape(img1[843:920,1043:1070],(1,(920-843)*(1070-1043)))) + \
			list(np.reshape(img1[1022:1078,1004:1030],(1,(1078-1022)*(1030-1004)))) + \
			list(np.reshape(img1[1287:1361,911:931],(1,(1361-1287)*(931-911))))
img2_sky = list(np.reshape(img2[667:830,1240:1256],(1,(830-667)*(1256-1240)))) + \
			list(np.reshape(img2[843:920,1043:1070],(1,(920-843)*(1070-1043)))) + \
			list(np.reshape(img2[1022:1078,1004:1030],(1,(1078-1022)*(1030-1004)))) + \
			list(np.reshape(img2[1287:1361,911:931],(1,(1361-1287)*(931-911))))
img1_sky = np.array(img1_sky)
img2_sky = np.array(img2_sky)

bins = np.arange(0,100,1)
plt.subplot(2,1,1)
plt.hist(img1_sky,bins=bins)
plt.axvline(np.mean(img1_sky),linestyle='--',color='red')
plt.subplot(2,1,2)
plt.hist(img2_sky,bins=bins)
plt.axvline(np.mean(img2_sky),linestyle='--',color='red')

rms1 = np.sqrt(np.mean(img1_sky**2.))
rms2 = np.sqrt(np.mean(img2_sky**2.))

print rms1,rms2

