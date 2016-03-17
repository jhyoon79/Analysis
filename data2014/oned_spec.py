import pyfits
import numpy as np
import matplotlib.pylab as plt
from subprocess import call
import pdb
from astropy.io import ascii
from matplotlib.backends.backend_pdf import PdfPages
fig_pdf = PdfPages('oned_spec.pdf')

flist = ['/Users/jhyoon/Tools/MosfireDRP-1.0/REDUX/mosR1v4/2014dec29/Y/mosR1v4_Y_R1_1009808',\
			'/Users/jhyoon/Tools/MosfireDRP-1.0/REDUX/mosR1v4/2014dec29/Y/mosR1v4_Y_R1_1015084',\
			'/Users/jhyoon/Tools/MosfireDRP-1.0/REDUX/mosR1v4/2014dec29/Y/mosR1v4_Y_R2_1015516',\
			'/Users/jhyoon/Tools/MosfireDRP-1.0/REDUX/mosR2v4/2014dec29/Y/mosR2v4_Y_R2_1015516',\
			'/Users/jhyoon/Tools/MosfireDRP-1.0/REDUX/mosR2v4/2014dec30/Y/mosR2v4_Y_R2_1015516',\
			'/Users/jhyoon/Tools/MosfireDRP-1.0/REDUX/mosR3v4/2014dec30/Y/mosR3v4_Y_R2_1015516',\
			'/Users/jhyoon/Tools/MosfireDRP-1.0/REDUX/mosR2v4/2014dec29/Y/mosR2v4_Y_R2_1010583',\
			'/Users/jhyoon/Tools/MosfireDRP-1.0/REDUX/mosR2v4/2014dec30/Y/mosR2v4_Y_R2_1010583',\
			'/Users/jhyoon/Tools/MosfireDRP-1.0/REDUX/mosR2v4/2014dec29/Y/mosR2v4_Y_R2_1010994',\
			'/Users/jhyoon/Tools/MosfireDRP-1.0/REDUX/mosR2v4/2014dec30/Y/mosR2v4_Y_R2_1010994',\
			'/Users/jhyoon/Tools/MosfireDRP-1.0/REDUX/mosR2v4/2014dec29/Y/mosR2v4_Y_R2_1013478',\
			'/Users/jhyoon/Tools/MosfireDRP-1.0/REDUX/mosR2v4/2014dec30/Y/mosR2v4_Y_R2_1013478',\
			'/Users/jhyoon/Tools/MosfireDRP-1.0/REDUX/mosR3v4/2014dec30/Y/mosR3v4_Y_R3_772773',\
			'/Users/jhyoon/Tools/MosfireDRP-1.0/REDUX/mosR3v4/2014dec30/Y/mosR3v4_Y_R3_775548']
y1_list = [22,28,39,36,32,90 ,52,50,85 ,84 ,53,50,46,41]
y2_list = [45,49,62,58,56,114,75,73,110,107,76,75,69,64]
plot_list = [1,1,0,0,0,1,0,1,0,1,0,1,1,1]

henry2013 = np.loadtxt('/Users/jhyoon/Dropbox/Research/mosfire/t_exp4obs/t_exp_Henry2013.dat',dtype=str,skiprows=1)
id_henry2013 = henry2013[:,0]
z_henry2013 = np.array(map(float,henry2013[:,3]))
pre_target_name = '9999999'

for j in range(len(flist)):
	fname = flist[j]
	target_name = fname[74:]
	y1 = y1_list[j]
	y2 = y2_list[j]
	print j,target_name
	#call_name = '%s_eps.fits'%(fname)
	#call(['ds9',call_name])

	eps,epshdr = pyfits.getdata(fname+'_eps.fits',0,header=True)
	sig,varhdr = pyfits.getdata(fname+'_sig.fits',0,header=True)
	snr,snrhdr = pyfits.getdata(fname+'_snrs.fits',0,header=True)
	var = sig**2.
	wave0 = epshdr['CRVAL1']
	dw = epshdr['CD1_1']
	Nw = eps.shape[1]
	wavelength = wave0+np.arange(Nw)*dw

	sub_match = list(id_henry2013).index(target_name)
	z = z_henry2013[sub_match]
	lambda_Ha = 6563.*(1.+z)
	lambda_NII1 = 6548.*(1.+z)
	lambda_NII2 = 6584.*(1.+z)

	Nw = len(eps[0,:])
	spec1d = np.zeros(Nw)
	spec1d_weight = np.zeros(Nw)
	var1d = np.zeros(Nw)
	snr1d = np.zeros(Nw)
	for i in range(Nw):
		spec1d[i] = sum(eps[y1:y2,i])
		spec1d_weight[i] = (y2-y1+1.)*sum(eps[y1:y2,i]/var[y1:y2,i])/sum(1./var[y1:y2,i])
		var1d[i] = sum(var[y1:y2,i])
		snr1d[i] = np.sqrt(sum(snr[y1:y2,i]**2.)/sum(snr[y1:y2,i]))

	if (pre_target_name != target_name):
		spec1d_coadd = spec1d
		var1d_coadd = var1d
	if (pre_target_name == target_name):
		spec1d_coadd += spec1d
		var1d_coadd += var1d

	if plot_list[j]==1:
		
		ascii.write([wavelength,spec1d_coadd,np.sqrt(var1d_coadd)],\
			target_name+'.dat',names=['# lambda','flux(e-/s)','stddev'],\
			formats={'# lambda':'%10.3f','flux(e-/s)':'%8.5f',\
			'stddev':'%8.5f'},Writer=ascii.FixedWidth,delimiter=None)

		plt.clf()
		plt.plot(wavelength,spec1d_coadd,color='black',drawstyle='steps-mid')
		#plt.plot(wavelength,spec1d_weight,color='green',linestyle='-',drawstyle='steps-mid')
		plt.suptitle(fname[61:])
		plt.plot(wavelength,np.sqrt(var1d_coadd),color='red',drawstyle='steps-mid')
		#plt.plot(snr1d,color='blue',drawstyle='steps-mid')
		xmin = 9750
		xmax = 11250
		plt.xlim(xmin,xmax)
		plt.axvline(lambda_Ha,color='black',linestyle=':')
		plt.axvline(lambda_NII1,color='black',linestyle=':')
		plt.axvline(lambda_NII2,color='black',linestyle=':')
		fig_pdf.savefig()
		plt.clf()
		plt.plot(wavelength,spec1d_coadd,color='black',drawstyle='steps-mid')
		#plt.plot(wavelength,spec1d_weight,color='green',linestyle='-',drawstyle='steps-mid')
		plt.plot(wavelength,np.sqrt(var1d_coadd),color='red',drawstyle='steps-mid')
		plt.suptitle(fname[61:])
		plt.axvline(lambda_Ha,color='black',linestyle=':')
		plt.axvline(lambda_NII1,color='black',linestyle=':')
		plt.axvline(lambda_NII2,color='black',linestyle=':')
		xmin = lambda_NII1-30
		xmax = lambda_NII2+30
		plt.xlim(xmin,xmax)
		fig_pdf.savefig()
	
	pre_target_name = target_name

#	continuum = np.array(list(spec1d[1067:1078]) + \
#								list(spec1d[952:990]) + \
#								list(spec1d[1088:1101]) + \
#								list(spec1d[1241:1261]) + \
#								list(spec1d[1313:1343]))
#	rms = np.sqrt(np.mean(continuum**2.))
#	Ha1 = 1161
#	Ha2 = 1174
#	Ha_flux = sum(spec1d[Ha1:Ha2])
#	Ha_sigma = np.sqrt(sum(var1d[Ha1:Ha2]**2.))
#	plt.text(xmin*1.1,3.0,'Ha=%6.3f'%(Ha_flux))
#	plt.text(xmin*1.1,2.5,'sigma_Ha=%6.3f'%(rms))
#	plt.text(xmin*1.1,2.0,'S/N(Ha)=%6.3f'%(Ha_flux/rms))
#	plt.text(xmin*1.1,1.5,'Ha*0.05=%6.3f'%(Ha_flux*0.05))
#	NII1 = 1193
#	NII2 = 1206
#	plt.axvline(NII1,color='black',linestyle=':')
#	plt.axvline(NII2,color='black',linestyle=':')
#	NII_flux = sum(spec1d[NII1:NII2])
#	NII_sigma = np.sqrt(sum(var1d[NII1:NII2]**2.))
#	plt.text(xmax*0.9,3.0,'NII=%6.3f'%(NII_flux))
#	plt.text(xmax*0.9,2.5,'sigma_NII=%6.3f'%(rms))
#	plt.text(xmax*0.9,2.0,'3*sigma_NII=%6.3f'%(3*rms))
#	plt.text(xmax*0.9,1.5,'S/N(NII)=%6.3f'%(NII_flux/rms))
#	NII1 = 1136
#	NII2 = 1149
#	plt.axvline(NII1,color='black',linestyle=':')
#	plt.axvline(NII2,color='black',linestyle=':')
#	plt.show()
#	pdb.set_trace()

fig_pdf.close()

