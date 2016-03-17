import pyfits
import numpy as np
import matplotlib.pylab as plt
from subprocess import call
import pdb
from astropy.io import ascii
from matplotlib.backends.backend_pdf import PdfPages
fig_pdf = PdfPages('oned_spec_plot.pdf')

flist = ['mosR1v4_Y_R1_1009808',\
			'mosR1v4_Y_R1_1015084',\
			'mosR1v4_Y_R2_1015516',\
			'mosR2v4_Y_R2_1015516',\
			'mosR2v4_Y_R2_1015516',\
			'mosR3v4_Y_R2_1015516',\
			'mosR2v4_Y_R2_1010583',\
			'mosR2v4_Y_R2_1010583',\
			'mosR2v4_Y_R2_1010994',\
			'mosR2v4_Y_R2_1010994',\
			'mosR2v4_Y_R2_1013478',\
			'mosR2v4_Y_R2_1013478',\
			'mosR3v4_Y_R3_772773',\
			'mosR3v4_Y_R3_775548']
plot_list = [1,1,0,0,0,1,0,1,0,1,0,1,1,1]

henry2013 = np.loadtxt('/Users/jhyoon/Dropbox/Research/mosfire/t_exp4obs/t_exp_Henry2013.dat',dtype=str,skiprows=1)
id_henry2013 = henry2013[:,0]
z_henry2013 = np.array(map(float,henry2013[:,3]))
pre_target_name = '9999999'

Ha1_list = np.array([10601,11008,11016,10747,10720,10730,10601,10956])
Ha2_list = np.array([10611,11025,11027,10763,10730,10740,10613,10968])
NII1_list = Ha1_list + 34
NII2_list = Ha2_list + 34
NII1_list[0] = 10638
NII2_list[0] = 10642
NII1_list[2] = 11056
NII2_list[2] = 11060
NII1_list[7] = 10993
NII2_list[7] = 10998
count = -1

for j in range(len(flist)):
	fname = flist[j]
	target_name = fname[13:]
	print j,target_name
	data = np.loadtxt(fname+'_1dspec.fits',dtype=float)
	wavelength = np.array(data[:,0])

	wavelength = np.array(data[:,0])
	spec1d = np.array(data[:,1])
	var1d = np.array(data[:,2])**2.

	sub_match = list(id_henry2013).index(target_name)
	z = z_henry2013[sub_match]
	lambda_Ha = 6563.*(1.+z)
	lambda_NII1 = 6548.*(1.+z)
	lambda_NII2 = 6584.*(1.+z)

	if (pre_target_name != target_name):
		spec1d_coadd = spec1d
		var1d_coadd = var1d
	if (pre_target_name == target_name):
		spec1d_coadd += spec1d
		var1d_coadd += var1d

	if plot_list[j]==1:
		plt.clf()
		plt.plot(wavelength,spec1d_coadd,color='black',drawstyle='steps-mid')
		plt.suptitle(fname)
		plt.plot(wavelength,np.sqrt(var1d_coadd),color='red',drawstyle='steps-mid')
		xmin = 9750
		xmax = 11250
		ymin = min(spec1d_coadd)
		ymax = max(spec1d_coadd)
		plt.xlim(xmin,xmax)
		plt.ylim(ymin,ymax)
		plt.axvline(lambda_Ha,color='black',linestyle='--')
		plt.axvline(lambda_NII1,color='black',linestyle='--')
		plt.axvline(lambda_NII2,color='black',linestyle='--')
		fig_pdf.savefig()
		plt.clf()
		plt.plot(wavelength,spec1d_coadd,color='black',drawstyle='steps-mid')
		#plt.plot(wavelength,spec1d_weight,color='green',linestyle='-',drawstyle='steps-mid')
		plt.plot(wavelength,np.sqrt(var1d_coadd),color='red',drawstyle='steps-mid')
		plt.suptitle(fname)
		plt.axvline(lambda_Ha,color='black',linestyle='--')
		plt.axvline(lambda_NII1,color='black',linestyle='--')
		plt.axvline(lambda_NII2,color='black',linestyle='--')
		xmin = lambda_NII1-30
		xmax = lambda_NII2+30
#		xmin = lambda_NII1-150
#		xmax = lambda_NII2+150
		plt.xlim(xmin,xmax)
		plt.ylim(ymin,ymax)
		xp1 = (xmax-xmin)*0.1 + xmin
		xp2 = (xmax-xmin)*0.7 + xmin
		yp1 = ymax*0.9
		yp2 = ymax*0.8
		yp3 = ymax*0.7
		yp4 = ymax*0.6

		count += 1
		Ha1 = Ha1_list[count]
		Ha2 = Ha2_list[count]
		NII1 = NII1_list[count]
		NII2 = NII2_list[count]
		sub_Ha = (wavelength > Ha1) & (wavelength < Ha2)
		sub_NII = (wavelength > NII1) & (wavelength < NII2)
		Ha_flux = sum(spec1d_coadd[sub_Ha])
		Ha_sigma = np.sqrt(sum(var1d_coadd[sub_Ha]))
		plt.text(xp1,yp1,'Ha=%6.3f'%(Ha_flux))
		plt.text(xp1,yp2,'sigma_Ha=%6.3f'%(Ha_sigma))
		plt.text(xp1,yp3,'S/N(Ha)=%6.3f'%(Ha_flux/Ha_sigma))
		plt.text(xp1,yp4,'Ha*0.05=%6.3f'%(Ha_flux*0.05))
		plt.axvline(Ha1,color='black',linestyle=':')
		plt.axvline(Ha2,color='black',linestyle=':')
		NII_flux = sum(spec1d_coadd[sub_NII])
		NII_sigma = np.sqrt(sum(var1d_coadd[sub_NII]))
		plt.text(xp2,yp1,'NII=%6.3f'%(NII_flux))
		plt.text(xp2,yp2,'sigma_NII=%6.3f'%(NII_sigma))
		plt.text(xp2,yp3,'3*sigma_NII=%6.3f'%(3*NII_sigma))
		plt.text(xp2,yp4,'S/N(NII)=%6.3f'%(NII_flux/NII_sigma))
		plt.axvline(NII1,color='black',linestyle=':')
		plt.axvline(NII2,color='black',linestyle=':')
		if count==0:
			sub_continuum = ((wavelength > 10475) & (wavelength < 10505)) | \
								((wavelength > 10530) & (wavelength < 10570)) | \
								((wavelength > 10615) & (wavelength < 10635)) | \
								((wavelength > 10660) & (wavelength < 10680))
		if count==1:
			sub_continuum = ((wavelength > 10865) & (wavelength < 10895)) | \
								((wavelength > 10905) & (wavelength < 10922)) | \
								((wavelength > 10955) & (wavelength < 10970)) | \
								((wavelength > 11097) & (wavelength < 11140)) | \
								((wavelength > 11160) & (wavelength < 11195))
		if count==2:
			sub_continuum = ((wavelength > 10905) & (wavelength < 10921)) | \
								((wavelength > 10933) & (wavelength < 10945)) | \
								((wavelength > 10955) & (wavelength < 10970)) | \
								((wavelength > 11095) & (wavelength < 11138)) | \
								((wavelength > 11160) & (wavelength < 11195))
		if count==3:
			sub_continuum = ((wavelength > 10595) & (wavelength < 10637)) | \
								((wavelength > 10660) & (wavelength < 10680)) | \
								((wavelength > 10808) & (wavelength < 10830)) | \
								((wavelength > 10905) & (wavelength < 10922))
		if count==4:
			sub_continuum = ((wavelength > 10595) & (wavelength < 10637)) | \
								((wavelength > 10660) & (wavelength < 10680)) | \
								((wavelength > 10780) & (wavelength < 10790)) | \
								((wavelength > 10810) & (wavelength < 10825))
		if count==5:
			sub_continuum = ((wavelength > 10600) & (wavelength < 10640)) | \
								((wavelength > 10660) & (wavelength < 10675)) | \
								((wavelength > 10780) & (wavelength < 10825)) | \
								((wavelength > 10865) & (wavelength < 10872))
		if count==6:
			sub_continuum = ((wavelength > 10455) & (wavelength < 10470)) | \
								((wavelength > 10479) & (wavelength < 10506)) | \
								((wavelength > 10535) & (wavelength < 10570)) | \
								((wavelength > 10615) & (wavelength < 10632)) | \
								((wavelength > 10660) & (wavelength < 10675))
		if count==7:
			sub_continuum = ((wavelength > 10780) & (wavelength < 10825)) | \
								((wavelength > 10863) & (wavelength < 10870)) | \
								((wavelength > 10905) & (wavelength < 10920)) | \
								((wavelength > 10980) & (wavelength < 10992)) | \
								((wavelength > 11035) & (wavelength < 11050)) | \
								((wavelength > 11095) & (wavelength < 11030))
		continuum = spec1d_coadd[sub_continuum]
		rms = np.sqrt(np.mean(continuum**2.))
		mean_sigma = np.sqrt(np.mean(var1d_coadd[sub_continuum]))
		plt.text(xp2,yp4*0.8,'rms=%6.3f,mean sigma=%6.3f'%(rms,mean_sigma))
		plt.plot(wavelength[sub_continuum],spec1d_coadd[sub_continuum],color='cyan',marker='.',linestyle='none')
		fig_pdf.savefig()
	
	pre_target_name = target_name


fig_pdf.close()

