import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
fig_pdf = PdfPages('error_check.pdf')
import matplotlib.pylab as plt

flist = ['mosR1v4_Y_R1_1009808_2014dec29_1dspec.dat',\
		'mosR1v4_Y_R1_1015084_2014dec29_1dspec.dat',\
		'mosR1v4_Y_R2_1015516_2014dec29_1dspec.dat',\
		'mosR2v4_Y_R2_1015516_2014dec29_1dspec.dat',\
		'mosR2v4_Y_R2_1015516_2014dec30_1dspec.dat',\
		'mosR3v4_Y_R2_1015516_2014dec30_1dspec.dat',\
		'mosR2v4_Y_R2_1010583_2014dec29_1dspec.dat',\
		'mosR2v4_Y_R2_1010583_2014dec30_1dspec.dat',\
		'mosR2v4_Y_R2_1010994_2014dec29_1dspec.dat',\
		'mosR2v4_Y_R2_1010994_2014dec30_1dspec.dat',\
#		'mosR2v4_Y_R2_1013478_2014dec29_1dspec.dat',\
#		'mosR2v4_Y_R2_1013478_2014dec30_1dspec.dat',\
		'mosR3v4_Y_R3_772773_2014dec30_1dspec.dat',\
		'mosR3v4_Y_R3_775548_2014dec30_1dspec.dat']
redshift = [0.616,0.678,0.679,0.679,0.679,0.679,0.639,0.639,0.634,0.634,0.616,0.670]

data_skyline = np.loadtxt('../oh_lines_Y.txt',dtype=float)
w_skyline = data_skyline[:,0]

rms_all = []
m_e1_all = []
m_e2_all = []

for i in range(len(flist)):
	spectrum = np.loadtxt(flist[i],dtype=float)
	lambda_spec = np.array(spectrum[:,0])
	flux_spec = np.array(spectrum[:,1])
	fluxErr_spec = np.array(spectrum[:,2])
	flux_spec2 = np.array(spectrum[:,3])
	fluxErr_spec2 = np.array(spectrum[:,4])

	Ha_obs = 6563*(1.+redshift[i])
	lambda_min = Ha_obs - 160
	lambda_max = Ha_obs + 160
	sub_continuum = (lambda_spec > lambda_min) & \
						(lambda_spec < lambda_max)
	lambda_continuum = list(lambda_spec[sub_continuum])
	flux_continuum = list(flux_spec[sub_continuum])
	fluxErr_continuum = list(fluxErr_spec[sub_continuum])
	fluxErr_continuum2 = list(fluxErr_spec2[sub_continuum])
	for j in range(len(w_skyline)):
		sub_near_skyline = abs(np.array(lambda_continuum)-w_skyline[j]) < 7
		if sum(sub_near_skyline) > 0:	
			lambda_near_skyline = np.array(lambda_continuum)[sub_near_skyline]
			for k in range(len(lambda_near_skyline)):
				sub_index = lambda_continuum.index(lambda_near_skyline[k])
				lambda_continuum.pop(sub_index)
				flux_continuum.pop(sub_index)
				fluxErr_continuum.pop(sub_index)
				fluxErr_continuum2.pop(sub_index)

	lambda_continuum = np.array(lambda_continuum)
	flux_continuum = np.array(flux_continuum)
	fluxErr_continuum = np.array(fluxErr_continuum)
	fluxErr_continuum2 = np.array(fluxErr_continuum2)
	sub_no_emission = (abs(lambda_continuum/(1+redshift[i])-6563) > 7) & \
							(abs(lambda_continuum/(1+redshift[i])-6548) > 7) & \
							(abs(lambda_continuum/(1+redshift[i])-6584) > 7) & \
							(abs(lambda_continuum/(1+redshift[i])-6717) > 7) & \
							(abs(lambda_continuum/(1+redshift[i])-6731) > 7)
	lambda_continuum = lambda_continuum[sub_no_emission]
	flux_continuum = flux_continuum[sub_no_emission]
	fluxErr_continuum = fluxErr_continuum[sub_no_emission]
	fluxErr_continuum2 = fluxErr_continuum2[sub_no_emission]
	mean_continuum = np.mean(flux_continuum)

	plt.clf()
	plt.plot(lambda_spec,flux_spec,drawstyle='steps-mid')
	plt.plot(lambda_spec,fluxErr_spec,drawstyle='steps-mid',color='red')
	plt.plot(lambda_spec,fluxErr_spec2,drawstyle='steps-mid',color='orange')
	plt.plot(lambda_continuum,flux_continuum,marker='+',ms=5,mew=2,\
				linestyle='none')
	plt.suptitle(flist[i][10:20])	

	rms = (np.sqrt(np.mean(flux_continuum**2.)))
	m_e1 = (np.sqrt(np.mean(fluxErr_continuum**2.)))
	m_e2 = (np.sqrt(np.mean(fluxErr_continuum2**2.)))
	rms_all.append(rms)
	m_e1_all.append(m_e1)
	m_e2_all.append(m_e2)
	sub_yp = abs(lambda_spec-Ha_obs) < 5
	yp = max(flux_spec[sub_yp])
	plt.text(Ha_obs,yp*0.8,'rms=%f'%rms,clip_on=True)
	plt.text(Ha_obs,yp*0.6,'m_e1=%f'%m_e1,clip_on=True)
	plt.text(Ha_obs,yp*0.4,'m_e2=%f'%m_e2,clip_on=True)
	plt.xlim(lambda_min,lambda_max)
	plt.ylim(min(flux_spec),yp)
	fig_pdf.savefig()

plt.clf()
xx = range(len(flist))
plt.plot(xx,rms_all,marker='o',label='rms')
plt.plot(xx,m_e1_all,marker='d',color='red',label='m_e1')
plt.plot(xx,m_e2_all,marker='^',color='orange',label='m_e2')
plt.xlim(-1,12)
plt.ylim(0,0.06)
plt.legend(loc=4,ncol=1,fontsize=20,numpoints=1)
fig_pdf.savefig()

fig_pdf.close()


