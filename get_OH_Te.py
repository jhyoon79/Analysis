import numpy as np
import matplotlib.pylab as plt
import pdb
from matplotlib.backends.backend_pdf import PdfPages
fig_pdf = PdfPages('get_OH_Te.pdf')

###
### compute 12 + log(O/H) using emission lines
###
def get_OH(OIII4363,OIII4959,OIII5007,Hb):
	Te = np.arange(5000,20000,1)
	t3 = Te/1e4
	ne = 100	# cm^-3
	x = 1e-4*ne*t3**(-0.5)
	C_T = (8.44-1.09*t3+0.5*t3**2.-0.08*t3**3.)*(1.+0.0004*x)/(1.+0.044*x)
	log_OIII_ratio = 1.432/t3+np.log10(C_T)
	log_OIII_ratio_obs = np.log10((OIII4959+OIII5007)/OIII4363)

	Te_obs = []
	plt.clf()
	plt.plot(Te,log_OIII_ratio,color='black',marker='.',linestyle='none')
	plt.yscale('log')
	for i in range(len(OIII4363)):
		plt.axhline(log_OIII_ratio_obs[i],linestyle='--')
		d_ratio = abs(log_OIII_ratio_obs[i]-log_OIII_ratio)
		min_d_ratio = min(d_ratio)
		min_sub = list(d_ratio).index(min_d_ratio)
		Te_obs.append(Te[min_sub])
	plt.xlim(10000,20000)
	plt.ylim(1.,3)
	plt.xlabel('Te')
	plt.ylabel('(OIII4959+5007)/OIII4363')

	Te_obs = np.array(Te_obs)
	t3_obs = Te_obs/1e4
	logOIIIH = np.log10((OIII4959+OIII5007)/Hb)+6.200+1.251+1.251/t3_obs - \
				5*np.log10(t3_obs)-0.014*t3_obs
	t2_obs = -0.577+t3*(2.065-0.498*t3)
	logOIIH = np.log10(OII3727/Hb)+5.961+1.676/t2_obs-0.4*np.log10(t2_obs) - \
				0.034*t2_obs+np.log10(1+1.35*x)
	OH = 10**(logOIIIH-12.)+10**(logOIIIH-12.)
	logOH = 12 + np.log10(OH)

	return Te_obs,logOIIH,logOIIIH,logOH

###
### read MOSFIRE data
###
emi = np.loadtxt('emission.dat',skiprows=1,dtype=str)
emi = emi[:12,:]
id = emi[:,0]
z = np.array(map(float,emi[:,1]))
Ha = np.array(map(float,emi[:,2]))
Ha_err = np.array(map(float,emi[:,3]))
OII3727 = np.array(map(float,emi[:,12]))
OII3727_err = np.array(map(float,emi[:,13]))
Hb = np.array(map(float,emi[:,14]))
Hb_err = np.array(map(float,emi[:,15]))
OIII4959 = np.array(map(float,emi[:,16]))
OIII4959_err = np.array(map(float,emi[:,17]))
OIII5007 = np.array(map(float,emi[:,18]))
OIII5007_err = np.array(map(float,emi[:,19]))

###
### read DEIMOS data
###
emi_OIII4363 = np.loadtxt('emission_OIII4363.dat',skiprows=1,dtype=str)
OIII4363 = np.array(map(float,emi_OIII4363[:,10]))
OIII4363_err = np.array(map(float,emi_OIII4363[:,11]))

###
### read stacked data
###
emi_stack = np.loadtxt('emission_stack.dat',skiprows=1,dtype=str)
Hb_stack = np.array(map(float,emi_stack[:,1]))
Hb_err_stack = np.array(map(float,emi_stack[:,2]))
OII3727_stack = np.array(map(float,emi_stack[:,4]))
OII3727_err_stack = np.array(map(float,emi_stack[:,5]))
OIII4363_stack = np.array(map(float,emi_stack[:,6]))
OIII4363_err_stack = np.array(map(float,emi_stack[:,7]))
OIII4959_stack = np.array(map(float,emi_stack[:,8]))
OIII4959_err_stack = np.array(map(float,emi_stack[:,9]))
OIII5007_stack = np.array(map(float,emi_stack[:,10]))
OIII5007_err_stack = np.array(map(float,emi_stack[:,11]))
OIII4363 = np.array(map(float,emi_OIII4363[:,10]))
OIII4363_err = np.array(map(float,emi_OIII4363[:,11]))

###
### compute 12+log(O/H)
###
Te,logOIIH,logOIIIH,logOH = get_OH(OIII4363,OIII4959,OIII5007,Hb)
fig_pdf.savefig()
fig_pdf.close()

