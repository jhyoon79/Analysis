#!/usr/bin/env python

###
### conduct polynomial and Gaussian fitting for spectrum data
###
import numpy as np
import matplotlib
import matplotlib.pylab as plt
import math
from scipy.optimize import curve_fit
from astropy.io import fits
from astropy.io import ascii
import pdb
import subroutine_flux_calib as subroutine

def gauss_func(x,a,mu,sigma):
	return a*np.exp(-(x-mu)**2./(2*sigma**2.))

###
### polynomial fit for continuum masking emission/absorption lines with 3*rms
###
def continuum_fit(lambda_input,flux_input,Ndeg):
	lambda_input_sub = lambda_input
	flux_input_sub = flux_input
	sub = np.arange(len(lambda_input))
	sub2 = [999]
	while sum(sub2) > 0:
		lambda_input_sub = lambda_input_sub[sub]
		flux_input_sub = flux_input_sub[sub]
		fit = np.polyfit(lambda_input_sub,flux_input_sub,Ndeg,full=True)
		rms_residual = np.sqrt(fit[1]/len(lambda_input_sub))
#		print 'rms type',type(rms_residual),'rms size',rms_residual.size
		if rms_residual.size == 0:
			rms_residual = 999999.
		func = np.poly1d(fit[0])
		sub = abs(flux_input_sub-func(lambda_input_sub)) < 3.*rms_residual
		sub2 = abs(flux_input_sub-func(lambda_input_sub)) > 3.*rms_residual
#	plt.plot(lambda_input_sub,flux_input_sub,'o',color='yellow')
	return func 

###############################
### compute emission line flux
###############################
def emission_flux(lambda_input,flux_input,fluxErr_input,lineI,lineF):
	sub_line = (lambda_input > lineI-6) & (lambda_input < lineF+6)

	dw_cont = 40
	sub_line_cont = ((lambda_input > lineI-dw_cont) & (lambda_input < lineI)) | \
					((lambda_input > lineF) & (lambda_input < lineF+dw_cont)) 
	sub_line_fit = np.copy(sub_line)
	p0 = [max(flux_input[sub_line]),(lineI+lineF)/2.,(lineF-lineI)/4.]
	if (galaxy_name=='p1_1632861') & (p0[1]==6563): 
		sub_line_fit = (lambda_input > 6552) & (lambda_input < 6562) 
	if (galaxy_name=='p1_1633625') & (p0[1]==6563): 
		sub_line_fit = (lambda_input > 6557) & (lambda_input < 6563) 
	if (galaxy_name=='p3_1629789') & (p0[1]==6563): 
		sub_line_fit = (lambda_input > 6558) & (lambda_input < 6566) 
	if (galaxy_name=='p1_1633414') & (p0[1]==6584): 
		sub_line_fit = (lambda_input > 6580) & (lambda_input < 6587) 
	if (galaxy_name=='p1_1633414') & (p0[1]==6548): 
		sub_line_err = (lambda_input > 6548-3*sigma_Ha) & \
							(lambda_input < 6548+3*sigma_Ha)
		fluxErr_line = np.sqrt(sum(fluxErr_input[sub_line_err]**2.))
		flux_limit_line = fluxErr_line*3.
		return -999,fluxErr_line,flux_limit_line,6548,sigma_Ha	
	if (galaxy_name=='p2_1622127') & (p0[1]==6584): 
		sub_line_err = (lambda_input > 6584-3*sigma_Ha) & \
							(lambda_input < 6584+3*sigma_Ha)
		fluxErr_line = np.sqrt(sum(fluxErr_input[sub_line_err]**2.))
		flux_limit_line = fluxErr_line*3.
		return -999,fluxErr_line,flux_limit_line,6584,sigma_Ha	
	if (galaxy_name=='R2_1010583') & (p0[1]==6731): 
		sub_line_err = (lambda_input > 6731-3*sigma_Ha) & \
							(lambda_input < 6731+3*sigma_Ha)
		fluxErr_line = np.sqrt(sum(fluxErr_input[sub_line_err]**2.))
		flux_limit_line = fluxErr_line*3.
		return -999,fluxErr_line,flux_limit_line,6731,sigma_Ha	

	popt, pcov = curve_fit(gauss_func,lambda_input[sub_line_fit],\
								flux_input[sub_line_fit],p0=p0)
	yfit = gauss_func(lambda_input[sub_line],*popt)
	line_center = popt[1]
	line_sigma = abs(popt[2])

	print 'line',p0[1],line_center,line_sigma

	plt.plot(lambda_input[sub_line],yfit,color='magenta')
	flux_line = sum(yfit)
	sub_line_err = (lambda_input > line_center-3*line_sigma) & \
						(lambda_input < line_center+3*line_sigma)
	fluxErr_line = np.sqrt(sum(fluxErr_input[sub_line_err]**2.))#sum(yfit-popt[3])
	flux_limit_line = fluxErr_line*3.
	return flux_line,fluxErr_line,flux_limit_line,line_center,line_sigma	

	
def find_redshift(galaxy_name):
	sub_match = [k for k in range(len(name_redshift)) \
					if name_redshift[k][30:]==galaxy_name]
	if galaxy_name=='F15': redshift_Hb = 0.674 
	if galaxy_name=='F17': redshift_Hb = 0.623 
	if galaxy_name=='F31': redshift_Hb = 0.621 
	if galaxy_name=='F32': redshift_Hb = 0.621 
	if galaxy_name[0]=='F': return redshift_Hb

	if len(sub_match) > 0:
		z_spec_Hb = z_spec[sub_match[0]]
		z_phot_Hb = z_phot[sub_match[0]]
	if len(sub_match) == 0:
		sub_match = [k for k in range(len(name_redshift)) \
						if name_redshift[k][28:]==galaxy_name]
		z_spec_Hb = z_spec[sub_match[0]]
		z_phot_Hb = z_phot[sub_match[0]]
	if galaxy_name[3:] == '1015084': redshift_Hb = 0.679
	if galaxy_name[3:] == '1015516': redshift_Hb = 11022/6563.-1
	if z_spec_Hb < 0: redshift_Hb = z_phot_Hb
	if z_spec_Hb >= 0: redshift_Hb = z_spec_Hb

	return redshift_Hb
###
### BODY
# requirement: define_group.py, subroutine_flux_calib.py
###

### read Y-band Vista filter data
global lambda_limit,name_redshift,z_spec,z_phot,galaxy_name,obs_date,\
		sigma_Ha

import commands

data_skyline = np.loadtxt('oh_lines_Y.txt',dtype=float)
w_skyline = data_skyline[:,0]

### read redshift
all_redshift = np.loadtxt('data2014/fname_mosR123_strict.list',\
					skiprows=1,dtype=str)
name_redshift = all_redshift[:,0]
z_spec = np.array(map(float,all_redshift[:,1]))
z_phot = np.array(map(float,all_redshift[:,2]))

redshift_all = []
flux_Ha_all = []
flux_Hb_all = []
flux_NII_1_all = []
flux_NII_2_all = []
flux_SII_1_all = []
flux_SII_2_all = []
fluxErr_Ha_all = []
fluxErr_Hb_all = []
fluxErr_NII_1_all = []
fluxErr_NII_2_all = []
fluxErr_SII_1_all = []
fluxErr_SII_2_all = []
logOH12_all = []
logOH12_limit_all = []
galaxy_name_all = []
NII_2_detection_all = []
SN_mosfire = []

fname = np.loadtxt('spec_calib/fname_all.list',dtype=str)
fname_deimos = np.loadtxt('/Users/jhyoon/Research/deimos/fname_deimos.list',dtype=str)

Ndeg = np.ones(len(fname))*6
Ndeg[8] = 4
Ndeg[12] = 5

i_deimos = 0
ii = 0
count = 1
fontsize = 9
dw = 5
Ha = 6563
nrow = 12#9
ncol = 7#9
OII3729 = 3729
OII = 3728	# for 3727 and 3729
Hb = 4861
OIII4959 = 4959
OIII5007 = 5007
NII_1 = 6548
NII_2 = 6584
SII_1 = 6717
SII_2 = 6731
dw = 10

fig = plt.figure(figsize=(10,8))
matplotlib.rcParams.update({'font.size':10})
matplotlib.rcParams['xtick.major.size'] = 2
matplotlib.rcParams['ytick.major.size'] = 2
flux_unit = '$10^{-18}$ erg s$^{-1}$ cm$^{-2} \AA^{-1}$'
labelsize = 5
ax0 = fig.add_subplot(1,1,1)
ax0.set_xlabel('$\lambda_{\\rm rest} (\AA)$',fontsize=15)
ax0.set_ylabel('flux ('+flux_unit+')',fontsize=15)
ax0.tick_params(labelcolor='w',top='off',bottom='off',left='off',right='off')
ax0.spines['top'].set_visible(False)
ax0.spines['bottom'].set_visible(False)
ax0.spines['left'].set_visible(False)
ax0.spines['right'].set_visible(False)
fig.subplots_adjust(hspace=0.05,wspace=0.25)

while ii < len(fname):
	### read MOSFIRE spectrum
	spectrum = np.loadtxt('spec_calib/'+fname[ii],skiprows=1,dtype=float)
	lambda_spec = spectrum[:,0]
	sep = fname[ii][13:].find('_')
	galaxy_name = fname[ii][10:13+sep]
	sep_date = fname[ii].find('201')
	obs_date = fname[ii][sep_date:sep_date+9]

	sub_overlap = []
	for j in range(len(fname)):
		sep2 = fname[j][13:].find('_')
		if fname[j][10:13+sep2]==galaxy_name:
			sub_overlap.append(j)

	if fname[ii][:6] == 'mosFv2':
		sep3 = fname[ii][9:].find('_')
		galaxy_name = fname[ii][9:9+sep3]

	if galaxy_name=='1583155': break
	galaxy_name_all.append(galaxy_name)
	print '### gal name ',ii,galaxy_name,obs_date


	###
	### coadd spectra of objects observed multiple times
	###
	flux_spec = 0
	fluxVar_spec = 0
	for k in sub_overlap:
		print '---overlap ',k
		spectrum2 = np.loadtxt('spec_calib/'+fname[k],skiprows=1,dtype=float)
#	spectrum2 = np.loadtxt('spec_calib/'+fname[ii],skiprows=1,dtype=float)
		flux_spec += spectrum2[:,1]
		fluxVar_spec += spectrum2[:,2]**2
	
	f_norm = 1e-18
	fluxErr_spec = np.sqrt(fluxVar_spec)/f_norm
	flux_spec = flux_spec/f_norm

	ii = ii + len(sub_overlap)
#	ii += 1

	###
	### find redshift measured for Hb or photo-z
	###
	redshift_Hb = find_redshift(galaxy_name)
	is_z_phot = False
	if redshift_Hb < 0: is_z_phot = True
	redshift_Hb = abs(redshift_Hb)

	###
	### read deimos spectrum
	###
	data_deimos = np.loadtxt(fname_deimos[i_deimos],dtype=float,skiprows=1)
	i_deimos += 1
	lambda_deimos = data_deimos[:,0]
	c = 2.99792e8/1e-10	# speed of light, Angstrom/sec
	flux_deimos = data_deimos[:,1]/(lambda_deimos**2./c)/f_norm
			#  convert f_nu to f_lambda
	fluxErr_deimos = data_deimos[:,3]/(lambda_deimos**2./c)/f_norm
	lambda_deimos = lambda_deimos/(1.+redshift_Hb)

	###
	### convert lambda to restframe
	###
	Ha_obs = 6563*(1.+redshift_Hb)
	lambda_spec = lambda_spec/(1.+redshift_Hb)
	print 'lambda_max',max(lambda_spec)

	###
	### zoom in plot at restframe, show line fluxes
	###
	fluxErr_deimos = fluxErr_deimos - 2.
	fluxErr_spec = fluxErr_spec - 2.

	sub_OII = abs(lambda_deimos-OII) < dw
	sub_Hb = abs(lambda_deimos-Hb) < dw
	sub_OIII4959 = abs(lambda_deimos-OIII4959) < dw
	sub_OIII5007 = abs(lambda_deimos-OIII5007) < dw
	ymin_OII = min(flux_deimos[sub_OII])*1.2
	ymax_OII = max(flux_deimos[sub_OII])*1.1
	ymin_Hb = min(list(flux_deimos[sub_OII])+\
				list(flux_deimos[sub_Hb])+\
				list(flux_deimos[sub_OIII4959])+\
				list(fluxErr_deimos[sub_Hb])+\
				list(fluxErr_deimos[sub_Hb])+\
				list(fluxErr_deimos[sub_OIII4959]))*1.2
	ymax_Hb = max(list(flux_deimos[sub_OII])+\
				list(flux_deimos[sub_Hb])+\
				list(flux_deimos[sub_OIII4959])+\
				list(fluxErr_deimos[sub_Hb])+\
				list(fluxErr_deimos[sub_Hb])+\
				list(fluxErr_deimos[sub_OIII4959]))*1.1
	ymin_OIII5007 = min(fluxErr_deimos[sub_OIII5007])*1.2
	ymax_OIII5007 = max(flux_deimos[sub_OIII5007])*1.1

	sub_Ha = abs(lambda_spec-Ha) < dw
	ymin_Ha = min(fluxErr_spec[sub_Ha])*1.2
	ymax_Ha = max(flux_spec[sub_Ha])*1.1
	sub_NII_1 = abs(lambda_spec-NII_1) < dw
	sub_NII_2 = abs(lambda_spec-NII_2) < dw
	ymin_NII = min(list(flux_spec[sub_NII_1])+\
				list(flux_spec[sub_NII_2])+\
				list(fluxErr_spec[sub_NII_1])+\
				list(fluxErr_spec[sub_NII_2]))*1.2
	ymax_NII = max(list(flux_spec[sub_NII_1])+\
				list(flux_spec[sub_NII_2])+\
				list(fluxErr_spec[sub_NII_1])+\
				list(fluxErr_spec[sub_NII_2]))*1.1
	sub_SII_1 = abs(lambda_spec-SII_1) < dw
	sub_SII_2 = abs(lambda_spec-SII_2) < dw

	if sum(sub_SII_1) == 0: sub_SII_1 = [-1]
	if sum(sub_SII_2) == 0: sub_SII_2 = [-1]
	ymin_SII = min(list(flux_spec[sub_SII_1])+\
				list(flux_spec[sub_SII_2]))*1.2
	ymax_SII = max(list(flux_spec[sub_SII_1])+\
				list(flux_spec[sub_SII_2]))*1.1

	if count==ncol*1+1: ymax_NII = 3

	ypos_Ha = (ymax_Ha-ymin_Ha)/5.+ymin_Ha
	ypos_NII = (ymax_NII-ymin_NII)/5.+ymin_NII
	ypos_SII = (ymax_SII-ymin_SII)/5.+ymin_SII

	labelbottom = 'off'
	if galaxy_name=='F32': labelbottom = 'on'
	ax1 = fig.add_subplot(nrow,ncol,count)
	plt.plot(lambda_deimos,flux_deimos,color='black',drawstyle='steps-mid')
	plt.plot(lambda_deimos,fluxErr_deimos,color='red',drawstyle='steps-mid')
	for j in range(len(w_skyline)):
		plt.plot([w_skyline[j]/(1+redshift_Hb),\
					w_skyline[j]/(1+redshift_Hb)],\
					[-10.,ypos_NII],color='blue',linestyle='-')
	ax1.tick_params(axis='both',labelsize=labelsize,labelbottom=labelbottom)
	plt.axvline(OII3729,color='orange',linestyle=':')
	plt.xlim(OII-dw,OII+dw)
	plt.ylim(ymin_Hb,ymax_Hb)

	ax2 = fig.add_subplot(nrow,ncol,count+1)
	plt.plot(lambda_deimos,flux_deimos,color='black',drawstyle='steps-mid')
	plt.plot(lambda_deimos,fluxErr_deimos,color='red',drawstyle='steps-mid')
	for j in range(len(w_skyline)):
		plt.plot([w_skyline[j]/(1+redshift_Hb),\
					w_skyline[j]/(1+redshift_Hb)],\
					[-10.,ypos_NII],color='blue',linestyle='-')
	ax2.tick_params(axis='both',labelsize=labelsize,labelbottom=labelbottom)
	plt.axvline(Hb,color='orange',linestyle=':')
	plt.xlim(Hb-dw,Hb+dw)
	plt.ylim(ymin_Hb,ymax_Hb)

	ax3 = fig.add_subplot(nrow,ncol,count+2)
	plt.plot(lambda_deimos,flux_deimos,color='black',drawstyle='steps-mid')
	plt.plot(lambda_deimos,fluxErr_deimos,color='red',drawstyle='steps-mid')
	for j in range(len(w_skyline)):
		plt.plot([w_skyline[j]/(1+redshift_Hb),\
					w_skyline[j]/(1+redshift_Hb)],\
					[-10.,ypos_NII],color='blue',linestyle='-')
	ax3.tick_params(axis='both',labelsize=labelsize,labelbottom=labelbottom)
	plt.axvline(OIII4959,color='orange',linestyle=':')
	plt.xlim(OIII4959-dw,OIII4959+dw)
	plt.ylim(ymin_Hb,ymax_Hb)

	ax4 = fig.add_subplot(nrow,ncol,count+3)
	plt.plot(lambda_deimos,flux_deimos,color='black',drawstyle='steps-mid')
	plt.plot(lambda_deimos,fluxErr_deimos,color='red',drawstyle='steps-mid')
	for j in range(len(w_skyline)):
		plt.plot([w_skyline[j]/(1+redshift_Hb),\
					w_skyline[j]/(1+redshift_Hb)],\
					[-10.,ypos_NII],color='blue',linestyle='-')
	ax4.tick_params(axis='both',labelsize=labelsize,labelbottom=labelbottom)
	plt.axvline(OIII5007,color='orange',linestyle=':')
	plt.xlim(OIII5007-dw,OIII5007+dw)
	plt.ylim(ymin_OIII5007,ymax_OIII5007)

	ax5 = fig.add_subplot(nrow,ncol,count+4)
	plt.plot(lambda_spec,flux_spec,color='black',drawstyle='steps-mid')
	plt.plot(lambda_spec,fluxErr_spec,color='red',drawstyle='steps-mid')
	ax5.tick_params(axis='both',labelsize=labelsize,labelbottom=labelbottom)
	plt.axvline(NII_1,color='orange',linestyle=':')
	plt.xlim(NII_1-dw,NII_1+dw)
	plt.ylim(ymin_NII,ymax_NII)

	ax6 = fig.add_subplot(nrow,ncol,count+5)
	plt.plot(lambda_spec,flux_spec,color='black',drawstyle='steps-mid')
	plt.plot(lambda_spec,fluxErr_spec,color='red',drawstyle='steps-mid')
	ax6.tick_params(axis='both',labelsize=labelsize,labelbottom=labelbottom)
	plt.axvline(Ha,color='orange',linestyle=':')
	plt.xlim(Ha-dw,Ha+dw)
	plt.ylim(ymin_Ha,ymax_Ha)

	ax7 = fig.add_subplot(nrow,ncol,count+6)
	plt.plot(lambda_spec,flux_spec,color='black',drawstyle='steps-mid')
	plt.plot(lambda_spec,fluxErr_spec,color='red',drawstyle='steps-mid')
	ax7.tick_params(axis='both',labelsize=labelsize,labelbottom=labelbottom)
	plt.axvline(NII_2,color='orange',linestyle=':')
	plt.xlim(NII_2-dw,NII_2+dw)
	plt.ylim(ymin_NII,ymax_NII)

	if count == 1:
		ax1.set_title('OII $\lambda\lambda$3727,3729')
		ax2.set_title('H$\\beta$')
		ax3.set_title('OIII $\lambda$4959')
		ax4.set_title('OIII $\lambda$5007')
		ax5.set_title('NII $\lambda$6548')
		ax6.set_title('H$\\alpha$')
		ax7.set_title('NII $\lambda$6584')

	plt.savefig('line_plot.eps')
	plt.savefig('line_plot.pdf')
	count += ncol

	if galaxy_name=='F32': break

