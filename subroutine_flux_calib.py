#!/usr/bin/env python
import numpy as np
from scipy.interpolate import interp1d
import pdb

###
### lambda interpolation
###
def lambda_interpol(lambda_input,lambda_set,throughput_set,min_lambda,max_lambda):
	sub_select = (lambda_input > min_lambda) & \
					(lambda_input < max_lambda)
	lambda_select = lambda_input[sub_select]
	f = interp1d(lambda_set,throughput_set,kind='linear')
	throughput_filter_new = f(lambda_select)
	return throughput_filter_new,sub_select

###
### read filter info.
###
def read_filter(filter_name):
	Y_filter = np.loadtxt('%s.res'%filter_name,dtype=float)
	lambda_filter = Y_filter[:,0]
	throughput_filter = Y_filter[:,1]
	min_lambda_filter = min(lambda_filter)
	max_lambda_filter = max(lambda_filter)
	return lambda_filter,throughput_filter,min_lambda_filter,max_lambda_filter

###
### convolve filter throughput, see Yi et al. 1995 PASP
###
def convolve_filter(lambda_input,flux_input,fluxErr_input,filter_name,below9750=False,optional_return=False):
	lambda_filter, throughput_filter, min_lambda_filter, max_lambda_filter = read_filter(filter_name)
	throughput_filter_new,sub_select = lambda_interpol(lambda_input,lambda_filter,throughput_filter,min_lambda_filter,max_lambda_filter)
#	if below9750==False:
#		sub_filter = lambda_input[sub_select] < 9750
#		throughput_filter_new[sub_filter] = 0

	lambda_select = lambda_input[sub_select]
	flux_select = flux_input[sub_select]
	fluxErr_select = fluxErr_input[sub_select]
#	subsub = (fluxErr_select >= 0) & (np.isfinite(flux_select) == True)
	subsub = np.isfinite(flux_select)
	

	#dlambda_input = np.mean(lambda_select[1:]-lambda_select[:-1])
	dlambda_input = lambda_select[1:]-lambda_select[:-1]
	tmp = dlambda_input.tolist()
	tmp.extend([tmp[0]])
	dlambda_input = np.array(tmp)
	numer = np.sum(flux_select[subsub]*throughput_filter_new[subsub]* \
						lambda_select[subsub]*dlambda_input[subsub])
	numer_err = np.sqrt(np.sum((fluxErr_select[subsub]* \
					throughput_filter_new[subsub]* \
					lambda_select[subsub]*dlambda_input[subsub])**2.))
	denom = np.sum(throughput_filter_new[subsub]*lambda_select[subsub]*dlambda_input[subsub])
	f_filter = numer/denom
	fErr_filter = numer_err/denom
	m_filter = -2.5*np.log10(f_filter)
	if optional_return==True:
		return f_filter,fErr_filter,m_filter,lambda_filter,throughput_filter
	return f_filter,fErr_filter,m_filter

###
### do relative flux calibration with the sensivity function
###
def flux_calib(lambda_input,flux_input,fluxErr_input):
#	sens_func = np.loadtxt('sensfunc_mean_Joo.dat',dtype=float,skiprows=1)
	sens_func = np.loadtxt('sensfunc_2013jan04_mean.dat',dtype=float,skiprows=1)
	lambda_sfunc = sens_func[:,0]
	flux_sfunc = sens_func[:,1]
	lambda_sub = lambda_input[2:-2]
	flux_sub = flux_input[2:-2]
	fluxErr_sub = fluxErr_input[2:-2]
	return lambda_sub,flux_sub*flux_sfunc,fluxErr_sub*flux_sfunc

def flux_calib_2014dec29(lambda_input,flux_input,fluxErr_input):
	sens_func = np.loadtxt('sensfunc_2014dec29_mean.dat',dtype=float,skiprows=1)
#	sens_func = np.loadtxt('sensfunc_2014dec29_mean_fit.dat',dtype=float,skiprows=0)
	lambda_sfunc = sens_func[:,0]
	flux_sfunc = sens_func[:,1]
	lambda_sub = lambda_input[2:-2]
	flux_sub = flux_input[2:-2]
	fluxErr_sub = fluxErr_input[2:-2]
	return lambda_sub,flux_sub*flux_sfunc,fluxErr_sub*flux_sfunc

def flux_calib_2014dec30(lambda_input,flux_input,fluxErr_input):
	sens_func = np.loadtxt('sensfunc_2014dec30_mean.dat',dtype=float,skiprows=1)
#	sens_func = np.loadtxt('sensfunc_2014dec29_mean_fit.dat',dtype=float,skiprows=0)
	lambda_sfunc = sens_func[:,0]
	flux_sfunc = sens_func[:,1]
	lambda_sub = lambda_input[2:-2]
	flux_sub = flux_input[2:-2]
	fluxErr_sub = fluxErr_input[2:-2]
	return lambda_sub,flux_sub*flux_sfunc,fluxErr_sub*flux_sfunc



