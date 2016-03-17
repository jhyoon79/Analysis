#!/usr/bin/env python

import numpy as np
from astropy.io import fits
from astropy.io import ascii
import pdb

A0_list = np.loadtxt('fname_std_iraf.list',dtype=str,skiprows=0)
#A0_list = np.loadtxt('fname_2014dec29.list',dtype=str,skiprows=0)
A0_2014dec29 = A0_list[:10]
fits_list = [129,130,210,211,217,218,227,228,320,321]
fits_list = map(str,fits_list)
name = []
UT = []
airmass = []
for i in range(len(fits_list)):
	fname = '/Users/jhyoon/Tools/MosfireDRP-1.0/Data_20141229/m141229_0'+\
				fits_list[i]+'.fits'
	A0_fits = fits.open(fname)
	fits_header = A0_fits[0].header
	UT.append(fits_header[-1])
	airmass.append(fits_header[-2])
	name.append(A0_2014dec29[i][12:])

ascii.write([name,airmass,UT],'airmass_2014dec29.dat',names=['#  name','   airmass','    UT'])


A0_2014dec30 = A0_list[10:20]
fits_list = [138,139,241,242,243,244,251,252,350,351]
fits_list = map(str,fits_list)
name = []
UT = []
airmass = []
for i in range(len(fits_list)):
	fname = '/Users/jhyoon/Tools/MosfireDRP-1.0/Data_20141230/m141230_0'+\
				fits_list[i]+'.fits'
	A0_fits = fits.open(fname)
	fits_header = A0_fits[0].header
	UT.append(fits_header[-1])
	airmass.append(fits_header[-2])
	name.append(A0_2014dec30[i][12:])

ascii.write([name,airmass,UT],'airmass_2014dec30.dat',names=['#  name','   airmass','    UT'])

A0_2013jan04 = A0_list[20:]
fits_list = [284,285,286,287,320,321,322,323,
				328,329,365,366,367,368,369,370,
				371,372,373,374,375,376]
fits_list = map(str,fits_list)
name = []
UT = []
airmass = []
for i in range(len(fits_list)):
	fname = '/Users/jhyoon/Tools/MosfireDRP-1.0/Data_20130104/m130104_0'+\
				fits_list[i]+'.fits'
	A0_fits = fits.open(fname)
	fits_header = A0_fits[0].header
	UT.append(fits_header[-1])
	airmass.append(fits_header[-2])
	name.append(A0_2013jan04[i][12:])

ascii.write([name,airmass,UT],'airmass_2013jan04.dat',names=['#  name','   airmass','    UT'])



