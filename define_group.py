#!/usr/bin/env python
#################################################
### classfy mosFv2 targets if it has a slit loss
#################################################
def define_group():
	# star
	gr0 = ['0005']
	gr0 = [4]
	# no mosfire slit loss, some Y-band slit loss
	gr1 = ['0001','0003','0008','0014','0015','0017']
	gr1 = [0,2,7,10,11,12]
	# no mosfire slit loss, poor Y-band image
	gr2 = []
	gr2 = []
	# a galaxy is larger than a slit
	gr3 = ['0002','0007','0009','0013','0019','0020','0021','0022']
	gr3 = [1,6,8,9,13,14,15,16]
	# two galaxies
	gr4 = ['0004','0006']
	gr4 = [3,5]
	return gr0,gr1,gr2,gr3,gr4

def define_group_2014dec():
	# star
	gr0 = []
	gr0 = []
	# no mosfire slit loss, some Y-band slit loss
	gr1 = ['1010583','1015516','1015516']
	gr1 = [3,2,6]
	# no mosfire slit loss, poor Y-band image
	gr2 = []
	gr2 = []
	# a galaxy is larger than a slit
	gr3 = [0,1,4,5]
	# two galaxies
	gr4 = []
	return gr0,gr1,gr2,gr3,gr4

