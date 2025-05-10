import math
import numpy as np
from numpy.linalg import norm


def coeFromSV( r, v ):

	mu = 398600
	rmag = norm(r)

	# Angular Momentum
	h = np.cross(r, v)
	hmag = norm(h)


	# Eccentricity
	e = np.cross(v, h)/mu - r/rmag
	emag = norm(e)


	# Inclination
	irad = math.acos(np.dot(h, [0, 0, 1])/hmag)
	ideg = np.rad2deg(irad)


	# Right ascension of the ascending node
	N = np.cross([0, 0, 1], h)
	Nmag = norm(N)

	if Nmag == 0:
		RArad = 0
	else:
		RAtemp = math.acos(np.dot([1, 0, 0], N)/Nmag)
		if np.dot(N, [0, 1, 0]) >= 0:
			RArad = RAtemp
		else:
			RArad = 2*math.pi - RAtemp

	RAdeg = np.rad2deg(RArad)	


	# Argument of Perigee
	if Nmag == 0:
		wdeg = 0
		wrad = 0
	else:
		wtemp = math.acos(np.dot(N, e)/(Nmag*emag))		
		if np.dot(e, [0, 0, 1]) >= 0:
			wrad = wtemp
		else:
			wrad = 2*math.pi - wtemp

	wdeg = np.rad2deg(wrad)		


	# True Anomaly
	if Nmag == 0:
		urad = 0		
	else:
		utemp = math.acos(np.dot(N, r)/(Nmag*rmag))		
		if np.dot(r, [0, 0, 1]) >= 0:
			urad = utemp
		else:
			urad = 2*math.pi - utemp


	TArad = urad - wrad
	if TArad < 0:
		TArad = TArad + 2*math.pi		
	TAdeg = np.rad2deg(TArad)					 	



	return hmag, emag, ideg, RAdeg, wdeg, TAdeg

def ellipticalPeriod( coes, mu ):
	a = ((coes[0]**2)/mu)/(1 - coes[1]**2)
	period = (np.pi*2/np.sqrt(mu))*a**(3/2)	
	return period