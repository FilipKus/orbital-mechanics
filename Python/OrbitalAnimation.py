import numpy as np
from Spacecraft import Spacecraft as SC

if __name__ == '__main__':

	# Initial conditions of orbit parameters
	# 	Alternative state vectors for testing
	# 	r0 = [ -3670 , -3870 , 4400 ]
	# 	v0 = [ 4.7 , -7.4 , 1 ])

	r_mag = 6378 + 500
	v_mag = np.sqrt(398600/r_mag)
	r0 = [r_mag, 0, 0]
	v0 = [0, v_mag, 0]

	sc = SC(
		{
		'r'				:r0,
		'v'				:v0,
		'tspan'			:"2",
		'dt'			:100.0,			
		} )

	sc.plot()


