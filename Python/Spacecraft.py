import matplotlib.pyplot as plt
import numpy as np
import Calculations as calc
import math

from numpy.linalg import norm
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import ode
from matplotlib.animation import FuncAnimation

plt.style.use('dark_background')

def init_config():
	return {
		'cbRadius'			: 6378,
		'mu'				: 398600,
		'tspan'				: 10*60,
		'dt'				: 100.0,
		'n_steps'			: [], 
		'r'					: [],
		'v'					: [],
		'rs'				: [],
		'coes'				: [],
		'propagator'		: 'lsoda',		
	}

class Spacecraft:

	def __init__( self, config ):

		# Presetting the orbital configuration
		self.config = init_config()
		# Updating the orbital configuration based on user input
		for key in config.keys():
			self.config[ key ] = config[ key ]		

		# Calculate coes from given state vectors
		if self.config [ 'r' and 'v' ]:
			self.config [ 'coes' ] = calc.coeFromSV( self.config[ 'r' ], self.config[ 'v' ])

		# Interpret 'tspan' as number of orbital periods if user input is a string
		if type(self.config[ 'tspan' ]) ==  str:
			self.config[ 'tspan' ] = int(self.config[ 'tspan' ]) * calc.ellipticalPeriod( self.config[ 'coes' ], self.config[ 'mu' ])

		# Calculate number of steps for ODE solver
		self.config[ 'n_steps' ] = int(np.ceil(self.config[ 'tspan' ] / self.config[ 'dt' ]))
		# Concatenate state vectors into one array
		self.config[ 'rs' ] = self.config[ 'r' ] + self.config[ 'v' ]

	def diffy_q( self, t, y, mu ):
		rx, ry, rz, vx, vy, vz = y
		r = np.array([rx, ry, rz])

		norm_r = np.linalg.norm(r)

		ax, ay, az = -r*self.config[ 'mu' ]/norm_r**3

		return [vx, vy, vz, ax, ay, az]		

	def propagator( self ):
		# Initialise arrays
		ys = np.zeros((self.config[ 'n_steps' ], 6))
		ts = np.zeros((self.config[ 'n_steps' ], 1))

		# Initialise conditions
		y0 = self.config[ 'r' ] + self.config[ 'v' ]		
		ys[0] = np.array(y0)
		step = 1

		# Initiate solver
		solver = ode(self.diffy_q)
		solver.set_integrator(self.config[ 'propagator' ])
		solver.set_initial_value(y0, 0)
		solver.set_f_params(self.config[ 'mu' ])

		# Propagate orbit
		while solver.successful() and step < self.config[ 'n_steps' ]:
			solver.integrate(solver.t + self.config[ 'dt' ])
			ts[step] = solver.t
			ys[step] = solver.y
			step += 1

		rs = ys[:, :3]	
		return rs
				
	def plot( self ):		
		
		rs = self.propagator()		

		# Set up figure
		fig = plt.figure(figsize = (18, 6), )			
		ax = fig.add_subplot(111, projection='3d')	
		ax.xaxis.pane.fill = False
		ax.yaxis.pane.fill = False
		ax.zaxis.pane.fill = False
		
		
		# Plot central body
		_u,_v = np.mgrid[0:2*math.pi:20j, 0:math.pi:10j]
		_x = self.config[ 'cbRadius' ]*np.cos(_u)*np.sin(_v)
		_y = self.config[ 'cbRadius' ]*np.sin(_u)*np.sin(_v)
		_z = self.config[ 'cbRadius' ]*np.cos(_v)
		ax.plot_surface(_x, _y, _z, cmap = 'Blues', alpha = 0.5)

		# Plot the x, y, z vectors
		l = self.config[ 'cbRadius' ]*2
		x, y, z = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
		u, v, w = [[l, 0, 0], [0, l, 0], [0, 0, l]]
		ax.quiver(x, y, z, u, v, w, color = 'w')

		# Labelling and setting limits to axes
		max_val = np.max(np.abs(rs))
		ax.set_xlim([-max_val, max_val])
		ax.set_ylim([-max_val, max_val])
		ax.set_zlim([-max_val, max_val])

		ax.set_xlabel(['X (km)'])
		ax.set_ylabel(['Y (km)'])
		ax.set_zlabel(['Z (km)'])	


		# Initialising the points of the spacecraft (particle) and its path (traj)
		particle, = plt.plot([], [], [], marker = 'o', color = 'r')
		traj, = plt.plot([], [], [], color = 'r', alpha = 1.0)

		# Function that will be called to update animation
		def animate(frame):
			# update plot
			particle.set_data( rs[frame, 0] , rs[frame, 1] )
			particle.set_3d_properties( rs[frame, 2] )
			traj.set_data( rs[:frame, 0], rs[:frame, 1] )
			traj.set_3d_properties( rs[:frame, 2] )

			return particle, traj

		ani = FuncAnimation(fig, animate, frames = range(self.config[ 'n_steps' ]), interval = 50)

		ax.set_aspect('auto')
		ax.set_title('Orbital Animation')				
		plt.show()	
