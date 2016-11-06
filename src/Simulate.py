"""
Jonathan Curtis
11/05/2016

This program is meant to simulate the dynamics of a driven/dissipative optical cavity mode subject to a non-linearity.

The method of simulation is the Monte-Carlo Wavefunction method (AKA quantum-jumps method). An excellent exposition of this method can be found in 

Monte Carlo wave-function method in quantum optics. Molmer, K. Castin, Y. Dalibard, J. Journal of the Optical Society of America B, 10 (3) 1993. DOI 10.1364/JOSAB.10.000524

We will assume we can describe our cavity's may be described by a single quantum mode of the electromagnetic field. 
We truncate the photon Hilbert space at Nmax photons which we take to be much larger than 1. 
Thus, we can describe the system with a single vector in C^{Nmax} where C is the complex number space. Throughout we take hbar = 1. 

We take our dynamics to have three main components. These are the coherent dynamics of the photon Hamiltonian (in a rotating frame), given by 

Hcoh = Delta*N + U(N)

where Delta = Cavity frequency - Drive frequency

N is the number operator of the photon mode. 
U is a non-linearity that is assumed to be approximated well by a simple function of the number operator. 
The second component is a coherent drive mechanism given by a drive potential (in the corotating frame) of 

Hdrive = Omega*(a- + a+)/sqrt(2)

where Omega>0 is the Rabi frequency of drive and a-,a+ are the lowering/raising operators of the photon field respectively. 
The final component is the single-photon loss mechanism given by quantum jump operator of 

C = sqrt(kappa)*a-

where kappa is the cavity decay rate."""


"""
As per the quantum-jump method, this can be simulated by the non-Hermitian Hamiltonian of 

H = Hcoh + Hdrive -ikappa N/2 

where i is the imaginary number. 
Starting with an initial wavefunction of state(0), we first evolve a time dt>0 by applying the operator 

(1 - iHdt)

This defines a new state of 

state(t + dt) = (1-iHdt)state(t)

which has a norm of 

1-dp

with (in our case)

dp = dt <state(t)|kappa N | state(t)>

which is proportional to the decay rate, time, and the expected number of photons in the state. We require dp <<1. 
After evolving by the Hamiltonian we then perform a probabilistic update of the state. With probability 1-dp we stay in state(t+dt).
We must then normalize this state so that after the jump operator, it has norm 1.

With probability dp we jump to state 

state(t+dt) -> sqrt(kappa) a- |state(t+dt)> / sqrt(dp/dt)

We then repeat until the desired precision is achieved. We then repeat this process to obtain an ensemble of states from which we can compute the ensemble/quantum average values of quantities. 
"""

"""
Parameter descriptions.

Nmax = max number of photons in Hilbert space
Delta = cavity frequency - drive frequency 
Omega = Rabi frequency of drive field 
kappa = decay rate of cavity 
U(N) = non-linearity function
	Typical non-linearities are 
	U = U/2(N)(N-1) (Kerr-type)
	U = g sqrt(N)	(pseudo-Jaynes-Cummings)

dt = time step size 
Nsteps = number of time steps computed for 
Ntraj = number of independent trajectories taken in ensemble 
"""

"""
Dynamical parameters.

Intial state: 
state(0) = ground state 

Quantities Caculated: 
<a-> = photon coherent field in steady state 
<N> = expected photon count in steady state 		 
"""

import numpy as np
import matplotlib as mpl

class Ket:
	"""This is a class for the ket of a photon state in the Hilber space."""
	def __init__(self,nmax):
		self.Nmax = nmax	
		#The size of the Hilbert space 
		
		self.components = np.zeros( self.Nmax, dtype=np.complex ) 	
		#This vector is the coefficients of the state in the photon number basis
		#That is, components[n] = < n | state(t) > where |n> is the n-photon state	
	
		self.components[0] = 1.0	
		#We initialize the state into the ground state

	def innerProduct(self, bra):
		"""This computes the inner product of the ket 'self' with the bra formed from 'bra'"""
		return np.inner(self, np.conjugate(bra.components))


	def norm(self):
		"""A special case of the inner product of self with itself."""
		return self.innerProduct(self)

class Operator:
	"""This class is for operators on the kets."""
	def __init__(self





 
