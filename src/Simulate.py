"""
Jonathan Curtis
11/05/2016

See the theory page for an explanation/overview of the program.
"""

import numpy as np
import matplotlib as mpl
import random as rand

class HilbertSpace:
	"""
	This class contains the basis information about the Hilber space, like the dimensionality.
	"""
	def __init__(self,nspin,nphot):
		"""
		Initializes a new HilbertSpace class with the parameter for the photon space given.
		The Hilbert space has the structure of a tensor product between two Hilbert spaces of size 
		nspin and nphot
		Thus, the Hilbert Space has a structure of the form C^{nspin} X C^{nphot}
		"""
		self.Nspin = nspin
		self.Nphot = nphot
		self.shape = (self.Nspin,self.Nphot)
		self.dimension = self.Nspin*self.Nphot

	def indexRange(self):
		"""
		Returns a list of all the possible indices in the basis of this Hilbert space
		Equivalent to np.arange but generalized to the two-component index set
		"""
		return np.arange(

class State:
	"""
	This is a class that details the data structure used to store the information about the cavity/atom state.
	It is a tensor product of two complex vectors, one with 2 components (the atom) and one with Nmax components (cavity). 
	Thus, it is an element of the 2*Nmax dimensional complex Hilbert space.
	"""

	def __init__(self, hilbertspace):
		"""
		Initializes a vacuum state in the Hilbert Space hilbertspace
		"""
		self.HSpace = hilbertspace
		self.ket = np.zeros(shape=self.HSpace.shape,dtype=np.complex)	#We structure this as a tensor product
		"""
		self.ket[0,i] = <i; cavity|<0; atom| state>
		self.ket[1,i] = <i; cavity|<1; atom| state>
		That is, the first index is the projection of the state onto the cavity ground state and the second index is the cavity state.
		"""
		self.ket[0,0] = 1.0	#We initialize the state into the ground state of both the cavity and atom 

	def innerProduct(self, state):
		"""
		Computes the inner product 
		<state | self>
		"""
		bra = np.conj(state.ket)

		return np.tensordot(bra,self.ket,axes=2)

	def getNorm(self):
		"""
		This computes the norm of a state as 
		<state|state>
		"""
		return self.innerProduct(self)

	def normalize(self):
		"""
		This normalizes the state, if it is not already normalized.
		If the state is 0 it returns 0
		Otherwise it just divides by sqrt(innerProduct(self,self))
		It will set the internal ket to the normed ket and return the value of the norm before normalizing.
		"""
		norm = self.getNorm()
		if not norm == 0.0:
			self.ket = self.ket/np.sqrt(norm)	
		return norm

	def vacuum(self):
		"""
		This sets the current state to the vacuum |0>|0>
		It also returns the amplitude <0|<0|state>
		"""
		amplitude = self.ket[0,0]
		self.ket = self.ket*0.0
		self.ket[0,0] = 1.0
		return amplitude

class Operator:
	"""
	This is a class that implements operators on states.
	It will mostly serve as a unifying namespace for specific instances of operators.
	They will be represented by matrices
	"""
	def __init__(self,hilbertspace):
		"""
		Initializes an operator on the Hilbert space with photon Hilbert space of size Nmax
		Initializes to identity operator
		"""
		self.HSpace = hilbertspace

		self.matrix = np.eye(shape=(self.HSpace.shape,self.HSpace.shape),dtype=np.complex)
		#We represent the matrix as (in the quantum matrix element notation)
		#self.represent[i,j,k,l]= <j; cavity|<i; atom| Operator |k; atom> |l; cavity>	
		#We initialize it to the identity operator 

	def onKet(self, state):
		"""
		This implements the action of the operator (self) on the ket (state)
		"""	
		state.ket = np.tensordot(self.matrix,state.ket,axes=2)
		 
	def __mul__(self,other):
		"""
		This implements the matrix multiplication of 
		self.matrix * other.matrix
		Note that the order matters and the self is the left operator in this case.
		"""
		prod = operator(self.Nmax)
		prod.matrix = np.tensordot(self.matrix,other.matrix,axes=([2,0],[3,1]))		
		return prod	

def main():
	HS = HilbertSpace(nspin = 2,nphot = 3)	#Create a Hilbert space with 2 spin levels and 3 photon levels
	vac = State(HS)
	print vac.ket
	
	state1 = State(HS)
	state1.ket[0,0] = 1.0
	state1.ket[1,0] = 1.0
	state1.ket[1,1] = 1.0
	state1.normalize()
	print state1.ket
	
	state2 = State(HS)
	state2.ket[0,0] = 1.0
	state2.ket[1,0] = -1.0
	state2.ket[1,1] = 1.0
	state2.normalize()
	print state2.ket
	
	print state1.innerProduct(state2)
	
	identity = Operator(HS)

	print identity.matrix
	
if __name__=="__main__":
	main()






		
