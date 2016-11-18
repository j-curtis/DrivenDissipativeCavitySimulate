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
	We will assume that our system is described by a tensor product of two finite dimensional Hilbert Spaces
	Thus, we will accept upon instantiation the two dimensionalities of the factor spaces
	We then construct the tensor product space with the appropriate dimensionality
	We also include a method that converts from the human-readable indexing of (n1,n2) 
	for tensor product spaces to the internal flattened indexing 
	"""
	def __init__(self,d_phot,d_spin):
		"""
		This generates a class that represents a Hilbert space with the given number of states
		For the purposes of computation, the states will be enumerated by a single index
		The structure is a tensor product structure of |n_phot> X |n_spin> 
		"""
		self.phot_dimension = d_phot
		self.spin_dimension = d_spin
		self.dimension = self.phot_dimension*self.spin_dimension
		
	def tensorIndex(self,index_phot,index_spin):
		"""
		This maps an index of the form (spin_index,phot_index)
		onto the flattened index enumerating the whole Hilbert space via the map
		index = (index_spin)*(d_phot) + index_phot
		This has the structure of 
		index_spin = 0 => (index = 0,1,....d_phot-1)
		index_spin = 1 => (index = d_phot,d_phot+1,.....2d_phot-1 )
		
		As an example, for d_phot = 3, d_spin = 2 then the states are 
		00 = 0 ,10 = 1 ,20 = 2
		01 = 4 ,11 = 5 ,21 = 6
		where ij = |i; phot> |j; spin>
		"""
		return self.phot_dimension*index_spin + index_phot

	def vacuumIndex(self):
		"""
		Returns the index of the vacuum state
		We choose the vacuum state to be the 0 index
		"""	
		return 0	

class State:
	"""
	This is a class representing a datastructure for the wavefunction.
	It is a complex-valued vector of size HilbertSpace.dimension
	If our system is a tensor product structure, it will be assumed that the index has been flattened into a single set 
	"""

	def __init__(self, hilb_space):
		"""
		Initializes a vacuum state in the Hilbert Space hilb_space
		"""
		self.HSpace = hilb_space
		self.ket = np.zeros(shape=self.HSpace.dimension,dtype=np.complex)	
		self.ket[self.HSpace.vacuumIndex()] = 1.0	#We initialize the state into the ground state of both the cavity and atom 

	def __repr__(self):	
		"""
		We represent the state by simply returning the representation of the ket (and use the numpy __repr__ scaffolding)
		"""

		return self.ket.__repr__()

	def innerProduct(self, other):
		"""
		Computes the inner product 
		<other | self>
		"""
		return np.vdot(other.ket,self.ket)

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


class Operator:
	"""
	This is a class that implements operators on states.
	It will mostly serve as a unifying namespace for specific instances of operators.
	They will be represented by matrices
	"""
	def __init__(self,hilb_space):
		"""
		Initializes an operator on the Hilbert space given at initialization
		Initializes to identity operator
		"""
		self.HSpace = hilb_space
		
		#We initialize it to the identity operator
		self.matrix = np.eye(self.HSpace.dimension,dtype=np.complex)
		#We represent the matrix as (in the quantum matrix element notation)
		#self.matrix[i,j]= <i| Operator | j>
		#Where whatever indexing scheme we use for the state in the HIlbert Space is also employed here 	
 
	def __repr__(self):
		"""
		Returns the representation of an operator as the representation of its internal matrix
		"""
		return self.matrix.__repr__()


	def onKet(self, state):
		"""
		This implements the action of the operator (self) on the ket (state)
		Returns the new |ket'> = operator|ket>
		"""	
		new_state = State(self.HSpace)
		new_state.ket =  np.inner(self.matrix,state.ket)
		return new_state
	
	def __mul__(self,other):
		"""
		Computes the matrix product 
		self*other
		returns a new operator as this product
		"""
		Prod = Operator(self.HSpace)
		Prod.matrix = np.dot(self.matrix,other.matrix)
		return Prod

	def __add__(self,other):
		"""
		Computes the sum of two matrices by their component-wise sum
		"""
		Sum = Operator(self.HSpace)
		Sum.matrix = self.matrix + other.matrix
		return Sum

	def scalarMul(self,num):
		"""
		This implements the scalar multiplication of the operator (self) by the constant (num)
		"""
		Prod = Operator(self.HSpace)
		Prod.matrix = self.matrix*num
		return Prod

	def setPhotonLower(self):
		"""
		This sets the current operator to the lowering operator for photons, tensor product with the identity for spins 
		We use the formula 
		<m; phot| PhotonLowering | n; phot> = sqrt(n) delta(n-1 == m)
		"""
		self.matrix *= 0.0	#first we zero out all entries
		for i in np.arange(1,self.HSpace.phot_dimension):
			#i goes from 1 to photon_dimension-1
			for j in np.arange(self.HSpace.spin_dimension):
				#j ranges over all the spin states 
				#these are the elements <i-1; phot|<j; atom| a | j; atom>|i; phot> == sqrt(i)
				self.matrix[self.HSpace.tensorIndex(i-1,j),self.HSpace.tensorIndex(i,j)] = np.sqrt(i) 	

	def setPhotonRaise(self):
		"""
		Sets the current operator to the photon raising operator tensorproduct with the spin identity
		We use the formula for the non-zero matrix elements of 
		<m; phot| PhotonRaising|n; phot> = sqrt(m) delta(n+1==m)				
		"""
		self.matrix *= 0.0	#first zero out all matrix elements
		for i in np.arange(self.HSpace.phot_dimension-1):
			#i ranges from 0 to phot_dimension - 2 (the last state is destroyed in the truncated scheme)
			for j in np.arange(self.HSpace.spin_dimension):
				self.matrix[self.HSpace.tensorIndex(i+1,j),self.HSpace.tensorIndex(i,j)] = np.sqrt(i+1)

	

def main():
	num_spin = 2	#number of spin states we consider 
	num_phot = 3	#number of photon states we consider 

	HS = HilbertSpace(num_phot,num_spin)	
	
	state1 = State(HS)
	state2 = State(HS)
	state2.ket[0] = 0.0
	state2.ket[HS.tensorIndex(0,1)]	= 1.0

	print state1
	print state2

	Am = Operator(HS)
	Ap = Operator(HS)

	Am.setPhotonLower()
	Ap.setPhotonRaise()
	
	print (Ap*Ap).onKet(state1)
	print (Ap*Ap).onKet(state2)
	print (Ap*Am*Ap).onKet(state1)
if __name__=="__main__":
	main()






		
