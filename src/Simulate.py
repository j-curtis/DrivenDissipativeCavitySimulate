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
	We will assume that there is a unique ground state and will assign it the index i=0
	"""
	def __init__(self,dim):
		"""
		This generates a class that represents a Hilbert space with the given number of states
		For the purposes of computation, the states will be enumerated by a single index
		"""
		self.dimension = dim	

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
		self.ket[0] = 1.0	#We initialize the state into the ground state of both the cavity and atom 

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
		prod = Operator(self.HSpace)
		prod.matrix = np.dot(self.matrix,other.matrix)
		return prod

def main():
	num_spin = 2	#number of spin states we consider 
	num_phot = 3	#number of photon states we consider 
	dimension = num_spin * num_phot	#the dimension of the Hilbert space consisting of the tensor product 

	HS = HilbertSpace(dimension)	
	
	state1 = State(HS)
	state2 = State(HS)

	state1.ket[1] = 1.0
	state2.ket[1] = -1.0

	state1.normalize()
	state2.normalize()
	
	print state1
	print state2

	op1 = Operator(HS)
	op2 = Operator(HS)

	op1.matrix *= 0.0
	op2.matrix *= 0.0
	
	op1.matrix[0,1] = 1.0
	op1.matrix[1,0] = 1.0

	op2.matrix[1,1] = 1.0
	op2.matrix[0,0] = -1.0
	
	print op1
	print 
	print op2
	print 
	print op1*op2
	

if __name__=="__main__":
	main()






		
