import xerus as xe

class InternalSolver :
	def __init__(self, A, x, b):
		self.A = A
		self.b = b
		self.x = x
		self.d = x.degree()
		self.solutionsNorm = b.frob_norm()
		
		self.leftAStack = [ xe.Tensor.ones([1,1,1]) ]
		self.leftBStack = [ xe.Tensor.ones([1,1]) ]
		self.rightAStack = [ xe.Tensor.ones([1,1,1]) ]
		self.rightBStack = [ xe.Tensor.ones([1,1]) ]
		
		self.maxIterations = 1000
	
	
	def push_left_stack(self, pos) :
		i1,i2,i3, j1,j2,j3, k1,k2 = xe.indices(8)
		Ai = self.A.get_component(pos)
		xi = self.x.get_component(pos)
		bi = self.b.get_component(pos)
		
		tmpA = xe.Tensor()
		tmpB = xe.Tensor()
		tmpA(i1, i2, i3) << self.leftAStack[-1](j1, j2, j3)\
				*xi(j1, k1, i1)*Ai(j2, k1, k2, i2)*xi(j3, k2, i3)
		self.leftAStack.append(tmpA)
		tmpB(i1, i2) << self.leftBStack[-1](j1, j2)\
				*xi(j1, k1, i1)*bi(j2, k1, i2)
		self.leftBStack.append(tmpB)
	
	def push_right_stack(self, pos) :
		i1,i2,i3, j1,j2,j3, k1,k2 = xe.indices(8)
		Ai = self.A.get_component(pos)
		xi = self.x.get_component(pos)
		bi = self.b.get_component(pos)
		
		tmpA = xe.Tensor()
		tmpB = xe.Tensor()
		tmpA(j1, j2, j3) << xi(j1, k1, i1)*Ai(j2, k1, k2, i2)*xi(j3, k2, i3) \
				* self.rightAStack[-1](i1, i2, i3)
		self.rightAStack.append(tmpA)
		tmpB(j1, j2) << xi(j1, k1, i1)*bi(j2, k1, i2) \
				* self.rightBStack[-1](i1, i2)
		self.rightBStack.append(tmpB)
	
	
	def calc_residual_norm(self) :
		i,j = xe.indices(2)
		return xe.frob_norm(self.A(i/2, j/2)*self.x(j&0) - self.b(i&0)) / self.solutionsNorm
	
	
	def solve(self) :
		# build right stack
		self.x.move_core(0, True)
		for pos in reversed(xrange(1, self.d)) :
			self.push_right_stack(pos)
		
		i1,i2,i3, j1,j2,j3, k1,k2 = xe.indices(8)
		residuals = [1000]*10
		
		for itr in xrange(self.maxIterations) :
			residuals.append(self.calc_residual_norm())
			if residuals[-1]/residuals[-10] > 0.99 :
				print("Done! Residual decreased from:", residuals[10], "to", residuals[-1], "in", len(residuals)-10, "sweeps")
				return
			
			print("Iteration:",itr, "Residual:", residuals[-1])
			
			# sweep left -> right
			for pos in xrange(self.d):
				op = xe.Tensor()
				rhs = xe.Tensor()
				
				Ai = self.A.get_component(pos)
				bi = self.b.get_component(pos)
				
				op(i1, i2, i3, j1, j2, j3) << self.leftAStack[-1](i1, k1, j1)*Ai(k1, i2, j2, k2)*self.rightAStack[-1](i3, k2, j3)
				rhs(i1, i2, i3) <<            self.leftBStack[-1](i1, k1) *   bi(k1, i2, k2) *   self.rightBStack[-1](i3, k2)
				
				tmp = xe.Tensor()
				tmp(i1&0) << rhs(j1&0) / op(j1/2, i1/2)
				self.x.set_component(pos, tmp)
				
				if pos+1 < self.d :
					self.x.move_core(pos+1, True)
					self.push_left_stack(pos)
					self.rightAStack.pop()
					self.rightBStack.pop()
			
			
			# right -> left, only move core and update stack
			self.x.move_core(0, True)
			for pos in reversed(xrange(1,self.d)) :
				self.push_right_stack(pos)
				self.leftAStack.pop()
				self.leftBStack.pop()


def simpleALS(A, x, b) :
	solver = InternalSolver(A, x, b)
	solver.solve()

if __name__ == "__main__":
	i,j,k = xe.indices(3)
	
	A = xe.TTOperator.random([4]*16, [2]*7)
	A(i/2,j/2) << A(i/2, k/2) * A(j/2, k/2)
	
	solution = xe.TTTensor.random([4]*8, [3]*7)
	b = xe.TTTensor()
	b(i&0) << A(i/2, j/2) * solution(j&0)
	
	x = xe.TTTensor.random([4]*8, [3]*7)
	simpleALS(A, x, b)
	
	print("Residual:", xe.frob_norm(A(i/2, j/2) * x(j&0) - b(i&0))/xe.frob_norm(b))
	print("Error:", xe.frob_norm(solution-x)/xe.frob_norm(x))
	
