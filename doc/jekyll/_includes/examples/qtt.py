import xerus as xe

# construct the stiffness matrix A using a fill function
def A_fill(idx):
	if idx[0] == idx[1] :
		return 2.0
	elif idx[1] == idx[0]+1 or idx[1]+1 == idx[0] :
		return -1.0
	else:
		return 0.0

A = xe.Tensor.from_function([512,512], A_fill)

# and dividing it by h^2 = multiplying it with N^2
A *= 512*512

# reinterpret the 512x512 tensor as a 2^18 tensor
# and create (Q)TT decomposition of it
A.reinterpret_dimensions([2,]*18)
ttA = xe.TTOperator(A)

# and verify its rank
print("ttA ranks:", ttA.ranks())

# the right hand side of the equation both as Tensor and in (Q)TT format
b = xe.Tensor.ones([2,]*9)
ttb = xe.TTTensor.ones(b.dimensions)

# construct a random initial guess of rank 3 for the ALS algorithm
ttx = xe.TTTensor.random([2,]*9, [3,]*8)

# and solve the system with the default ALS algorithm for symmetric positive operators
xe.ALS_SPD(ttA, ttx, ttb)

# to perform arithmetic operations we need to define some indices
i,j,k = xe.indices(3)

# calculate the residual of the just solved system to evaluate its accuracy
# here i^9 denotes a multiindex named i of dimension 9 (ie. spanning 9 indices of the respective tensors)
residual = xe.frob_norm( ttA(i^9,j^9)*ttx(j^9) - ttb(i^9) )
print("residual:", residual)

# as an comparison solve the system exactly using the Tensor / operator
x = xe.Tensor()
x(j^9) << b(i^9) / A(i^9, j^9)

# and calculate the Frobenius norm of the difference
print("error:", xe.frob_norm(x - xe.Tensor(ttx)))

