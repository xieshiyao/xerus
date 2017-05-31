import xerus as xe

MAX_NUM_PER_SITE = 64

def create_M():
	M = -1 * xe.Tensor.identity([MAX_NUM_PER_SITE, MAX_NUM_PER_SITE])
	for i in xrange(MAX_NUM_PER_SITE-1) :
		M[[i+1, i]] = 1.0
	return M

def create_L():
	L = xe.Tensor([MAX_NUM_PER_SITE, MAX_NUM_PER_SITE])
	for i in xrange(MAX_NUM_PER_SITE) :
		L[[i,i]] = i / (i+5.0)
	return L

def create_S():
	S = xe.Tensor([MAX_NUM_PER_SITE, MAX_NUM_PER_SITE])
	
	# set diagonal
	for i in xrange(MAX_NUM_PER_SITE) :
		S[[i,i]] = -i
	
	# set offdiagonal
	for i in xrange(MAX_NUM_PER_SITE-1) :
		S[[i,i+1]] = i+1
	
	return 0.07*S


def create_operator(degree):
	i,j,k,l = xe.indices(4)
	
	# create matrices
	M = create_M()
	S = create_S()
	L = create_L()
	Sstar = 0.7*M + S;
	I = xe.Tensor.identity([MAX_NUM_PER_SITE, MAX_NUM_PER_SITE])
	
	# create empty TTOperator
	A = xe.TTOperator(2*degree)
	
	# create first component
	comp = xe.Tensor()
	comp(i, j, k, l) << \
		Sstar(j, k) * xe.Tensor.dirac([1, 3], 0)(i, l) \
		+   L(j, k) * xe.Tensor.dirac([1, 3], 1)(i, l) \
		+   I(j, k) * xe.Tensor.dirac([1, 3], 2)(i, l)
	
	A.set_component(0, comp)
	
	# create middle components
	comp(i, j, k, l) << \
		  I(j, k) * xe.Tensor.dirac([3, 3], [0, 0])(i, l) \
		+ M(j, k) * xe.Tensor.dirac([3, 3], [1, 0])(i, l) \
		+ S(j, k) * xe.Tensor.dirac([3, 3], [2, 0])(i, l) \
		+ L(j, k) * xe.Tensor.dirac([3, 3], [2, 1])(i, l) \
		+ I(j, k) * xe.Tensor.dirac([3, 3], [2, 2])(i, l)
	
	for c in xrange(1, degree-1) :
		A.set_component(c, comp)
	
	# create last component
	comp(i, j, k, l) << \
		  I(j, k) * xe.Tensor.dirac([3, 1], 0)(i, l) \
		+ M(j, k) * xe.Tensor.dirac([3, 1], 1)(i, l) \
		+ S(j, k) * xe.Tensor.dirac([3, 1], 2)(i, l)
	
	A.set_component(degree-1, comp)
	
	return A


def one_norm(x):
	i = xe.Index()
	return float(x(i&0) * xe.TTTensor.ones(x.dimensions)(i&0))

def implicit_euler(A, x, stepSize, n):
	op = xe.TTOperator.identity(A.dimensions) - stepSize*A
	
	j,k = xe.indices(2)
	ourALS = xe.ALS_SPD
	ourALS.convergenceEpsilon = 1e-4
	ourALS.numHalfSweeps = 100
	
	results = [x]
	nextX = xe.TTTensor(x)
	
	for i in xrange(n) :
		ourALS(op, nextX, x)
		
		# normalize
		norm = one_norm(nextX)
		nextX /= norm
		
		print("done itr", i, \
			"residual:", xe.frob_norm(op(j/2,k/2)*nextX(k&0) - x(j&0)), \
			"one-norm:", norm)
		
		x = xe.TTTensor(nextX) # ensure it is a copy
		results.append(x)
	
	return results


def get_mean_concentration(x, i):
	k,l = xe.indices(2)
	result = xe.TensorNetwork(x)
	weights = xe.Tensor.from_function([MAX_NUM_PER_SITE], lambda idx: idx[0])
	ones = xe.Tensor.ones([MAX_NUM_PER_SITE])
	
	for j in xrange(x.degree()) :
		if j == i :
			result(l&0) << result(k, l&1) * weights(k)
		else :
			result(l&0) << result(k, l&1) * ones(k)
	
	# at this point the degree of 'result' is 0, so there is only one entry
	return result[[]]

def print_mean_concentrations_to_file(results):
	f = open("mean.dat", 'w')
	for res in results :
		for k in xrange(res.degree()) :
			f.write(str(get_mean_concentration(res, k))+' ')
		f.write('\n')
	f.close()


numProteins = 10
numSteps = 200
stepSize = 1.0
rankX = 3

A = create_operator(numProteins)

start = xe.TTTensor.dirac([MAX_NUM_PER_SITE]*numProteins, 0)
for i in xrange(numProteins) :
	start.set_component(i, start.get_component(i).dense_copy())

start += 1e-14 * xe.TTTensor.random(start.dimensions, [rankX-1]*(start.degree()-1))

results = implicit_euler(A, start, stepSize, numSteps)

print_mean_concentrations_to_file(results)
