#                     Chem 210A  C.Z. Pereyra
#                         Jan  2020
#
#        Homework 2: integration help for this problem set.
#
from sympy import *
from sympy import exp, oo
from sympy.plotting import plot
from matplotlib.pyplot import figure

#==============================================================================
# time-independent symbolic variables
#==============================================================================
#
L = symbols(' L ',real=True, positive=True)
x, t= symbols(' x t ', real=True)
n, m  = symbols('n m', integer=True, positive=True)  # note this must be integer, not int
M = symbols(' M ',real=True, positive= True)
h_bar = symbols(' h_bar ',real=True,positive=true)
alpha = symbols(' alpha ',real=True,positive=true)

#==============================================================================
# time-independent wave functions
#==============================================================================
#
def psi(n,L,x):
    result = sqrt(2/L)*sin(n*pi*x/L)
    return result

def psiPrime(n,L,x):
    result = diff(psi(1,L,x),x,1)
    return result

def psiL_Half(n,L,x):
    result = sqrt(4/L)*sin(n*pi*x/L/2)
    return result

def psiL_Full(n,L,x):
    result = sqrt(2/L)*sin(n*pi*x/L)
    return result

#==============================================================================
# integration operation
#==============================================================================
# problem 1. d.)
#i1_1 = integrate(psi(1,L,x) * x * diff(psi(1,L,x),x,1), (x,0,L))
i1_1 = integrate(psi(1,L,x) * x * psi(1,L,x),(x,0,L))
i2_2 = integrate(psi(2,L,x) * x * psi(2,L,x),(x,0,L))
i1_2 = integrate(psi(1,L,x) * x * psi(2,L,x),(x,0,L))

# problem 1. e.)
px_avg1 = integrate(psi(1,L,x) * h_bar / 1.0j * psiPrime(1,L,x), (x,0,L))
px_avg2 = integrate(psi(2,L,x) * h_bar / 1.0j * psiPrime(2,L,x), (x,0,L))
px_avg3 = integrate(psi(2,L,x) * psiPrime(1,L,x), (x,0,L))
px_avg4 = integrate(psi(1,L,x) * h_bar / 1.0j * psiPrime(2,L,x), (x,0,L))

#==============================================================================
# show or print result
#==============================================================================
# problem 1. d.)
'''
print("integration of n=1, n=1")
pprint(i1_1)
print("integration of n=2, n=2")
pprint(i2_2)
print("integration of n=1, n=2")
pprint(i1_2)
'''

#problem 1. e.)
'''
print("integration of phi(n=1), dphi(n=1)")
pprint(px_avg1)
print("integration of phi(n=2), dphi(n=2)")
pprint(px_avg2)
print("integration of phi(n=2), dphi(n=1)")
pprint(px_avg3)
print("integration of phi(n=1), dphi(n=2)")
pprint(px_avg4)
'''

#problem 2. a.)
c1 = integrate(psiL_Half(2,L,x) * psiL_Full(2,L,x), (x,0,L/2))
'''
print("integration of C1 -> phi1 * psi1")
pprint(c1)
'''
#problem 3. a.)
#
h = []
h.append(1)
h.append(2*x)

def H(n,x):
	##
	# can only use this in sequencial order 1,2,3,4 ..., etc
	if (n==0): return h[n]
	elif (n==1): return h[n]
	elif (n>1):
		result = 2*x*h[n-1] - 2*(n-1)*h[n-2]
		h.append(result)

	return result

for i in range(0, 4):
	print("h[{}] = {}".format(i, H(i,x)))

#print("h[1]")
#print("====")
#pprint(H(1,x))

#problem 3. b.)
#
def psi1_HMO(n,L,x):
#    result = exp(-alpha**2 * x**2 / 2) * x
    result = exp(-alpha**2 * x**2 / 2) * H(1,x)
    return result

def psi2_HMO(n,L,x):
#    result = exp(-alpha**2 * x**2 / 2) * (2 * (alpha * x)**2 - 1)
    result = exp(-alpha**2 * x**2 / 2) * H(2,x)
    return result

print("H1")
pprint(H(1,x))

print("H2")
pprint(H(2,x))

print("integration of <psi1 * psi2>")
orthogonal = integrate(psi1_HMO(1,L,x) * psi2_HMO(2,L,x), (x,-oo,oo))
pprint(orthogonal)

alpha = 1
#figure(num=None, figsize=(7, 5), dpi=70, facecolor='w', edgecolor='k')
#p1 = plot( psi1_HMO(1,L,x)*psi2_HMO(2,L,x), psi1_HMO(1,L,x), psi2_HMO(2,L,x), show=False, xlabel= 'x', ylabel= '$\psi$(x)', legend = True)
#p1 = plot( psi1_HMO(1,L,x)*psi2_HMO(2,L,x), show=False, xlabel= 'x', ylabel= '$\psi$(x)', legend = True)
p1 = plot( psi1_HMO(1,L,x), psi2_HMO(2,L,x), show=False, xlabel= 'x', ylabel= '$\psi$(x)', legend = True)

'''
p1[0].line_color='black'
p1[0].label = '$<\psi_1 \cdot \psi_2> $'
p1[0].style='solid'
'''

p1[0].line_color='blue'
p1[0].label = '$<\psi_1>$'
p1[0].style='wireframe'

p1[1].line_color='red'
p1[1].style='wireframe'
p1[1].label = '$<\psi_2>$'

p1.show()


