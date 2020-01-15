#                     Chem 210A  C.Z. Pereyra
#                         Jan  2020
#
#        Homework 1, problem 6
#
from sympy import *

#==============================================================================
# time-independent symbolic variables
#==============================================================================
#
L = symbols(' L ',real=True, positive=True)
x, t= symbols(' x t ', real=True)
n, m  = symbols('n m', integer=True, positive=True)  # note this must be integer, not int
M = symbols(' M ',real=True, positive= True)
h_bar = symbols(' h_bar ',real=True,positive=true)

#==============================================================================
# time-independent wave functions
#==============================================================================
#
def psi(n,L,x):
    result = sqrt(2/L)*sin(n*pi*x/L)
    return result

#==============================================================================
# commutation operation
#==============================================================================
#
comm = (h_bar / 1.0j) * (x * diff(psi(n,L,x),x,1) - diff(x * psi(n,L,x),x,1))
#comm = Commutator(x, d) * psi(n,L,x)

#==============================================================================
# show or print result
#==============================================================================
#
pprint(comm)


