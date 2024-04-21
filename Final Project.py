# Note: Indexing starts at 0 in Python.

# Import the required libraries.
import numpy as np

# Set the boundary conditions, if any.
def setBoundaryConditions(x1, y1, x2, y2):
    A[x1, x1] = y1
    b[x1] = y1
    A[x2, x2] = y2
    b[x2] = y2

#####################################################################
################### Finite Difference Method: 1D ####################
#####################################################################

def get1DCoefficientMatrix(n):
    A = np.zeros((n + 1, n + 1))
    A[0, 0] = 1
    for i in range(1, n): # Note: range() includes the first argument but not the second argument
        A[i, i - 1] = 1
        A[i, i] = -2
        A[i, i + 1] = 1
    A[n, n] = 1
    return A

def get1DbVector(n, h):
    b = np.zeros(n + 1)
    b[0] = 0
    b[1:n] = -9.8 * h**2 # Note: ** in Python represents ^, and n is not included in the range
    b[n] = 50
    return b

def solve1DSystem(A, b):
    return np.linalg.solve(A, b)

#####################################################################
################### Finite Difference Method: 2D ####################
#####################################################################

n = 10
Nx = n - 1
h_top = 0

def get2DCoefficientMatrix():
    Ddiag  = -4 * np.eye(Nx - 1)
    Dupper = np.diag([1] * (Nx - 2), 1)
    Dlower = np.diag([1] * (Nx - 2), -1)
    D = Bdiag + Bupper + Blower
    Ds = [D] * (Nx - 1)
    A = lin.block_diag(*Ds)
    I = np.ones((Nx - 1) * (Nx - 2))
    Iupper = np.diag(I, Nx - 1)
    Ilower = np.diag(I, -Nx + 1)
    A += Iupper + Ilower
    return A

def get2DbVector():
    b = np.zeros((Nx - 1) ** 2)
    b[-Nx + 1:] = -h_top
    return b

def solve2DSystem(A, b):
    solved = np.linalg.solve(A, b)
    return solved


#####################################################################
########################## Test Functions ###########################
#####################################################################
n = 100
h = 10**-4

A = get1DCoefficientMatrix(n)
print(A)
b = get1DbVector(n, h)
print(b)
print(solve1DSystem(A, b))