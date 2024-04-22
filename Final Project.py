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
############### Finite Difference Method: 1D Harmonic ###############
#####################################################################

def get1DHarmonicCoefficientMatrix(n):
    A = np.zeros((n + 1, n + 1))
    A[0, 0] = 1
    for i in range(1, n): # Note: range() includes the first argument but not the second argument
        A[i, i - 1] = 1
        A[i, i] = -2
        A[i, i + 1] = 1
    A[n, n] = 1
    return A

def get1DHarmonicVector(n, h):
    b = np.zeros(n + 1)
    b[0] = 0
    b[1:n] = -9.8 * h**2 # Note: ** in Python represents ^, and n is not included in the range
    b[n] = 50
    return b

def solve1DHarmonicSystem(A, b):
    return np.linalg.solve(A, b)

#####################################################################
############## Finite Difference Method: 1D Biharmonic ##############
#####################################################################

def get1DBiharmonicCoefficientMatrix(n):
    A = np.zeros((n + 1, n + 1))
    A[0, 0] = 1
    A[1, 0] = -4
    A[1, 1] = 6
    A[1, 2] = -4
    A[1, 3] = 1
    for i in range(2, n - 1): # Note: range() includes the first argument but not the second argument
        A[i, i - 2] = 1
        A[i, i - 1] = -4
        A[i, i] = 6
        A[i, i + 1] = -4
        A[i, i + 2] = 1
    A[n - 1, n - 3] = 1
    A[n - 1, n - 2] = -4
    A[n - 1, n - 1] = 6
    A[n - 1, n] = -4
    A[n, n] = 1
    return A

# Use f(x) = cos() this time, where f''''(x) = cos().
def get1DBiharmonicVector(n, h):
    b = np.zeros(n + 1)
    b[0] = np.cos(0)
    for i in range(1, n):
        b[1:n] = np.cos(i * h) # Note: ** in Python represents ^, and n is not included in the range
    b[n] = np.cos(n)
    return b

def solve1DBiharmonicSystem(A, b):
    return np.linalg.solve(A, b)

#####################################################################
############### Finite Difference Method: 2D Harmonic ###############
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

np.set_printoptions(linewidth = 100) # Set threshold = np.inf to see the entire matrix

"""
# Test the 1D harmonic implementation.
n = 100
h = 10**-4

A = get1DHarmonicCoefficientMatrix(n)
print("A:", A)
b = get1DHarmonicVector(n, h)
print("b:", b)
print("solved:", solve1DHarmonicSystem(A, b))
"""

# Test the 1D biharmonic implementation.
n = 5
h = 1

A = get1DBiharmonicCoefficientMatrix(n)
print("A:", A)
b = get1DBiharmonicVector(n, h)
print("b:", b)

print()
for a in range(0, 5):
    print(np.round(np.cos(a), 8))
print("solved:", solve1DBiharmonicSystem(A, b))