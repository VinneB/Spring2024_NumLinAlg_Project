# Note: Indexing starts at 0 in Python.

# Import the required libraries.
import numpy as np
import sys


def main():
    if sys.argv[1].upper() == "TEMP":
        temp()



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

###################################################################################################
############## Finite Difference Method: 1D Biharmonic (2nd Convergence)(Neumann BC) ##############
###################################################################################################

def bi_dim1_mat(n, h):
    A = np.diag(np.full(n, -2, dtype=float)) + \
    np.diag(np.ones(n-1, dtype=float)*-4, -1) + \
    np.diag(np.ones(n-1, dtype=float)*-4, 1) + \
    np.diag(np.ones(n-2, dtype=float), -2) + \
    np.diag(np.ones(n-2, dtype=float), 2)
    # dtype unncessary

    # Neumann Bounds
    A[0] = np.zeros(n)
    A[0][0] = float(h ** 4)
    A[n-1] = np.zeros(n)
    A[n-1][n-1] = float(h ** 4)

    # Non-Centered Diff for 1 and n-2
    A[1] = np.concatenate((np.array([3, -14, 26, -24, 11, -2]), np.zeros(n-6)))
    A[n-2] = np.flip(A[1])

    return A

# Use f(x) = cos() this time, where f''''(x) = cos().
def bi_dim1_vec(mesh, func, bc):
    b = np.array(list(map(func, mesh)))
    b[0] = bc[0]
    b[-1] = bc[1]
    return b

# Returns 
# If a stepsize (h) and number of spacial values are specified, stepsize takes precedence
def bi_dim1_solve(func, range, bc, n=-1, h=-1):
    if (n==-1 and h==-1):
        raise Exception("must specify n or h")
    
    range[0] = float(range[0])
    range[1] = float(range[1])
    #bc[0] = float(bc[0])
    #bc[1] = float(bc[1])

    if h != -1:
        mesh = np.arange(range[0], range[1]+h, h)
        n = int(((range[1] - range[0])/h)+1)
    else:
        h = (range[1] - range[0])/float(n)
        mesh = np.linspace(range[0], range[1], n)
    
    #print(mesh)

    # Generate b
    b = bi_dim1_vec(mesh, func, bc)

    # Generate difference matrix
    if (n == -1):
        n = (range[1] - range[0])/h
    elif (h == -1):
        h = (range[1] - range[0])/n
    D = bi_dim1_mat(n, h)
    return mesh, b, D, [n,h], np.linalg.solve(D, b * np.dot(np.ones(n), (h ** 4)))

    


"""     b = np.zeros(n + 1)
    b[0] = np.cos(0)
    for i in range(1, n):
        b[1:n] = np.cos(i * h) # Note: ** in Python represents ^, and n is not included in the range
    b[n] = np.cos(n)
    return b """

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



def temp():
    np.set_printoptions(precision=3)
    np.set_printoptions(linewidth = 100) # Set threshold = np.inf to see the entire matrix
    np.set_printoptions(suppress=True)
    mesh, b, D, dim, res = bi_dim1_solve(np.cos, [0., 3.], [np.cos(0), np.cos(3)], h = .25)
    #h, n = bi_dim1_solve(np.cos, [0, 3], [np.cos(0), np.cos(3)], h = .25)
    #print("h = {} ; n = {}".format(h, n))
    #h, n = bi_dim1_solve(np.cos, [0, 3], [np.cos(0), np.cos(3)], n = 12)
    #print("h = {} ; n = {}".format(h, n))
    #return
    print("Result")
    print(res)
    print("b")
    print(b)
    print("mesh")
    print(mesh)
    print("n: {}  ;  h: {}".format(dim[0], dim[1]))
    np.savetxt("1dmatrix.txt", D, delimiter=" ", fmt='%2.0f')

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

""" # Test the 1D biharmonic implementation.
n = 5
h = 1

A = get1DBiharmonicCoefficientMatrix(n)
print("A:", A)
b = get1DBiharmonicVector(n, h)
print("b:", b)

print()
for a in range(0, 5):
    print(np.round(np.cos(a), 8))
print("solved:", solve1DBiharmonicSystem(A, b)) """

if __name__ == "__main__":
    main()
