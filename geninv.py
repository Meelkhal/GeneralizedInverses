# this is a python script which calculates the Moore-Penrose and the Drazin Inverse of a matrix

# first import the necessary libraries
from sympy.physics.quantum.dagger import Dagger
from sympy import *

# First we create the function which generate F and G the full column rank and row rank matrices respectively.

def createFullColumnRank(A):
    # Takes a sympymatrix A and calculates the full column rank matrix of A 
    F = 1*A 
    # We define F = 1*A, if we had tried doing F = A, then any changes made to F would be applied to A, which we don't want as 
    # it may be useful to retain the original matrix A
    rows,columns = shape(F)
    RREF,pivot_indexs = F.rref()
    for col in range(columns):
        if col not in pivot_indexs:
            F.col_del(col)
    return F

def createFullRowRank(A):
    # Takes a sympymatrix A and creates a full row rank matrix based on the RREF of A
    G = (1*A).rref()[0]
    rows,columns = shape(G)
    for r in range(rows):
        if G.row(r) == Matrix([[0]*columns]):
            G.row_del(r)
    return G

def Drazin(A):
    # Calculates the drazin inverse of A
    F,G = createFullColumnRank(A),createFullRowRank(A)
    X_D = F*(G*F*G*F).inv()*G
    return X_D

def MoorePenrose(A):
    # Calculates the Moore Penrose inverse of A
    X_M = Dagger(G) * (G * Dagger(G)).inv() * (Dagger(F)*F).inv()*Dagger(F)
    return X_M