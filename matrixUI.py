import streamlit as st
import numpy as np
import pandas as pd

from sympy.physics.quantum.dagger import Dagger
from sympy import *

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
    F,G = createFullColumnRank(A),createFullRowRank(A)
    # Calculates the Moore Penrose inverse of A
    X_M = Dagger(G) * (G * Dagger(G)).inv() * (Dagger(F)*F).inv()*Dagger(F)
    return X_M


st.set_page_config(
    page_title="Generalized Inverse Calculator",
    page_icon="📏",
    layout="wide",
    initial_sidebar_state="expanded"
)

#mono = st.text_input("Input a monomial")
#A = sympy.Matrix([[1,3,1,4],[2,7,3,9],[1,5,3,1],[1,2,0,8]])

#result = sympy.parsing.sympy_parser.parse_expr(A)

st.header("3x3 Generalized Inverse Solver")
st.text("Type in desired column elements and press Generate and obtain desired matrices")

#st.latex("A="+sympy.latex(A))
col1, col2, col3 = st.columns(3)

with col1:
   A_11 = col1.text_input("A_11")
   A_21 = col1.text_input("A_21")
   A_31 = col1.text_input("A_31")

with col2:
   A_12 = col2.text_input("A_12")
   A_22 = col2.text_input("A_22")
   A_32 = col2.text_input("A_32")

with col3:
    A_13 = col3.text_input("A_13")
    A_23 = col3.text_input("A_23")
    A_33 = col3.text_input("A_33")
  
generate = st.button("Generate")

if generate:
    A = Matrix([[int(A_11),int(A_12),int(A_13)],
                     [int(A_21),int(A_22),int(A_23)],
                     [int(A_31),int(A_32),int(A_33)]])
    st.header("Desired Matrix")
    st.latex("A="+latex(A))

    F = createFullColumnRank(A)

    st.header("Full Column Rank(F)")
    st.latex("F="+latex(F))

    G = createFullRowRank(A)

    st.header("Full Row Rank(G)")
    st.latex("G="+latex(G))

    st.header("Moore Penrose Inverse")
    X_M = MoorePenrose(A)
    st.latex("A^{\dagger}="+latex(X_M))

    st.header("Drazin Inverse")
    X_D = Drazin(A)
    st.latex("A^D="+latex(X_D))