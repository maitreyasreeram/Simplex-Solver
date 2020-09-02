"""
code to solve LP formulation provided by user input. Uses Tableau Method and Two Phase method to solve the simplex
"""

import math
import numpy as np
from numpy.linalg import *
import sys

# Function 1: Use this function to take inputs from user as standard form (do not input the artificial variables as yet- the code will detect that)
def take_initialize_inputs():
    n1 = int(input("Enter number of decision variables: "))
    n2 = int(input("Enter number of slack variables: "))
    n = n1 + n2
    m = int(input("Enter number of constraints: "))
    A = np.zeros((m,n))
    b = np.zeros(m)
    C = np.zeros(n)
    print("Enter the cost matrix (including slack): ")
    for i in range(n):
        C[i] = float(input("C(x{})".format(i)))
    print("Entered cost matrix is: ", C)
    print("Enter the matrix of coefficients: ")
    for i in range(m):
        for j in range(n):
            
            A[i][j] = float(input("A({},{})".format(i,j)))
    print("Entered coefficient matrix is: ", A)
    print("Enter the RHS Matrix: ")
    for i in range(m):
        b[i] = float(input("b({})".format(i)))
    print("Entered RHS Martix: ", b)
    return A, b, C

# Use this function if artificial variables exit - and update coefficient matrix and cost matrix
def input_artificial_variables(A, C):

    a = int(input("Initial Identity basis does not exist. Enter number of artifical variables: "))
    m = len(A)
    n = len(A[0])
    n = n + a
    A_new = np.zeros((m,n))
    C_new = np.zeros(n)
    print("Enter the new coefficient matrix (with slack and artificial): ")
    for i in range(m):
        for j in range(n): 
            A_new[i][j] = float(input("A({},{})".format(i,j)))
    print("Entered coefficient matrix is: ", A)
    print("Enter the new cost matrix: ")
    for i in range(n): 
        C_new[i] = float(input("C(x{})".format(i)))
    print("Entered cost matrix is: ", C_new)
    return A_new, C_new

# checking for identity matrix in the initial coefficient matrix

def check_identity_basis(A):
    identity = False
    ind = np.arange(np.shape(A)[1])
    I = np.identity(len(A))
    loc = []
    idd = {}
    basis = []
    for i in range(np.shape(A)[1]):
        for j in range(len(I)):
            if (A[:,i] == I[:,j]).all(): # column by column comparison
                idd[j] = i
                loc.append(i)

    keys = list(idd.keys()) # maintaining list of matches
    vals = list(idd.values())
    if len(set(vals)) == len(I): # if matches == dimension of B
        print("Identity exists")
        sorted_index = sorted(keys)
        for i in sorted_index:
            basis.append(idd[i])
        print("Basis is: ", basis)
        identity = True
    else:
        print("Identity does not exist")
    non_basis = list(set(ind) - set(basis))
    return basis, non_basis, identity, idd

# Check initial feasibility of coefficient matrix based on rank and b

def check_initial_feasibility(A, b):
    feas1 = True
    feas2 = True
    
    aug_A = np.c_[A,b]
    #print(aug_A)
    #print(matrix_rank(A), matrix_rank(aug_A))
    if matrix_rank(A) < matrix_rank(aug_A):
        feas1 = False
        print("Infeasible due to rank")
    #print(feas1)
    for i in range(len(b)):
        if b[i] < 0: 
            print("Infeasible solution due to b")

            feas2 = False
            break
    #print(feas2)
    feas = feas1 or feas2
    return feas

# Redundancy check  - if redundant - delete a row at random and check if matrix is still redundant 

def redundancy_check(A):
    b = len(A)
    while(matrix_rank(A) < min(len(A), len(A[0]))) and b > 1:

        print("Redundant constraints")
        A = np.delete(A, b-1, 0)
        #print("A", A)
        b = b - 1
        #print("B", b)
    return A

# checking degeneracy based on whether a 0 exists in basis - if yes, True

def degeneracy_check(tab): 
    degen_ind = []
    degen = False
    for i,j in enumerate(tab[1:,-1]):
        if j == 0: 
            degen = True
            degen_ind.append(i)
    
    return degen, np.asarray(degen_ind)
        
# deconstructing the coefficient matrix into B and N

def deconstruct(A, C):
    A_t = A[:]
    m = len(A_t)
    B = np.empty((m,m))
    
    basis_index, non_basis_index, check, idd= check_identity_basis(A_t)
    for i,j in enumerate(basis_index):
        B[:,i] = A_t[:,j]


    N = np.delete(A_t, [basis_index], 1)


    return B, N

# constructing initial tableau 

def construct_tableau(A, b, C, Ib, In):
    m = len(A)
    n = len(C)
    tableau = np.zeros((m+1, n+1))
    #print(tableau)
    B,N = deconstruct(A, C)
    tableau[1:,:n] = A
    tableau[0,:n] = -C

    tableau[1:,n] = np.transpose(b)
    
    print(tableau)
    return tableau
        
# Initialize tableau to make Reduced costs of basis variables 0
def tableau_initial(tab, Ib):

    for i, j in enumerate(Ib):
        if tab[0][j] !=0: 
            #print(i,j)
            tab[0][:] = -tab[0][j]*tab[i+1][:] + tab[0][:]
        else: 
            pass

    return tab

# Ratio test for exit basis

def ratio_test(col, end):
    ratio = {}
    degen = {}
    if (col<0).all(): 
        print("UNBOUNDED LP")
        sys.exit()
        pass
    else: 
        for i,j in enumerate(col):
            if j >0:
                ratio[end[i]/col[i]] = i+1
            elif j == 0: 
                degen[j] = i+1
            else:
                pass
    return ratio
        
# choosing which indexes will enter and exit - Degeneracy check also done- anti-cycling rule (Blands Rule)

def enter_exit_basis(tab, In, Ib):
    max_RC = max(tab[0][In])
    if max_RC <=0: 
        print("TABLEAU IS OPTIMAL", tab)
        sys.exit()
    else: 
        enter_basis = np.where(tab[0][:] == max_RC)[0][0]

        print("Entering basis is ", enter_basis)
        degen, degen_ind = degeneracy_check(tab)
        if degen: 
            print("DEGENERATE LP - USING BLANDS RULE ")
            if len(degen_ind) == 0 and tab[1:,enter_basis][degen_ind] > 0:
                exit_basis_value = min(np.array(Ib)[degen_ind][np.array(Ib)[degen_ind]>0])
                exit_basis = Ib.index(exit_basis_value)+1
            else: 
                rat = ratio_test(tab[1:, enter_basis], tab[1:,-1])
                exit_basis = rat[min(list(rat.keys()))]
                print("Exiting basis is ", exit_basis)
                
        else: 

            rat = ratio_test(tab[1:, enter_basis], tab[1:,-1])
            exit_basis = rat[min(list(rat.keys()))]
            print("Exiting basis is ", exit_basis)
    return exit_basis, enter_basis

# Pivoting rules based on enter-exit rules

def pivot(tab, Ib, In, ext, ent):
    print("Ib is ", Ib)
    print("In is ", In)
    print(tab)
    tab[ext][:] = tab[ext][:]/tab[ext][ent]
    print("After First step \n", tab)
    for i in (set(np.arange(len(tab))) - set([ext])):
        print("Currently on ",i)
        print("R{} = {}R{} + R{}".format(i,-tab[i][ent],ext,i))
        tab[i][:] = -tab[i][ent]*tab[ext][:] + tab[i][:]
        print(tab)
    print(tab)
    xx = Ib[ext-1]
    del Ib[ext-1]
    Ib.insert(ext-1,ent)
    In.remove(ent)
    In.append(xx)
    return tab, Ib, In

# generic function to solve tableau with pivoting and ratio test 

def solve_tableau(tab, Ib, In):

    print("Solving tableau")
    print(tab)
    ext, ent = enter_exit_basis(tab, In, Ib)
    tab, Ib, In = pivot(tab, Ib, In, ext, ent)
    return tab, Ib, In
    
# function for phase 1: different cost function, artificial variables included

def Phase1(A, b, Ib, In,n, n_a):
    C_phase1 = np.zeros(n_a)
    for i in range(n,n_a):
        C_phase1[i] = 1
    print("C for phase 1 is: ", C_phase1)
    tab=construct_tableau(A, b, C_phase1, Ib, In)
    tab = tableau_initial(tab, Ib)
    # convert basis reduced cost to zero first
    while not (tab[0][-1] == 0 and (tab[0][:]<=0).all()):  
        tab, Ib, In = solve_tableau(tab, Ib, In)
    print("*********************PHASE 1 COMPLETE********************* ")
    return tab, Ib, In

# Function for Phase 2
def Phase2(tab, Ib, In, C):
    In2 = In[:]
    tab  = np.delete(tab, In[-1], 1)
    del In2[-1]
    C_phase2 = -C[:-1]
    tab[0,:-1] = C_phase2
    tab = tableau_initial(tab, Ib)
    while not ((tab[0][:-1]<=0).all()):  
        tab, Ib, In2 = solve_tableau(tab, Ib, In2)
    print("*********************PHASE 2 COMPLETE********************* ")
    print("Optimal Obj function value for minimization is: ", tab[0,-1])
    print("Optimal solution is: {} with basis index: {}".format(tab[1:,-1], Ib ))
    return tab, Ib, In2

# main function - workflow for solver

def main():
    A, b, C = take_initialize_inputs()
    m = len(A)
    n = len(A[0])
    Ib, In, check, idd = check_identity_basis(A)
    print(Ib, In)
    while not check:
        A, C = input_artificial_variables(A, C)
        print(A)
        n_a = len(A[0])
        Ib, In, check, idd = check_identity_basis(A)
    # feasibility check 
    feas = check_initial_feasibility(A,b)
    if not feas: 
        print("INFEASIBLE")
        
    a = n_a - n
    print("a is ",a)
    # begin phase 1 if artificial exist
    print("Entering PHASE 1")
    tab, Ib, In = Phase1(A, b, Ib, In, n, n_a)
    tab, Ib, In = Phase2(tab, Ib, In, C)
    print(tab, Ib, In)

main() # this function is run to start the solver
