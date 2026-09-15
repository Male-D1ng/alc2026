#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Eliminacion Gausianna
"""
import numpy as np

'''
def elim_gaussiana(A):
    cant_op = 0
    m=A.shape[0]
    n=A.shape[1]
    Ac = A.copy()
    
    if m!=n:
        print('Matriz no cuadrada')
        return
    
    ## desde aqui -- CODIGO A COMPLETAR





                
    ## hasta aqui, calculando L, U y la cantidad de operaciones sobre 
    ## la matriz Ac
            
    
    return L, U, cant_op

'''

############### res_tri: Mision Cumplida

def res_tri(L, b, inferior=True):

    n = L.shape[0]
    x = np.zeros(n)

    if inferior:
        # Sustitución hacia adelante
        for i in range(n):
            if L[i, i] == 0:
                return None

            x[i] = b[i]

            for j in range(i):
                x[i] = x[i] - L[i, j] * x[j]

            x[i] = x[i] / L[i, i]

    else:
        # Sustitución hacia atrás
        for i in range(n - 1, -1, -1):
            if L[i, i] == 0:
                return None

            x[i] = b[i]

            for j in range(i + 1, n):
                x[i] = x[i] - L[i, j] * x[j]

            x[i] = x[i] / L[i, i]

    return x


###########inversa: Mision Cumplida

def inversa(A):

    resultado = calculaLU(A)

    if resultado is None:
        return None

    L, U, nops = resultado
    n = A.shape[0]

    Ainv = np.zeros((n, n))

    for j in range(n):

        # Columna j de la matriz identidad
        b = np.zeros(n)
        b[j] = 1

        # L y = e_j
        y = res_tri(L, b)

        if y is None:
            return None

        # U x = y
        x = res_tri(U, y, inferior=False)

        if x is None:
            return None

        # La solución es la columna j de A^(-1)
        for i in range(n):
            Ainv[i, j] = x[i]

    return Ainv

def calculaLU(A):

    if A is None:
        return None

    n = A.shape[0]

    U = A.copy()
    nops = 0

    for k in range(n - 1):

        pivote = U[k, k]

        if pivote == 0:
            return None

        for i in range(k + 1, n):

            m = U[i, k] / pivote
            nops += 1

            U[i, k] = m

            U[i, k + 1:n] = U[i, k + 1:n] - m * U[k, k + 1:n]
            nops += 2 * (n - k - 1)

    L = np.eye(n)
    U_final = np.zeros((n, n))

    for i in range(n):
        for j in range(n):
            if i > j:
                L[i, j] = U[i, j]
            else:
                U_final[i, j] = U[i, j]

    return L, U_final, nops


def res_tri(L, b, inferior=True):

    n = L.shape[0]
    x = np.zeros(n)

    if inferior:
        # Sustitución hacia adelante
        for i in range(n):
            if L[i, i] == 0:
                return None

            x[i] = b[i]

            for j in range(i):
                x[i] = x[i] - L[i, j] * x[j]

            x[i] = x[i] / L[i, i]

    else:
        # Sustitución hacia atrás
        for i in range(n - 1, -1, -1):
            if L[i, i] == 0:
                return None

            x[i] = b[i]

            for j in range(i + 1, n):
                x[i] = x[i] - L[i, j] * x[j]

            x[i] = x[i] / L[i, i]

    return x



###########calculaLDV : Mision Cumplida
def calculaLDV(A):
    """
    Calcula la factorización LDV de A:
    
        A = L D V
    
    donde L es triangular inferior,
    D diagonal y V triangular superior.
    """

    L, U, nops = calculaLU(A)

    if L is None:
        return None, None, None

    n = A.shape[0]

    D = np.zeros((n, n))
    V = np.zeros((n, n))

    for i in range(n):

        if U[i, i] == 0:
            return None, None, None

        D[i, i] = U[i, i]

        for j in range(i, n):
            V[i, j] = U[i, j] / D[i, i]

    return L, D, V


def calculaLU(A):

    if A is None:
        return None

    n = A.shape[0]

    U = A.copy()
    nops = 0

    for k in range(n - 1):

        pivote = U[k, k]

        if pivote == 0:
            return None

        for i in range(k + 1, n):

            m = U[i, k] / pivote
            nops += 1

            U[i, k] = m

            U[i, k + 1:n] = U[i, k + 1:n] - m * U[k, k + 1:n]
            nops += 2 * (n - k - 1)

    L = np.eye(n)
    U_final = np.zeros((n, n))

    for i in range(n):
        for j in range(n):
            if i > j:
                L[i, j] = U[i, j]
            else:
                U_final[i, j] = U[i, j]

    return L, U_final, nops



############ esSDP: No se pueden correr los test

def esSDP(A, atol=1e-8):
    """
    Comprueba si A es simétrica definida positiva
    utilizando la factorización LDV.
    """

    n = A.shape[0]

    # Verificar simetría
    for i in range(n):
        for j in range(i + 1, n):

            if abs(A[i, j] - A[j, i]) > atol:
                return False

    # Calcular LDV
    L, D, V = calculaLDV(A)

    if L is None:
        return False

    # Verificar diagonal positiva
    for i in range(n):

        if D[i, i] <= atol:
            return False

    return True

def calculaLDV(A):
    """
    Calcula la factorización LDV de A:
    
        A = L D V
    
    donde L es triangular inferior,
    D diagonal y V triangular superior.
    """

    L, U, nops = calculaLU(A)

    if L is None:
        return None, None, None

    n = A.shape[0]

    D = np.zeros((n, n))
    V = np.zeros((n, n))

    for i in range(n):

        if U[i, i] == 0:
            return None, None, None

        D[i, i] = U[i, i]

        for j in range(i, n):
            V[i, j] = U[i, j] / D[i, i]

    return L, D, V


def calculaLU(A):

    if A is None:
        return None

    n = A.shape[0]

    U = A.copy()
    nops = 0

    for k in range(n - 1):

        pivote = U[k, k]

        if pivote == 0:
            return None

        for i in range(k + 1, n):

            m = U[i, k] / pivote
            nops += 1

            U[i, k] = m

            U[i, k + 1:n] = U[i, k + 1:n] - m * U[k, k + 1:n]
            nops += 2 * (n - k - 1)

    L = np.eye(n)
    U_final = np.zeros((n, n))

    for i in range(n):
        for j in range(n):
            if i > j:
                L[i, j] = U[i, j]
            else:
                U_final[i, j] = U[i, j]

    return L, U_final, nops


############################ cholesky: Mision Cumplida
def calculaCholesky(A, atol=1e-10):

    if not esSDP(A, atol):
        return None

    resultado = calculaLDV(A)

    if resultado is None:
        return None

    L, D, V = resultado

    n = np.shape(A)[0]

    R = np.zeros((n,n))

    for j in range(n):

        raiz = np.sqrt(D[j,j])

        for i in range(n):

            R[i,j] = L[i,j] * raiz

    return R

def esSDP(A, atol=1e-8):
    """
    Comprueba si A es simétrica definida positiva
    utilizando la factorización LDV.
    """

    n = A.shape[0]

    # Verificar simetría
    for i in range(n):
        for j in range(i + 1, n):

            if abs(A[i, j] - A[j, i]) > atol:
                return False

    # Calcular LDV
    L, D, V = calculaLDV(A)

    if L is None:
        return False

    # Verificar diagonal positiva
    for i in range(n):

        if D[i, i] <= atol:
            return False

    return True

def calculaLDV(A):
    """
    Calcula la factorización LDV de A:
    
        A = L D V
    
    donde L es triangular inferior,
    D diagonal y V triangular superior.
    """

    L, U, nops = calculaLU(A)

    if L is None:
        return None, None, None

    n = A.shape[0]

    D = np.zeros((n, n))
    V = np.zeros((n, n))

    for i in range(n):

        if U[i, i] == 0:
            return None, None, None

        D[i, i] = U[i, i]

        for j in range(i, n):
            V[i, j] = U[i, j] / D[i, i]

    return L, D, V


def calculaLU(A):

    if A is None:
        return None

    n = A.shape[0]

    U = A.copy()
    nops = 0

    for k in range(n - 1):

        pivote = U[k, k]

        if pivote == 0:
            return None

        for i in range(k + 1, n):

            m = U[i, k] / pivote
            nops += 1

            U[i, k] = m

            U[i, k + 1:n] = U[i, k + 1:n] - m * U[k, k + 1:n]
            nops += 2 * (n - k - 1)

    L = np.eye(n)
    U_final = np.zeros((n, n))

    for i in range(n):
        for j in range(n):
            if i > j:
                L[i, j] = U[i, j]
            else:
                U_final[i, j] = U[i, j]

    return L, U_final, nops









###########################################

def main():
    n = 7
    B = np.eye(n) - np.tril(np.ones((n,n)),-1) 
    B[:n,n-1] = 1
    print('Matriz B \n', B)
    
    L,U,cant_oper = elim_gaussiana(B)
    
    print('Matriz L \n', L)
    print('Matriz U \n', U)
    print('Cantidad de operaciones: ', cant_oper)
    print('B=LU? ' , 'Si!' if np.allclose(np.linalg.norm(B - L@U, 1), 0) else 'No!')
    print('Norma infinito de U: ', np.max(np.sum(np.abs(U), axis=1)) )

if __name__ == "__main__":
    main()
    
    
