from enum import nonmember
from math import sqrt
from copy import deepcopy as dcpy

def separate_augmented(Aaug):
    """
    Splits matrix into A and b
    :param Aaug: augmented matrix
    :return: (A, b)
    """

    n = len(Aaug)
    A = []
    b = []
    for i in range(n):
        A.append([Aaug[i][j] for j in range(n)])
        b.append(Aaug[i][n])
    return A, b

def LUFactorization(A):
    """
    Lower upper factorization part of Doolittles method
    :param A: nxn matrix
    :return: a tuple w/ (L, U)
    """

    n = len(A)

    U = [([0 for c in range(n)] if not r == 0 else [a for a in A[0]]) for r in range(n)]
    L = [[(1 if c == r else (A[r][0] / U[0][0] if c == 0 else 0)) for c in range(n)] for r in range(n)]

    for j in range(1, n):
        for k in range(j, n):
            U[j][k] = A[j][k]
            for s in range(j):
                U[j][k] -= L[j][s] * U[s][k]

            for i in range(k + 1, n):
                sig = 0
                for s in range(k):
                    sig += L[i][s] * U[s][k]
                L[i][k] = (1 / (U[k][k])) * (A[i][k] - sig)
    return (L, U)

def BackSolve(A, b, UT=True):
    """
    backsolving triangular matrix
    :param A: triangular matrix
    :param b: vector
    :param UT: upper triangular matrix
    :return: solution vector
    """

    nRows = len(b)
    x = [0] * nRows

    if UT:
        for nR in range(nRows - 1, -1, -1):
            s = 0
            for nC in range(nR + 1, nRows):
                s += A[nR][nC] * x[nC]
            x[nR] = 1 / A[nR][nR] * (b[nR] - s)

    else:
        for nR in range(nRows):
            s = 0
            for nC in range(nR):
                s += A[nR][nC] * x[nC]
            x[nR] = 1 / A[nR][nR] * (b[nR] - s)

    return x

def Doolittle(Aaug):
    """
    Doolittle method for solving ax=b using LU factorization
    :param Aaug: augmented matrix
    :return: solution vector
    """
    A, b = separate_augmented(Aaug)
    L, U = LUFactorization(A)
    y = BackSolve(A, b, UT=False)
    x = BackSolve(A, b, UT=True)
    return x

def is_symmetric(A, tol=1e-10):
    """
    check if matrix is symmetric
    :param A: square matrix
    :param tol: tolerance
    :return:
    """

    n = len(A)
    for i in range(n):
        for j in range(i + 1, n):
            if abs(A[i][j] - A[j][i]) > tol:
                return False
    return True

def cholesky_factor(A):
    n = len(A)
    L = [[0.0 for _ in range(n)] for _ in range(n)]

    for i in range(n):
        for j in range(i + 1):
            s = 0.0
            for k in range(j):
                s += L[i][k] * L[j][k]

            if i == j:
                val = A[i][i] - s
                if val <= 0.0:
                    return None
                L[i][j] = sqrt(val)
            else:
                if L[j][j] == 0.0:
                    return None
                L[i][j] = (A[i][j] - s) / L[j][j]

    return L
def forward_sub(L, b):
    """
    Solves forward sub
    :param L: lower triangular matrix
    :param b: vector
    :return: y solution vector
    """

    n = len(L)
    y = [0.0] * n
    for i in range(n):
        rhs = b[i]
        for j in range(i):
            rhs -= L[i][j] * y[j]
        y[i] = rhs / L[i][i]
    return y

def back_sub(U, y):
    """
    Solves back sub
    :param U: upper triangular matrix
    :param y: vector
    :return: x solution vector
    """

    n = len(U)
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        rhs = y[i]
        for j in range(i + 1, n):
            rhs -= U[i][j] * x[j]
        x[i] = rhs / U[i][i]
    return x

def transpose(A):
    """
    Transpose matrix
    :param A: matrix
    :return: transposed matrix
    """

    n = len(A)
    m = len(A[0])
    AT = [[0 for _ in range(n)] for _ in range(m)]
    for i in range(n):
        for j in range(m):
            AT[j][i] = A[i][j]
    return AT

def solve_cholesky_or_none(A, b):
    """
    if positive and definate, solve using cholesky
    :param A: matrix
    :param b: vector
    :return: x solution vector or none
    """

    if not is_symmetric(A):
        return None

    L = cholesky_factor(A)
    if L is None:
        return None

    y = forward_sub(L, b)
    LT = transpose(L)
    x = back_sub(LT, y)
    return x

def pretty_print(x):
    """
    Pretty printing function
    :param x: vector
    :return: none
    """

    for i, val in enumerate(x):
        print("x{} = {:0.6f}".format(i + 1, val))

def make_aug(A, b):
    """
    Make augmented matrix
    :param A: matric
    :param b: vector
    :return: augmented matrix
    """

    n = len(A)
    Aaug = []
    for i in range(n):
        Aaug.append([A[i][j] for j in range(n)] + [b[i]])

    return Aaug

def main():
    """
    Solves using cholesky, else, solves using doolittle
    :return:
    """

    print("HW3c: Cholesky if SPD, else Doolittle\n")

    A1 = [
        [1, -1, 3, 2],
        [-1, 5, -5, -2],
        [3, -5, 19, 3],
        [2, -2, 3, 21]
    ]
    b1 = [15, -35, 94, 1]

    x1 = solve_cholesky_or_none(A1, b1)
    if x1 is not None:
        print("System 1 used Cholesky")
        pretty_print(x1)
    else:
        print("System 1 used Doolittle")
        x1 = Doolittle(make_aug(A1, b1))
        pretty_print(x1)

    print("\n-------------------------\n")

    A2 = [
        [4, 2, 4, 0],
        [2, 2, 3, 2],
        [4, 3, 6, 3],
        [0, 2, 3, 9]
    ]
    b2 = [20, 36, 60, 122]

    x2 = solve_cholesky_or_none(A2, b2)
    if x2 is not None:
        print("System 2 used Cholesky")
        pretty_print(x2)
    else:
        print("System 2 used Doolittle")
        x2 = Doolittle(make_aug(A2, b2))
        pretty_print(x2)


if __name__ == "__main__":
    main()

#Dr. Kolla's code provided a basis and ChatGPT was consulted for logic and correctness