import numpy as np

from numpy import identity, diagonal
from math import sqrt


def max_elem(a):  # Find largest off-diag. element a[k,l]

    n = len(a)
    a_max = 0.0
    k, l = 0, 0
    for i in range(n - 1):
        for j in range(i + 1, n):
            if abs(a[i][j]) >= a_max:
                a_max = abs(a[i][j])
                k = i
                l = j
    return a_max, k, l


def rotate(a, p, k, l):  # Rotate to make a[k,l] = 0

    n = len(a)
    a_diff = a[l][l] - a[k][k]

    phi = a_diff / (2.0 * a[k][l])
    t = 1.0 / (abs(phi) + sqrt(phi**2 + 1.0))
    if phi < 0.0:
        t = -t

    c = 1.0 / sqrt(t**2 + 1.0)
    s = t * c
    tau = s / (1.0 + c)
    temp = a[k][l]

    a[k][l] = 0.0
    a[k][k] -= t * temp
    a[l][l] += t * temp

    for j in range(k):  # Case of i < k
        temp = a[j][k]
        a[j][k] = temp - s * (a[j][l] + tau * temp)
        a[j][l] += s * (temp - tau * a[j][l])

    for j in range(k + 1, l):  # Case of k < i < l
        temp = a[k][j]
        a[k][j] = temp - s * (a[j][l] + tau * a[k][j])
        a[j][l] += s * (temp - tau * a[j][l])

    for j in range(l + 1, n):  # Case of i > l
        temp = a[k][j]
        a[k][j] = temp - s * (a[l][j] + tau * temp)
        a[l][j] += s * (temp - tau * a[l][j])

    for j in range(n):  # Update transformation matrix
        temp = p[j][k]
        p[j][k] = temp - s * (p[j][l] + tau * p[j][k])
        p[j][l] += s * (temp - tau * p[j][l])
    return c, s, t, phi


def jacobi(a, eps=1.0e-5):  # Jacobi method

    n = len(a)

    # Set limit on number of rotations
    rotation_limit = 5 * (n**2)

    # Initialize transformation matrix
    p = identity(n) * 1.0
    a_ = p.copy()
    # Jacobi rotation loop
    for i in range(rotation_limit):
        print('\n\n\tITERATION {}\n'.format( i+1))
        print('Matrix: \n', a, '\n')

        w2 = 0
        d = 0
        for i1 in range(len(a)):
            for j1 in range(len(a)):
                d += a[i1][j1]**2 if i1 == j1 else 0
                w2 += a[i1][j1]**2 if i1 < j1 else 0

        print('d = {:20}, w2 = {:20}, d + w2 = {:20}'.format(d, w2, d + 2*w2))
        a_max, k, l = max_elem(a)
        print('i_max = {:4}, j_max = {:4}, Aij = {:20}, Aii = {:20}, Ajj = {:20}'.format(k + 1, l + 1, a[k][l], a[k][k], a[l][l]))
        if a_max < eps:
            return diagonal(a), p
        c, s, t, ksi = rotate(a, p, k, l)
        if i == 1:
            a_ = p.copy()
        else:
            a_ = a_.dot(p)
        print('s = {:20}, c = {:20}, t = {:20}, ksi = {:20}, c^2+s^2 = {:2}'.format(s, c, t, ksi, c**2 + s**2))

    print('Jacobi method did not converge')


def main():

    operator = np.array(np.loadtxt("lab6_input", comments="#", delimiter=" ", unpack=False))
    operator_ = operator.copy()
    eigenvalues, eigenvectors = jacobi(operator_)

    print('eigenvalues: \n', eigenvalues, '\neigenvectors: ')
    eigenvector = np.ndarray((len(eigenvectors), len(eigenvectors)))

    for i in range(len(eigenvectors)):
        for j in range(len(eigenvectors)):
            eigenvector[i][j] = eigenvectors[j][i]

    print(eigenvector)
    print('Residuals :')
    for vec, val in list(zip(eigenvector, eigenvalues)):
        print(operator.dot(vec) - vec * val)

if __name__ == '__main__':
    main()