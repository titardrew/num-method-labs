import numpy as np


def power_iteration(operator, eps=1.0e-5):

    b_k = np.ones(operator.shape[0])
    b_k1 = np.dot(operator, b_k)

    eigenvalue = b_k1[0] / b_k[0]

    # calculate the norm
    b_k1_norm = np.linalg.norm(b_k1)

    # re normalize the vector
    b_k = b_k1 / b_k1_norm

    iteration = 0

    while True:
        iteration += 1

        # calculate the matrix-by-vector product Ab
        b_k1 = np.dot(operator, b_k)

        # calculate the norm
        b_k1_norm = np.linalg.norm(b_k1)

        # re normalize the vector
        b_k = b_k1 / b_k1_norm

        print('{:3}\t eigenvalue = {:20},\terror = {:22}, residual = {}'
              .format(iteration,
                      eigenvalue,
                      abs(eigenvalue - b_k1[0] / b_k[0]),
                      (np.dot(operator, b_k) - b_k * b_k1[0] / b_k[0])
                      )
              )

        cond = True
        for i in range(len(b_k)):
            cond = cond and -eps < (np.dot(operator, b_k) - b_k * b_k1[0] / b_k[0])[i] < eps

        if abs(eigenvalue - b_k1[0] / b_k[0]) < eps and cond:
            break
        else:
            eigenvalue = b_k1[0] / b_k[0]

    return b_k, b_k1[0] / b_k[0]


def main():

    operator = np.array(np.loadtxt("lab6_input", comments="#", delimiter=" ", unpack=False))

    # calc max eigenvalue & eigenvector
    eigenvector, eigenvalue = power_iteration(operator)
    print('\nresidual = {}\neigenvector = {}\neigenvalue = {}\n'
          .format(np.dot(operator, eigenvector)
                  - eigenvalue * eigenvector,
                  eigenvector,
                  eigenvalue)
          )

    # calc min eigenvalue & eigenvector
    eigenvector, eigenvalue = power_iteration(np.linalg.inv(operator))
    eigenvalue = 1 / eigenvalue
    print('\nresidual = {}\neigenvector = {}\neigenvalue = {}\n'
          .format(np.dot(operator, eigenvector)
                  - eigenvalue * eigenvector,
                  eigenvector,
                  eigenvalue)
          )


if __name__ == '__main__':
    main()
