import numpy as np

a = 1.878
b = -1.606
c = 2.071
d = -2.839
e = -0.320
boundary_coef11 = 1
boundary_coef12 = -0.3
boundary_coef21 = 0
boundary_coef22 = 1
left = 0.6
right = 0.9


def y(x):
    return a*x**2 + b*x + c + 1 / (d*x + e)


def dydx(x):
    return 2*a*x + b - d/(d*x + e)**2


def dy2dx(x):
    return 2*a + 2*d**2/(d*x + e)**3


def p(x):
    return 1 / x


def q(x):
    return -0.4


def f(x):
    return dy2dx(x) + p(x) * dydx(x) + q(x) * y(x)


def triagonal_solve(a_, b_, c_, d_, n_):
    x = np.ndarray(n_)
    c_[0] /= b_[0]
    d_[0] /= b_[0]

    for i in range(1, n_):
        id_ = 1 / (b_[i] - c_[i - 1] * a_[i])
        c_[i] *= id_
        d_[i] = (d_[i] - d_[i-1] * a_[i]) * id_

    x[n_ - 1] = d_[n_ - 1]
    for i in reversed(range(0, n_ - 1)):
        x[i] = d_[i] - c_[i] * x[i + 1]

    return x


def finite_diff(n):

    left_cond = boundary_coef11 * y(left) + boundary_coef12 * dydx(left)
    right_cond = boundary_coef21 * y(right) + boundary_coef22 * dydx(right)

    h = (right - left) / n

    x0 = left - h/2
    x_net = np.ndarray(n + 2)
    y_net = np.ndarray(n + 2)

    f_net = np.ndarray(n + 2)
    for i in range(n + 2):
        x_net[i] = x0
        y_net[i] = y(x_net[i])
        f_net[i] = f(x_net[i])
        x0 += h

    a_ = np.ndarray(n + 2)
    b_ = np.ndarray(n + 2)
    c_ = np.ndarray(n + 2)
    d_ = np.ndarray(n + 2)

    for i in range(1, n + 1):
        a_[i] = 1 / h / h - p(x_net[i]) / 2 / h
        b_[i] = q(x_net[i]) - 2 / h / h
        c_[i] = 1 / h / h + p(x_net[i]) / 2 / h
        d_[i] = f_net[i]

    a_[0] = 0
    b_[0] = boundary_coef11 / 2 - boundary_coef12 / h
    c_[0] = boundary_coef11 / 2 + boundary_coef12 / h
    d_[0] = left_cond
    a_[n + 1] = boundary_coef21 / 2 - boundary_coef22 / h
    b_[n + 1] = boundary_coef21 / 2 + boundary_coef22 / h
    c_[n + 1] = 0
    d_[n + 1] = right_cond

    r = np.ndarray(n + 2)
    r = triagonal_solve(a_, b_, c_, d_, n + 2)

    e = np.ndarray(n + 2)
    norm = 0
    for i in range(n + 2):
        e[i] = abs(y_net[i] - r[i])
        norm = abs(e[i]) if norm < abs(e[i]) else norm
        print('x[i] = {:20}, numeric[i] = {:20}, exact[i] = {:20}, error[i] = {:20}'
              .format(x_net[i], r[i], y_net[i], e[i]))
    print('||e|| = {:20}'.format(norm))
    return norm


def main():
    n = 30
    finite_diff(n)

    # n = 220
    # err_h = []
    # n_s = []
    # left = 0.6
    # right = 0.9
    # for i in range(10):
    #     n_s.append((right - left) / n)
    #     err_h.append(finite_diff(n))
    #     n -= 20
    # from matplotlib import pyplot as plt
    # plt.plot(n_s, err_h, marker='o')
    # plt.show()

if __name__ == '__main__':
    main()