import numpy as np
import numpy.linalg as lg

# Simple iteration method & Newton iteration method for both systems (1), (2)


def simple_iteration(pt0, coef, eps, iteration=1):
    pt = np.array([coef[5] - np.sin(pt0[1] + coef[1]), (coef[4] - np.cos(pt0[0] + coef[0])) / coef[3]])
    residual = np.array([pt[0] - coef[5] + np.sin(pt[1] + coef[1]),
                         coef[3] * pt[1] - coef[4] + np.cos(pt0[0] + coef[0])])
    log(iteration, pt, residual)
    return simple_iteration(pt , coef, eps, iteration=iteration+1) if lg.norm(residual) > eps \
        and lg.norm(pt - pt0) > eps else pt


def log(iteration, pt, residual):
    print("iteration ", iteration, ": point = ", pt, "\tresidual = ", residual, "\t norm = ", lg.norm(residual))


def log_system1(coef):
    print("\ncos( x + ", coef[0], ") + ", coef[3], "* y = ", coef[4], "\n",
          "x + sin( y + ", coef[1], ") = ", coef[5])


def newton_iteration(pt0, coef, eps, iteration=1):
    if lg.det(dFdxdy(pt0, coef)) == 0:
        raise Exception
    pt = pt0 - lg.inv(dFdxdy(pt0, coef)).dot(F(pt0, coef))
    log(iteration, pt, F(pt, coef))
    return newton_iteration(pt, coef, eps, iteration=iteration+1) if lg.norm(F(pt, coef)) > eps \
        and lg.norm(pt - pt0) > eps else pt


def log_system2(coef):
    print("\nsin( x + y ) + ", coef[4], "* x = ", coef[5], "\n",
          "( x^2 + y^2 ) = 1")


def F(pt, coef):
    return np.array([np.sin(pt[0] + pt[1]) + coef[4] * pt[0] - coef[5], pt[0]**2 + pt[1]**2 - 1])


def dFdxdy(pt, coef):
    return np.array([[np.cos(pt[0] + pt[1]) + coef[4], np.cos(pt[0] + pt[1])],
            [2 * pt[0], 2 * pt[1]]])


log_system1([0.037, -0.232,	0.808, -1.480, 1.817, 0.224])
print("\n\n\t\t\t\t\t\t\t\t\t\t\t\t\tSIMPLE ITERATION METHOD\n")
x, y = simple_iteration(np.array([1, -0.5]), [0.037, -0.232, 0.808, -1.480, 1.817, 0.224], 0.00001)
print("\nx =", x, "y =", y, end='\n\n\n\n\n')

log_system2([0.037, -0.232,	0.808, -1.480, 1.817, 0.224])
print("\n\n\t\t\t\t\t\t\t\t\t\t\t\t\tNEWTON ITERATION METHOD\n")
try:
    print("point 1: ")
    x1, y1 = newton_iteration(np.array([0, 2]), [0.037, -0.232, 0.808, -1.480, 1.817, 0.224], 0.00001)
    print("\npoint 2: ")
    x2, y2 = newton_iteration(np.array([1, -1]), [0.037, -0.232, 0.808, -1.480, 1.817, 0.224], 0.00001)
    print("\nx1 =", x1, "y1 =", y1)
    print("x2 =", x2, "y2 =", y2)
except Exception:
    print("Can't be processed...")