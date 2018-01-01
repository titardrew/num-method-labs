import numpy as np
from scipy import special as sc
from itertools import chain
from collections import defaultdict
from matplotlib import mlab
import bisect, pylab

a = 5
b = 10
knot_step = 0.5
error_step = knot_step / 5
sup = 1.67253
knot_points = mlab.frange(a, b, knot_step)
splines = defaultdict(lambda: Tuple())
error_set = mlab.frange(a, b, error_step)


def major(x):
    return sup / 720 * np.product([abs(x - p) for p in knot_points])


def f(x):
    return np.log(np.log(np.log(x)**3)**2)


def straight_difference(k, r, knot_points):
    return sum([(-1)**j * sc.binom(r, j) * f(knot_points[k+r-j]) for j in range(r+1)])


def lagrange_poly(x, knot_points):
    return sum(f(knot_points[i]) * np.product([(x - knot_points[j]) / (knot_points[i] - knot_points[j])
                                               for j in chain(range(0, i), range(i + 1, len(knot_points)))])
               for i in range(len(knot_points)))


def newton_poly_straight(x, knot_points):
    return f(knot_points[0]) + sum([sum([f(knot_points[j]) / np.product([((knot_points[j] - knot_points[k])
                                                                          if k != j else 1) for k in range(0, i+1)])
                                   for j in range(0, i+1)]) * np.prod([x - knot_points[k] for k in range(i)])
                                    for i in range(1, len(knot_points))])


def newton_poly_backward(x, knot_points):
    knot_points = [knot_points[len(knot_points) - 1 - i] for i in range(len(knot_points))]
    return f(knot_points[0]) + sum([sum([f(knot_points[j]) / np.product([((knot_points[j] - knot_points[k])
                                                                          if k != j else 1) for k in range(0, i+1)])
                                   for j in range(0, i+1)]) * np.prod([x - knot_points[k] for k in range(i)])
                                    for i in range(1, len(knot_points))])


class Tuple:
    a, b, c, d, x = [0., 0., 0., 0., 0.]


def build_spline(knot_points):
    for i in range(len(knot_points)):
        splines[i].x, splines[i].a = knot_points[i], f(knot_points[i])

    alpha, beta = [defaultdict(lambda: 0.), defaultdict(lambda: 0.)]

    for i in range(1, len(knot_points)-1):
        C = 4. * knot_step
        F = 6. * ((f(knot_points[i + 1]) - f(knot_points[i])) / knot_step - (f(knot_points[i]) - f(knot_points[i - 1])) / knot_step)
        z = (knot_step * alpha[i - 1] + C)
        alpha[i] = -knot_step / z
        beta[i] = (F - knot_step * beta[i - 1]) / z

    for i in reversed(range(1, len(knot_points) - 1)):
        splines[i].c = alpha[i] * splines[i+1].c + beta[i]

    for i in reversed(range(1, len(knot_points))):
        hi = knot_points[i] - knot_points[i - 1]
        splines[i].d = (splines[i].c - splines[i-1].c) / hi
        splines[i].b = hi * (2.0 * splines[i].c + splines[i - 1].c) / 6.0 + (f(knot_points[i]) - f(knot_points[i - 1])) / hi


def cubic_spline(x):
    distribution = sorted([t[1].x for t in splines.items()])
    indx = bisect.bisect_left(distribution, x)
    if indx == len(distribution): return 0
    dx = x - splines[indx].x
    return splines[indx].a + splines[indx].b * dx + splines[indx].c * dx**2 / 2. + splines[indx].d * dx**3 / 6.


def log(i, value, expected):
    print("x =", str(i) + ";", "F(x) =", expected, "P(x) =", value,  "Error =", abs(value - expected), "Major =", major(i))


def main():
    build_spline(knot_points)

    print("\n\n\t\t\t\t\tLAGRANGE METHOD\n")
    for i in range(len(error_set)):
        log(error_set[i], lagrange_poly(error_set[i], knot_points), f(error_set[i]))

    print("\n\n\t\t\t\t\tNEWTON STRAIGHT METHOD\n")
    for i in range(len(error_set)):
        log(error_set[i], newton_poly_straight(error_set[i], knot_points), f(error_set[i]))

    print("\n\n\t\t\t\t\tNEWTON BACKWARD METHOD\n")
    for i in range(len(error_set)):
        log(error_set[i], newton_poly_backward(error_set[i], knot_points), f(error_set[i]))

    print("\n\n\t\t\t\t\tCUBIC SPLINE METHOD\n")
    for i in range(len(error_set)):
        log(error_set[i], cubic_spline(error_set[i]), f(error_set[i]))


def plot_error():
    pylab.plot(error_set, [major(x) for x in error_set])
    pylab.plot(error_set, [abs(f(x) - cubic_spline(x)) for x in error_set])
    # pylab.plot(error_set, [abs(f(x) - lagrange_poly(x, knot_points)) for x in error_set])
    # pylab.plot(error_set, [abs(f(x) - newton_poly_straight(x, knot_points)) for x in error_set])
    # pylab.plot(error_set, [f(x) - newton_poly_backward(x, knot_points) for x in error_set])
    pylab.plot()
    pylab.show()


def plot_interpolation():
    pylab.plot(error_set, [f(x) for x in error_set])
    pylab.plot(error_set, [cubic_spline(x) for x in error_set])
    pylab.plot(error_set, [lagrange_poly(x, knot_points) for x in error_set])
    pylab.plot()
    pylab.show()

main()
plot_error()
# plot_interpolation()



