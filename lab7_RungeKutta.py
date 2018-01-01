from numpy import log


def rk4(f, x0, y0, x1, n):
    vx = [0] * (n + 1)
    vy = [0] * (n + 1)
    h = (x1 - x0) / float(n)
    vx[0] = x = x0
    vy[0] = y = y0
    for i in range(1, n + 1):
        k1 = h * f(x, y)
        k2 = h * f(x + 0.5 * h, y + 0.5 * k1)
        k3 = h * f(x + 0.5 * h, y + 0.5 * k2)
        k4 = h * f(x + h, y + k3)
        vx[i] = x = x0 + i * h
        vy[i] = y = y + (k1 + k2 + k2 + k3 + k3 + k4) / 6
    return vx, vy


def f(x, y):  # Andrew
    return (x + 2*y - 3) / (x - 1)

def sol(x):
    return (x - 1)**2 - x + 2

# def f(x, y):  # Inna
#     return (2*x + y - 3) / (x-1)
#
#
# def sol(x):
#     return (2 * log(x - 1) + 1) * (x-1) + 1






# vx, vy = rk4(f_inna, 2, 2, 2+1, 20)
# for x, y in list(zip(vx, vy)):
#     print("%4.2f %10.5f %+12.4e" % (x, y, y - ((2*log(x - 1) + 1) * (x-1) + 1) ))