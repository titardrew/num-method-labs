import numpy as np
from numpy import sin, pi, exp, ndarray, log, linspace
from openpyxl import Workbook
from matplotlib import pyplot as plt
import os

A = 1
q = 1
L = 1
dx = 0.1
x_size = int(L/dx + 1)
dt = 0.1
t_size = int(1/dt + 1)
u = ndarray((t_size, x_size))
nplot = 1


def u0(x):
    return A * sin(pi * x / L)


def find_alpha():
    s = 0
    c = 0
    for i_t in range(1, t_size):
        j_x = 1
        while j_x < (x_size+1)/2:
            s += get_alpha(i_t, j_x)
            c += 1
            j_x += 1
    return s / c


def get_alpha(t, x):
    if (t > 0) and (x != x_size - 1) and (x != 0):
        return -1/float(t)/dt*log(u[t][x]/A/sin(pi*x*dx/L))
    else:
        return 0


def e(t, x, alpha):
    return dx*dt*(A * exp(-alpha * t * dt) * sin(x * dx * pi / L) - u[t][x])


def triagonal_solve2(n, a, b, c, d):
    p = ndarray(n)
    q = ndarray(n)
    x = ndarray(n)

    p[0] = b[0]/a[0]
    q[0] = d[0]/a[0]

    for i in range(1, n):
        p[i] = b[i] / (a[i] - c[i] * p[i-1])
        q[i] = (c[i]*q[i-1]+d[i])/(a[i]-c[i]*p[i-1])

    x[n-1] = q[n - 1]
    for i in reversed(range(n - 1)):
        x[i] = p[i] * x[i+1] + q[i]

    return x


def fill_u():
    a = ndarray(x_size)
    b = ndarray(x_size)
    c = ndarray(x_size)
    d = ndarray(x_size)
    c_ = 0
    x = linspace(0, L, x_size)

    for i in range(x_size):
        u[0][i] = u0(dx * i)

    fig = plt.figure()
    plt.plot(x, u[0], linewidth=2)
    filename = 'foo000.jpg'
    fig.set_tight_layout(True)
    plt.ylim([0, 1100])
    plt.xlabel("x")
    plt.ylabel("u")
    plt.title("t = 0")
    plt.savefig(filename)
    plt.clf()

    for i in range(1, t_size):
        a[0], b[0], c[0] = 1, 0, 0
        d[0] = u[i - 1][0] * (1 + q * dt)

        for j in range(1, x_size - 1):
            a[j] = 1 + 2 * A * dt / (dx**2)
            b[j] = A * dt / (dx**2)
            c[j] = A * dt / (dx**2)
            d[j] = u[i - 1][j] * (1 + q * dt)
            print(a[j], b[j], c[j], d[j])

        a[x_size - 1] = 1
        b[x_size - 1] = 0
        c[x_size - 1] = 0
        d[x_size - 1] = u[i - 1][x_size - 1]*(1 + q * dt)

        u[i] = triagonal_solve2(x_size, a, b, c, d)

        if i % nplot == 0:  # plot results every nplot timesteps

            plt.plot(x, u[i], linewidth=2)
            plt.ylim([0, 1100])
            filename = 'foo' + str(c_ + 1).zfill(3) + '.jpg'
            plt.xlabel("x")
            plt.ylabel("u")
            plt.title("t = %2.2f" % (dt * (i + 1)))
            plt.savefig(filename)
            plt.clf()
            c_ += 1


fill_u()

wb = Workbook()

# grab the active worksheet
ws = wb.active
out_rows = []
row = [' ']
for i in range(x_size):
    if (i) % int(1/ dx / 10) == 0:
        row.append(dx * i)
        print('{}'.format(dx * i), end='\t')

out_rows.append(row)

print()
for i in range(t_size):
    row = []
    row.append(i*dt)
    print(i * dt, end=' ')
    for j in range(x_size):
        if (j) % int(1/ dx / 10) == 0:
            row.append(u[i][j])
            print(u[i][j], end=' ')
    out_rows.append(row)

    print()

print()

for row_ in out_rows:
    ws.append(row_)

ws.append([' '])

out_rows = []
row = [' ']

for i in range(x_size):
    if (i) % int(1/ dx / 10) == 0:
        row.append(dx * i)
        print('{}'.format(dx * i), end='\t')
print()
print(find_alpha())
e_graph = []
out_rows.append(row)

for i in range(t_size):
    row = []
    row.append(i * dt)
    print(i * dt, end=' ')
    print(find_alpha())
    err = 0
    for j in range(x_size):
        if (j) % int(1/ dx / 10) == 0:
            row.append(e(i, j, find_alpha()))
            print(e(i, j, find_alpha()), end=' ')
    e_graph.append(np.sqrt(sum([(e(i, j, find_alpha())*dx*dt)**2 for j in range(x_size)])))
    out_rows.append(row)
    print()

for row_ in out_rows:
    ws.append(row_)
wb.save('output_9.xls')

plt.plot(linspace(0, 1, 1 / dt + 1), e_graph)
plt.show()


os.system("ffmpeg -y -i 'foo%03d.jpg' heat_equation.m4v")
os.system("rm -f *.jpg")
