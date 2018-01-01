from lab7_RungeKutta import rk4, f, sol
from pandas import DataFrame

def kutta_light(f, x, y, h):
    k1 = h * f(x, y)
    k2 = h * f(x + 0.5 * h, y + 0.5 * k1)
    k3 = h * f(x + 0.5 * h, y + 0.5 * k2)
    k4 = h * f(x + h, y + k3)
    return y + (k1 + k2 + k2 + k3 + k3 + k4) / 6


def adams_bashford1(f, x0, y0, n):

    x = vx
    y = []
    nev = []
    h = 1. / n
    df_adams = DataFrame(columns=['x','y', 'y_exact', 'e(x)'])
    for i in range(4):
        y.append(vy[i])
        df_adams.loc[i] = [vx[i], y[i], sol(vx[i]), y[i] - sol(vx[i])]
        print("x = {:4}, y = {:20}, y_exact = {:20}, e(x) = {:25}".format(x[i], y[i], (x[i] - 1) ** 2 - x[i] + 2,
                                                                          y[i] - sol(x[i])))
        if i < 3:
            nev.append(y[i])
    for i in range(3, n):
        y.append(y[i] + h * (55 * f(x[i], y[i]) - 59 * f(x[i-1], y[i-1]) + 37 * f(x[i-2], y[i-2]) - 9 * f(x[i - 3], y[i - 3])) / 24)
        print("x = {:4}, y = {:20}, y_exact = {:20}, e(x) = {:25}".format(x[i], y[i], (x[i] - 1) ** 2 - x[i] + 2, y[i] - sol(x[i])))
        df_adams.loc[i] = [x[i], y[i], sol(x[i]), y[i] - sol(x[i])]
        nev.append(y[i] - sol(x[i]))
    nev.append(y[n] - sol(x[n]))

    print("x = {:4}, y = {:20}, y_exact = {:20}, e(x) = {:25}".format(x[n], y[n], (x[n] - 1) ** 2 - x[n] + 2,
                                                                      y[n] - sol(x[n])))
    df_adams.loc[n] = [x[n], y[n], sol(x[n]), y[n] - sol(x[n])]
    return y, nev, df_adams


errh1 = []
errh2 = []
# for h in [10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80]:
h = 200
vx, vy = rk4(f, 2, 1, 2+1, h)
err = []
df_rk = DataFrame(columns=['x','y', 'y_exact', 'e(x)'])
i = 0
for x, y in list(zip(vx, vy)):
    print("%4.2f %10.5f %+12.4e" % (x, y, y - sol(x)))
    df_rk.loc[i] = [x, y, sol(x), y - sol(x)]
    i += 1
    err.append(y - sol(x))

y, nev, df_adams = adams_bashford1(f, 2, 1, h)
errh1.append(max([abs(e) for e in err]))
errh2.append(max([abs(e) for e in nev]))
# df_adams.to_csv('Adams.csv')
# df_rk.to_csv('RungeKutta.csv')

from matplotlib import pyplot as plt
import seaborn

# plt.plot([10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80], errh1, label='RungeKutta')
# plt.legend(loc='best')
# plt.title('E(h)')
# plt.figure(2)
# plt.title('E(h)')
# plt.plot([10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80], errh2, label='Adams')
# plt.legend(loc='best')
# plt.show()


plt.plot(vx, vy, label='RungeKutta')
plt.plot(vx, [sol(x) for x in vx], label='Exact')
plt.plot(vx, y, label='Adams')
plt.legend(loc='best')
plt.title('Functions')

plt.figure(2)
plt.plot(vx, err, label='RungeKutta')
plt.plot(vx, nev, label='Adams')
plt.legend(loc='best')
plt.title('E(x)')

plt.show()

#
# err = []
# h = []
# for i in range(1, 100):
#     err.append(adams_bashford(f, 2, 1, i*0.005))
#     h.append(i*0.005)
#
# print(err, h)
#
# import matplotlib.pyplot as plt
# plt.plot(h, err)
# plt.show()