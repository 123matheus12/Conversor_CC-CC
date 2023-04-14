from matplotlib import pyplot as plt
from scipy.integrate import odeint
import numpy as np


def f(u, x):
    return u[1], - 2 * u[1] - 2 * u[0] + np.cos(2 * x)


y0 = [0, 0]
xs = np.linspace(1, 10, 200)
us = odeint(f, y0, xs)
ys = us[:, 0]

plt.plot(xs, ys, '-')
plt.plot(xs, ys, 'r*')
plt.xlabel('x values')
plt.ylabel('y values')
plt.title('(D**2-2*D+2)y = cos(2*x)')
plt.show()
