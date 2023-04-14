import numpy as np


def eq(x):
    return 1 - x / np.pi + np.sin(2 * x) / (2 * np.pi)


vec = np.linspace(0., 3.14, 10000)

esp = np.linspace(1., 100., 100)

y = eq(vec)

res = [(np.pi, 180, 0)]

for n in esp:
    i = 1
    while i < 10000:
        # print(y[i] / n)
        if 0.0099 < (y[i] / n) < 0.00999:
            res.append((vec[i], vec[i] * 180 / np.pi, y[i]))
            break
        i += 1

res.append((0, 0, 100))

z = 0

for n in res:
    print(f'[{z}],{n[1]:.2f},{n[0]:.5f}')
    z += 1
