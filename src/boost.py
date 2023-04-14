from __future__ import unicode_literals  # accents texts in matplotlib
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
import time


def model(x, t, u1, d1, R, L, C):
    il1 = x[0]
    vc1 = x[1]

    dil1dt = 0 * il1 + (-(1 - d1) / L) * vc1 + (1 / L) * Vg
    dvc1dt = ((1 - d1) / C) * il1 + (-1 / (C * R)) * vc1 + 0 * Vg

    dxdt = [dil1dt, dvc1dt]

    return dxdt


def Solve_Dif_equations():
    # initial condition
    x0 = [0, 0]

    # inputs
    u1 = Vg * np.ones(len(t))

    # Razão cíclica constante
    d1 = D * np.ones(len(t))

    # # Razão cíclica com degrau
    # D1 = 0.8
    # D2 = 0.5
    # d_temp = []
    # for n in t:
    #     if n < (t[-1] / 2):
    #         d_temp.append(D1)
    #     else:
    #         d_temp.append(D2)
    # d1 = d_temp

    ## Razão cíclica variável
    # f = 60.0
    # w = 2*np.pi*f
    # Vp = Vg/2
    # Vc1 = Vp + Vp*np.sin(w*t)
    # d1 = Vc1/Vg

    # store solution
    il1 = np.empty_like(t)
    vc1 = np.empty_like(t)

    # record initial conditions
    il1[0] = x0[0]
    vc1[0] = x0[1]

    # solve ODE
    for i in range(1, len(t)):
        # span for next time step
        tspan = [t[i - 1], t[i]]
        # solve for next step
        x = odeint(model, x0, tspan, args=(u1[i], d1[i], R, L, C))
        # store solution for plotting
        il1[i] = x[1][0]
        vc1[i] = x[1][1]
        # next initial condition
        x0 = x[1]

    plt.figure()

    plt.subplot(211)  # create window plot with 2 rows and 2 columns
    # plt.subplots_adjust(hspace=0.5)
    # plt.plot(t, il1, 'r', label = '${i_{L}}_{Conversor}$')
    plt.plot(t, il1, 'r')
    plt.title('Corrente no Indutor')
    # plt.legend(loc='best')
    # plt.xlabel('t (ms)')
    plt.ylabel('i (A)')
    plt.tight_layout()
    plt.grid(True)

    plt.subplot(212)
    plt.plot(t, vc1, 'b')
    # plt.legend(loc='best')
    plt.title('Tensão no Capacitor')
    plt.xlabel('t (ms)')
    plt.ylabel('V (V)')
    plt.tight_layout()
    plt.grid(True)


Vg = 12.0
fs = 20000.0

L = 850.0e-6
C = 100.0e-6
R = 18.0
D = 0.5

stop_time = 0.05
t = np.arange(0, stop_time - 1 / fs, 1.0 / fs)
Solve_Dif_equations()

plt.show()
