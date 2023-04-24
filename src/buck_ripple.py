from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt


def pwm(t, P, D):
    pwm_vec = np.empty_like(t)

    p = D * P
    Cycle = 0

    for i in range(0, len(pwm_vec)):

        if Cycle == P:
            Cycle = 0

        if Cycle < p:
            pwm_vec[i] = 1
        else:
            pwm_vec[i] = 0

        Cycle += 1

    return pwm_vec


def model_open(x, t, u1, d1, R, L, C):
    il1 = x[0]
    vc1 = x[1]

    dil1dt = 0 * il1 + (-1 / L) * vc1 + 0 * Vg
    dvc1dt = (1 / C) * il1 + (-1 / (C * R)) * vc1

    dxdt = [dil1dt, dvc1dt]

    return dxdt


def model_closed(x, t, u1, d1, R, L, C):
    il1 = x[0]
    vc1 = x[1]

    dil1dt = 0 * il1 + (-1 / L) * vc1 + (d1 / L) * Vg
    dvc1dt = (1 / C) * il1 + (-1 / (C * R)) * vc1

    dxdt = [dil1dt, dvc1dt]

    return dxdt


def Solve_Dif_equations():
    # initial condition
    x0 = [0, 0]

    # inputs
    u1 = Vg * np.ones(len(t))

    # Razão cíclica constante
    d1 = D * np.ones(len(t))

    # store solution
    il1 = np.empty_like(t)
    vc1 = np.empty_like(t)
    pwm_vec = pwm(t, P, D)

    # record initial conditions
    il1[0] = x0[0]
    vc1[0] = x0[1]
    # solve ODE
    for i in range(1, len(t)):
        # span for next time step
        tspan = [t[i - 1], t[i]]
        # solve for next step
        if pwm_vec[i] == 1:
            x = odeint(model_closed, x0, tspan, args=(u1[i], d1[i], R, L, C))
        else:
            x = odeint(model_open, x0, tspan, args=(u1[i], d1[i], R, L, C))
        # store solution for plotting
        il1[i] = x[1][0]
        vc1[i] = x[1][1]
        # next initial condition
        x0 = x[1]

    plt.figure()

    plt.subplot(211)  # create window plot with 2 rows and 2 columns
    # plt.subplots_adjust(hspace=0.5)
    # plt.plot(t, il1, 'r', label = '${i_{L}}_{Conversor}$')
    plt.plot(t, pwm_vec, 'r')
    plt.title('PWM')
    plt.ylabel('PWM')
    plt.tight_layout()
    plt.grid(True)

    plt.subplot(212)
    plt.plot(t, vc1, 'b')
    # plt.legend(loc='best')
    plt.title('Tensão no Capacitor')
    plt.xlabel('t (s)')
    plt.ylabel('V (V)')
    plt.tight_layout()
    plt.grid(True)


Vg = 12.0
fs = 10000.0
P = 100  # quantidade de pontos por período de chaveamento

L = 729.0e-6
C = 20.0e-6
R = 2.5
D = 0.5

stop_time = 0.004
t = np.arange(0, (stop_time - 1 / (fs * P)), 1.0 / (fs * P))
Solve_Dif_equations()

plt.show()
