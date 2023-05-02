import PySimpleGUI as sg
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# First the window layout in 2 columns
file_list_column = [
    [
        sg.Text("Indutor       "),
        sg.In(default_text="729", size=(25, 1), enable_events=True, key="-IND-"),
        sg.Text("uH"),
    ],
    [
        sg.Text("Capacitor   "),
        sg.In(default_text="20", size=(25, 1), enable_events=True, key="-CAP-"),
        sg.Text("uC"),
    ],
    [
        sg.Text("Resistência"),
        sg.In(default_text="2.5", size=(25, 1), enable_events=True, key="-RES-"),
        sg.Text("Ohm"),
    ],
    [
        sg.Text("Frequência "),
        sg.In(default_text="10000", size=(25, 1), enable_events=True, key="-F-"),
        sg.Text("Hz"),
    ],
    [
        sg.Text("Passo        "),
        sg.In(default_text="100", size=(25, 1), enable_events=True, key="-P-"),
    ],
    [
        sg.Text("Duty Clycle"),
        sg.In(default_text="0.413", size=(25, 1), enable_events=True, key="-D-"),
    ],
    [
        sg.Text("Tensão      "),
        sg.In(default_text="12", size=(25, 1), enable_events=True, key="-Vg-"),
        sg.Text("Volts"),
    ],
    [
        sg.Button("Simulate"),
    ],
]

# ----- Full layout -----
layout = [
    [
        sg.Column(file_list_column),
    ]
]

window = sg.Window("Conversor CC-CC", layout)


def pwm(t, P, D):  # gerar um vetor pwm, com P pontos por período na razão cíclica de D.
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


def model_open(x, t, u1, Vg, R, L, C):  # modelo para chave aberta
    il1 = x[0]
    vc1 = x[1]

    dil1dt = 0 * il1 + (-1 / L) * vc1 + 0 * Vg
    dvc1dt = (1 / C) * il1 + (-1 / (C * R)) * vc1

    dxdt = [dil1dt, dvc1dt]

    return dxdt


def model_closed(x, t, u1, Vg, R, L, C):  # modelo para chave fechada
    il1 = x[0]
    vc1 = x[1]

    dil1dt = 0 * il1 + (-1 / L) * vc1 + (1 / L) * Vg
    dvc1dt = (1 / C) * il1 + (-1 / (C * R)) * vc1

    dxdt = [dil1dt, dvc1dt]

    return dxdt


def Solve_Dif_equations(t, Vg, P, D, R, L, C):
    # initial condition
    x0 = [0, 0]

    # inputs
    u1 = Vg * np.ones(len(t))

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
            x = odeint(model_closed, x0, tspan, args=(u1[i], Vg, R, L, C))
        else:
            x = odeint(model_open, x0, tspan, args=(u1[i], Vg, R, L, C))
        # store solution for plotting
        il1[i] = x[1][0]
        vc1[i] = x[1][1]
        # next initial condition
        x0 = x[1]

    plt.figure()

    plt.subplot(311)  # create window plot with 2 rows and 2 columns
    plt.subplots_adjust(hspace=0.5)
    plt.plot(t, il1, 'r', label='${i_{L}}_{Conversor}$')
    plt.title('Corrente no Indutor')
    plt.xlabel('t (s)')
    plt.ylabel('I (A)')
    plt.grid(True)

    plt.subplot(312)
    plt.plot(t, vc1, 'b')
    # plt.legend(loc='best')
    plt.title('Tensão no Capacitor')
    plt.xlabel('t (s)')
    plt.ylabel('V (V)')
    plt.grid(True)

    plt.subplot(313)
    plt.plot(t, pwm_vec, 'g')
    plt.title('PWM')
    plt.xlabel('t (s)')
    plt.grid(True)


# Run the Event Loop
while True:
    event, values = window.read()
    if event == "Exit" or event == sg.WIN_CLOSED:
        break
    # Folder name was filled in, make a list of files in the folder
    if event == "Simulate":
        try:
            Vg = float(values["-Vg-"])
            fs = float(values["-F-"])
            P = float(values["-P-"])
            L = float(values["-IND-"]) * 10 ** (-6)
            C = float(values["-CAP-"]) * 10 ** (-6)
            R = float(values["-RES-"])
            D = float(values["-D-"])

            stop_time = 0.004

            t = np.arange(0, (stop_time - 1 / (fs * P)), 1.0 / (fs * P))

            Solve_Dif_equations(t, Vg, P, D, R, L, C)
            plt.show()
        except:
            pass

window.close()
