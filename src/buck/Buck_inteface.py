import PySimpleGUI as sg
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
import os.path


def make_window(theme):
    sg.theme(theme)
    menu_def = [['&Application', ['&Exit']],
                ['&Help', ['&About']]]

    right_click_menu_def = [[], ['Exit']]

    # First the window layout in 2 columns
    input_layout = [
        [sg.Text("Indutor:", size=(10, 1)), sg.In(default_text="729", size=(25, 1), enable_events=True, key="-IND-"),
         sg.Text("uH"), ],
        [sg.Text("Capacitor:", size=(10, 1)), sg.In(default_text="20", size=(25, 1), enable_events=True, key="-CAP-"),
         sg.Text("uC"), ],
        [sg.Text("Resistência:", size=(10, 1)), sg.In(default_text="2.5", size=(25, 1), enable_events=True, key="-RES-"),
         sg.Text("Ohm"), ],
        [sg.Text("Frequência:", size=(10, 1)), sg.In(default_text="10000", size=(25, 1), enable_events=True, key="-F-"),
         sg.Text("Hz"), ],
        [sg.Text("Passo:", size=(10, 1)), sg.In(default_text="100", size=(25, 1), enable_events=True, key="-P-"), ],
        [sg.Text("Duty Clycle:", size=(10, 1)), sg.In(default_text="0.413", size=(25, 1), enable_events=True, key="-D-"), ],
        [sg.Text("Tensão:", size=(10, 1)), sg.In(default_text="12", size=(25, 1), enable_events=True, key="-Vg-"),
         sg.Text("Volts"), ],
        [sg.Button("Simulate"), ],
    ]
    selection_layout = [
        [sg.Text("Selecione o tipo de conversor: "), sg.OptionMenu(values=('Buck', 'Boost', 'Buck-Boost'), default_value='Buck', k='-CONV-'), ],
        [sg.Image(key="-IMAGE-"), ],
        [sg.Text("Chave:", size=(8, 1)), sg.Checkbox('Tensão', default=True, k='-SW_V-'), sg.Checkbox('Corrente', default=True, k='-SW_I-'), ],
        [sg.Text("Fonte:", size=(8, 1)), sg.Checkbox('Tensão', default=True, k='-SW_V-'), sg.Checkbox('Corrente', default=True, k='-SW_I-'), ],
        [sg.Text("Diodo:", size=(8, 1)), sg.Checkbox('Tensão', default=True, k='-SW_V-'), sg.Checkbox('Corrente', default=True, k='-SW_I-'), ],
        [sg.Text("Indutor:", size=(8, 1)), sg.Checkbox('Tensão', default=True, k='-SW_V-'), sg.Checkbox('Corrente', default=True, k='-SW_I-'), ],
        [sg.Text("Capacitor:", size=(8, 1)), sg.Checkbox('Tensão', default=True, k='-SW_V-'), sg.Checkbox('Corrente', default=True, k='-SW_I-'), ],
        [sg.Text("Carga:", size=(8, 1)), sg.Checkbox('Tensão', default=True, k='-SW_V-'), sg.Checkbox('Corrente', default=True, k='-SW_I-'), ],

    ]

    # ----- Full layout -----
    layout = [[sg.MenubarCustom(menu_def, key='-MENU-', font='Courier 15', tearoff=True)],
              [sg.Text('Conversores CC-CC', size=(38, 1), justification='center', font=("Helvetica", 16),
                       relief=sg.RELIEF_RIDGE, k='-TEXT HEADING-', enable_events=True)]]

    layout += [[sg.TabGroup([[sg.Tab('Seleção', selection_layout),
                              sg.Tab('Parâmetros', input_layout), ]], key='-TAB GROUP-', expand_x=True,
                            expand_y=True),
                ]]

    layout[-1].append(sg.Sizegrip())
    window = sg.Window('Conversores', layout, right_click_menu=right_click_menu_def,
                       right_click_menu_tearoff=True, grab_anywhere=True, resizable=True, margins=(0, 0),
                       use_custom_titlebar=True, finalize=True, keep_on_top=True)
    window.set_min_size(window.size)
    return window


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
def main():
    window = make_window(sg.theme())

    while True:
        event, values = window.read(timeout=100)
        if event == "Exit" or event == sg.WIN_CLOSED:
            break

        try:
            filename = os.path.join("/home/matheus/PycharmProjects/pythonProject/src/buck/imagens/", values["-CONV-"] + ".png")
            window["-IMAGE-"].update(filename=filename)
        except:
            pass

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
                sg.popup("Coloque valores válidos!", keep_on_top=True)

    window.close()
    exit(0)


if __name__ == '__main__':
    sg.theme('green')
    main()
