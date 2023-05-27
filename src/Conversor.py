import PySimpleGUI as sg
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
import os.path
import tools


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
        [sg.Text("Rl:", size=(10, 1)), sg.In(default_text="0.0025", size=(25, 1), enable_events=True, key="-Rl-"),
        sg.Text("Ohm"), ],
        [sg.Text("Rds:", size=(10, 1)), sg.In(default_text="0.0025", size=(25, 1), enable_events=True, key="-Rds-"),
        sg.Text("Ohm"), ],
        [sg.Text("Rd:", size=(10, 1)), sg.In(default_text="0.0025", size=(25, 1), enable_events=True, key="-Rd-"),
        sg.Text("Ohm"), ],
        [sg.Text("Vd:", size=(10, 1)), sg.In(default_text="0.7", size=(25, 1), enable_events=True, key="-Vd-"),
        sg.Text("Volts"), ],
        [sg.Button("Simulate"), ],
        [sg.Text("Cálculo do valor RMS:"), ],
        [sg.Slider(range=(0, 0.004), resolution=0.0001, orientation='h', enable_events=True, key='-tempo_rms_min-'),
         sg.Slider(range=(0, 0.004), resolution=0.0001, orientation='h', enable_events=True, key='-tempo_rms_max-'), ],
        [sg.Text("Valor RMS na carga:"), sg.Text("0", key='-rms-'), sg.Text("Valor médio na carga:"), sg.Text("0", key='-mean-'), ],
    ]

    selection_layout = [
        [sg.Text("Selecione o tipo de conversor: "), sg.OptionMenu(values=('Buck', 'Boost', 'Buck-Boost'), default_value='Buck', k='-CONV-'), ],
        [sg.Image(key="-IMAGE-"), ],
        [sg.Text("Chave:", size=(8, 1)), sg.Checkbox('Tensão', default=True, k='-sw_v-'), sg.Checkbox('Corrente', default=True, k='-sw_i-'), ],
        [sg.Text("Fonte:", size=(8, 1)), sg.Checkbox('Tensão', default=True, k='-source_v-'), sg.Checkbox('Corrente', default=True, k='-source_i-'), ],
        [sg.Text("Diodo:", size=(8, 1)), sg.Checkbox('Tensão', default=True, k='-dio_v-'), sg.Checkbox('Corrente', default=True, k='-dio_i-'), ],
        [sg.Text("Indutor:", size=(8, 1)), sg.Checkbox('Tensão', default=True, k='-ind_v-'), sg.Checkbox('Corrente', default=True, k='-ind_i-'), ],
        [sg.Text("Capacitor:", size=(8, 1)), sg.Checkbox('Tensão', default=True, k='-cap_v-'), sg.Checkbox('Corrente', default=True, k='-cap_i-'), ],
        [sg.Text("Carga:", size=(8, 1)), sg.Checkbox('Tensão', default=True, k='-res_v-'), sg.Checkbox('Corrente', default=True, k='-res_i-'), ],

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


# Run the Event Loop
def main():
    window = make_window(sg.theme())
    Solver = tools.Solver()

    while True:
        event, values = window.read(timeout=100)
        if event == "Exit" or event == sg.WIN_CLOSED:
            break
        
        filename = os.path.dirname(os.path.abspath(__file__))
        try:
            filename = filename + "\\imagens\\" + values["-CONV-"] + "_286X117.png"
            window["-IMAGE-"].update(filename=filename)
        except:
            pass

        # Folder name was filled in, make a list of files in the folder
        if event == "Simulate":
            try:
                Vg = float(values["-Vg-"])
                Vd = float(values["-Vd-"])
                fs = float(values["-F-"])
                P = float(values["-P-"])
                L = float(values["-IND-"]) * 10 ** (-6)
                C = float(values["-CAP-"]) * 10 ** (-6)
                R = float(values["-RES-"])
                Rl = float(values["-Rl-"])
                Rds = float(values["-Rds-"])
                Rd = float(values["-Rd-"])
                D = float(values["-D-"])

                stop_time = 0.004

                t = np.arange(0, (stop_time - 1 / (fs * P)), 1.0 / (fs * P))

                if values["-CONV-"] == "Buck":
                    Solver.Solve_Dif_equations_buck(t, Vg, Vd, P, D, R, Rl, Rds, Rd, L, C)
                elif values["-CONV-"] == "Boost":
                    Solver.Solve_Dif_equations_boost(t, Vg, Vd, P, D, R, Rl, Rds, Rd, L, C)
                elif values["-CONV-"] == "Buck-Boost":
                    Solver.Solve_Dif_equations_buck_boost(t, Vg, Vd, P, D, R, Rl, Rds, Rd, L, C)
                plt.show()
            except:
                sg.popup("Valor Inválido!", keep_on_top=True)

        elif event == "-tempo_rms_min-" or event == "-tempo_rms_max-":
            try:
                P = float(values["-P-"])
                fs = float(values["-F-"])
                stop_time = 0.004
                t = np.arange(0, (stop_time - 1 / (fs * P)), 1.0 / (fs * P))

                mean_rms = Solver.calculate_mean_rms(t, float(values["-tempo_rms_min-"]), float(values["-tempo_rms_max-"]))
                if type(mean_rms == float):
                    window["-rms-"].update(str(mean_rms[1]))
                    window["-mean-"].update(str(mean_rms[0]))
            except:
                pass

    window.close()
    exit(0)


if __name__ == '__main__':
    sg.theme('green')
    main()
