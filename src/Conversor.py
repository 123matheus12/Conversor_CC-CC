import tkinter as tk
from tkinter import ttk, messagebox
from tkinter import PhotoImage
from PIL import Image, ImageTk
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
import os.path
import tools


class ConverterApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Conversores CC-CC")
        self.root.geometry("900x600")
        self.solver = tools.Solver()

        self.setup_menu()
        self.setup_tabs()
        self.setup_selection_tab()
        self.setup_parameters_tab()
        self.setup_events()

    def setup_menu(self):
        menubar = tk.Menu(self.root)
        app_menu = tk.Menu(menubar, tearoff=0)
        app_menu.add_command(label="Exit", command=self.root.quit)
        menubar.add_cascade(label="Application", menu=app_menu)

        help_menu = tk.Menu(menubar, tearoff=0)
        help_menu.add_command(label="About", command=lambda: messagebox.showinfo("Sobre", "Conversores CC-CC"))
        menubar.add_cascade(label="Help", menu=help_menu)

        self.root.config(menu=menubar)

    def setup_tabs(self):
        self.notebook = ttk.Notebook(self.root)
        self.tab_selection = ttk.Frame(self.notebook)
        self.tab_parameters = ttk.Frame(self.notebook)
        self.notebook.add(self.tab_selection, text="Seleção")
        self.notebook.add(self.tab_parameters, text="Parâmetros")
        self.notebook.pack(fill='both', expand=True)

    def setup_selection_tab(self):
        self.image_label = tk.Label(self.tab_selection)
        self.image_label.grid(row=1, column=0, columnspan=3, pady=10)

        self.converter_var = tk.StringVar(value="Buck")
        tk.Label(self.tab_selection, text="Tipo de Conversor:").grid(row=0, column=0, sticky="w")
        tk.OptionMenu(self.tab_selection, self.converter_var, "Buck", "Boost", "Buck-Boost", command=self.update_image).grid(row=0, column=1)

        self.check_vars = {}
        elements = ['sw', 'source', 'dio', 'ind', 'cap', 'res']
        for i, label in enumerate(['Chave', 'Fonte', 'Diodo', 'Indutor', 'Capacitor', 'Carga']):
            tk.Label(self.tab_selection, text=label + ":").grid(row=i + 2, column=0, sticky="w")
            self.check_vars[f"{elements[i]}_v"] = tk.BooleanVar(value=True)
            self.check_vars[f"{elements[i]}_i"] = tk.BooleanVar(value=True)
            tk.Checkbutton(self.tab_selection, text="Tensão", variable=self.check_vars[f"{elements[i]}_v"]).grid(row=i + 2, column=1)
            tk.Checkbutton(self.tab_selection, text="Corrente", variable=self.check_vars[f"{elements[i]}_i"]).grid(row=i + 2, column=2)

    def setup_parameters_tab(self):
        labels = [
            ("Indutor (uH):", "729"), ("Capacitor (uC):", "20"), ("Resistência (Ohm):", "2.5"),
            ("Frequência (Hz):", "10000"), ("Passo:", "100"), ("Duty Cycle:", "0.413"),
            ("Tensão (V):", "12"), ("Rl (Ohm):", "0.0025"), ("Rds (Ohm):", "0.0025"),
            ("Rd (Ohm):", "0.0025"), ("Vd (V):", "0.7")
        ]
        self.entries = {}
        for i, (label, default) in enumerate(labels):
            tk.Label(self.tab_parameters, text=label).grid(row=i, column=0, sticky="e")
            ent = tk.Entry(self.tab_parameters)
            ent.insert(0, default)
            ent.grid(row=i, column=1)
            self.entries[label] = ent

        self.simulate_button = tk.Button(self.tab_parameters, text="Simular", command=self.simulate)
        self.simulate_button.grid(row=len(labels), column=0, columnspan=2, pady=10)

        tk.Label(self.tab_parameters, text="Cálculo do valor RMS:").grid(row=len(labels)+1, column=0, columnspan=2)
        self.rms_min = tk.DoubleVar(value=0)
        self.rms_max = tk.DoubleVar(value=0)
        tk.Scale(self.tab_parameters, from_=0, to=0.004, resolution=0.0001, orient='horizontal', variable=self.rms_min, label="Tempo Min", command=self.update_rms).grid(row=len(labels)+2, column=0, sticky="ew")
        tk.Scale(self.tab_parameters, from_=0, to=0.004, resolution=0.0001, orient='horizontal', variable=self.rms_max, label="Tempo Max", command=self.update_rms).grid(row=len(labels)+2, column=1, sticky="ew")

        self.rms_label = tk.Label(self.tab_parameters, text="Valor RMS: 0 | Valor médio: 0")
        self.rms_label.grid(row=len(labels)+3, column=0, columnspan=2)

    def setup_events(self):
        self.update_image()

    def update_image(self, *args):
        conv = self.converter_var.get()
        image_path = os.path.join(os.path.dirname(__file__), "imagens", f"{conv}_286X117.png")
        try:
            image = Image.open(image_path)
            image = image.resize((286, 117), Image.ANTIALIAS)
            self.imgtk = ImageTk.PhotoImage(image)
            self.image_label.configure(image=self.imgtk)
        except Exception:
            self.image_label.configure(image="")

    def simulate(self):
        try:
            Vg = float(self.entries["Tensão (V):"].get())
            Vd = float(self.entries["Vd (V):"].get())
            fs = float(self.entries["Frequência (Hz):"].get())
            P = float(self.entries["Passo:"].get())
            L = float(self.entries["Indutor (uH):"].get()) * 1e-6
            C = float(self.entries["Capacitor (uC):"].get()) * 1e-6
            R = float(self.entries["Resistência (Ohm):"].get())
            Rl = float(self.entries["Rl (Ohm):"].get())
            Rds = float(self.entries["Rds (Ohm):"].get())
            Rd = float(self.entries["Rd (Ohm):"].get())
            D = float(self.entries["Duty Cycle:"].get())

            stop_time = 0.004
            t = np.arange(0, (stop_time - 1 / (fs * P)), 1.0 / (fs * P))

            conv = self.converter_var.get()
            if conv == "Buck":
                self.solver.Solve_Dif_equations_buck(t, Vg, Vd, P, D, R, Rl, Rds, Rd, L, C)
            elif conv == "Boost":
                self.solver.Solve_Dif_equations_boost(t, Vg, Vd, P, D, R, Rl, Rds, Rd, L, C)
            elif conv == "Buck-Boost":
                self.solver.Solve_Dif_equations_buck_boost(t, Vg, Vd, P, D, R, Rl, Rds, Rd, L, C)
            plt.show()
        except Exception as e:
            messagebox.showerror("Erro", "Valor Inválido!")

    def update_rms(self, *args):
        try:
            P = float(self.entries["Passo:"].get())
            fs = float(self.entries["Frequência (Hz):"].get())
            stop_time = 0.004
            t = np.arange(0, (stop_time - 1 / (fs * P)), 1.0 / (fs * P))
            mean_rms = self.solver.calculate_mean_rms(t, self.rms_min.get(), self.rms_max.get())
            if isinstance(mean_rms, tuple) and len(mean_rms) == 2:
                self.rms_label.config(text=f"Valor RMS: {mean_rms[1]:.4f} | Valor médio: {mean_rms[0]:.4f}")
        except:
            pass


if __name__ == '__main__':
    root = tk.Tk()
    app = ConverterApp(root)
    root.mainloop()
