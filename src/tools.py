import PySimpleGUI as sg
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt


class Solver:
    def __init__(self):
        self.v_vector = list()

    def calc_average(self, List):
        
        result = 0.0
        
        for n in (List):
            result = result + n
        
        result = result/len(List)
            
        return result


    def calc_rms(self, List):
            
        result = 0.0
        
        List = List[0:int(len(List)/2)]
        
        for n in (List):
            result = result + n**2
        
        result = result/len(List)
        result = np.sqrt(result)
        
        return result


    def pwm(self, t, P, D):  # gerar um vetor pwm, com P pontos por período na razão cíclica de D.
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


    def buck_model_open(self, x, t, u1, Vg, Vd, R, Rl, Rd, L, C):  # modelo para chave aberta
        il1 = x[0]
        vc1 = x[1]

        dil1dt = (-(Rl + Rd) / L) * il1 + (-1 / L) * vc1 + (-1 / L) * Vd + 0 * Vg
        dvc1dt = (1 / C) * il1 + (-1 / (C * R)) * vc1 + 0 * Vg

        dxdt = [dil1dt, dvc1dt]
        return dxdt


    def buck_model_closed(self, x, t, u1, Vg, R, Rl, Rds, L, C):  # modelo para chave fechada
        il1 = x[0]
        vc1 = x[1]

        dil1dt = (-(Rl + Rds) / L) * il1 + (-1 / L) * vc1 + (1 / L) * Vg
        dvc1dt = (1 / C) * il1 + (-1 / (C * R)) * vc1 + 0 * Vg

        dxdt = [dil1dt, dvc1dt]

        return dxdt
    
    def boost_model_open(self, x, t, u1, Vg, Vd, R, Rl, Rd, L, C):  # modelo para chave aberta
        il1 = x[0]
        vc1 = x[1]

        dil1dt = (-(Rl + Rd) / L) * il1 + (-1 / L) * vc1 + (-1 / L) * Vd + (1 / L) * Vg
        dvc1dt = (1 / C) * il1 + (-1 / (C * R)) * vc1 + 0 * Vg

        dxdt = [dil1dt, dvc1dt]
        return dxdt


    def boost_model_closed(self, x, t, u1, Vg, R, Rl, Rds, L, C):  # modelo para chave fechada
        il1 = x[0]
        vc1 = x[1]

        dil1dt = (-(Rl + Rds) / L) * il1 + 0 * vc1 + (1 / L) * Vg
        dvc1dt = 0 * il1 + (-1 / (C * R)) * vc1 + 0 * Vg

        dxdt = [dil1dt, dvc1dt]

        return dxdt
    
    def buck_boost_model_open(self, x, t, u1, Vg, Vd, R, Rl, Rd, L, C):  # modelo para chave aberta
        il1 = x[0]
        vc1 = x[1]

        dil1dt = (-(Rl + Rd) / L) * il1 + (-1 / L) * vc1 + (1 / L) * Vd + 0 * Vg
        dvc1dt = (1 / C) * il1 + (-1 / (C * R)) * vc1 + 0 * Vg

        dxdt = [dil1dt, dvc1dt]
        return dxdt


    def buck_boost_model_closed(self, x, t, u1, Vg, R, Rl, Rds, L, C):  # modelo para chave fechada
        il1 = x[0]
        vc1 = x[1]

        dil1dt = (-(Rl + Rds) / L) * il1 + 0 * vc1 + (1 / L) * Vg
        dvc1dt = 0 * il1 + (-1 / (C * R)) * vc1 + 0 * Vg

        dxdt = [dil1dt, dvc1dt]

        return dxdt


    def Solve_Dif_equations_buck(self, t, Vg, Vd, P, D, R, Rl, Rds, Rd, L, C):
        # initial condition
        x0 = [0, 0]

        # inputs
        u1 = Vg * np.ones(len(t))

        # store solution
        il1 = np.empty_like(t)
        vc1 = np.empty_like(t)
        pwm_vec = self.pwm(t, P, D)

        # record initial conditions
        il1[0] = x0[0]
        vc1[0] = x0[1]
        # solve ODE
        for i in range(1, len(t)):
            # span for next time step
            tspan = [t[i - 1], t[i]]
            # solve for next step
            if pwm_vec[i] == 1:
                x = odeint(self.buck_model_closed, x0, tspan, args=(u1[i], Vg, R, Rl, Rds, L, C))
            else:
                x = odeint(self.buck_model_open, x0, tspan, args=(u1[i], Vg, Vd, R, Rl, Rd, L, C))
            # store solution for plotting
            il1[i] = x[1][0]
            vc1[i] = x[1][1]
            self.v_vector = vc1
            # next initial condition
            x0 = x[1]

        plt.figure()

        plt.subplot(211)  # create window plot with 2 rows and 2 columns
        plt.subplots_adjust(hspace=0.5)
        plt.plot(t, il1, 'r', label='Indutor')
        plt.title('Corrente no Indutor')
        plt.xlabel('t (s)')
        plt.ylabel('I (A)')
        plt.grid(True)

        plt.legend()

        plt.subplot(212)
        plt.plot(t, vc1, 'b', label='Capacitor')
        plt.title('Tensão no Capacitor')
        plt.xlabel('t (s)')
        plt.ylabel('V (V)')
        plt.grid(True)

        plt.legend()

    def Solve_Dif_equations_boost(self, t, Vg, Vd, P, D, R, Rl, Rds, Rd, L, C):
        # initial condition
        x0 = [0, 0]

        # inputs
        u1 = Vg * np.ones(len(t))

        # store solution
        il1 = np.empty_like(t)
        vc1 = np.empty_like(t)
        pwm_vec = self.pwm(t, P, D)

        # record initial conditions
        il1[0] = x0[0]
        vc1[0] = x0[1]
        # solve ODE
        for i in range(1, len(t)):
            # span for next time step
            tspan = [t[i - 1], t[i]]
            # solve for next step
            if pwm_vec[i] == 1:
                x = odeint(self.boost_model_closed, x0, tspan, args=(u1[i], Vg, R, Rl, Rds, L, C))
            else:
                x = odeint(self.boost_model_open, x0, tspan, args=(u1[i], Vg, Vd, R, Rl, Rd, L, C))
            # store solution for plotting
            il1[i] = x[1][0]
            vc1[i] = x[1][1]
            self.v_vector = vc1
            # next initial condition
            x0 = x[1]

        plt.figure()

        plt.subplot(211)  # create window plot with 2 rows and 2 columns
        plt.subplots_adjust(hspace=0.5)
        plt.plot(t, il1, 'r', label='Indutor')
        plt.title('Corrente no Indutor')
        plt.xlabel('t (s)')
        plt.ylabel('I (A)')
        plt.grid(True)

        plt.legend()

        plt.subplot(212)
        plt.plot(t, vc1, 'b', label='Capacitor')
        plt.title('Tensão no Capacitor')
        plt.xlabel('t (s)')
        plt.ylabel('V (V)')
        plt.grid(True)

        plt.legend()
    
    def Solve_Dif_equations_buck_boost(self, t, Vg, Vd, P, D, R, Rl, Rds, Rd, L, C):
        # initial condition
        x0 = [0, 0]

        # inputs
        u1 = Vg * np.ones(len(t))

        # store solution
        il1 = np.empty_like(t)
        vc1 = np.empty_like(t)
        pwm_vec = self.pwm(t, P, D)

        # record initial conditions
        il1[0] = x0[0]
        vc1[0] = x0[1]
        # solve ODE
        for i in range(1, len(t)):
            # span for next time step
            tspan = [t[i - 1], t[i]]
            # solve for next step
            if pwm_vec[i] == 1:
                x = odeint(self.buck_boost_model_closed, x0, tspan, args=(u1[i], Vg, R, Rl, Rds, L, C))
            else:
                x = odeint(self.buck_boost_model_open, x0, tspan, args=(u1[i], Vg, Vd, R, Rl, Rd, L, C))
            # store solution for plotting
            il1[i] = x[1][0]
            vc1[i] = x[1][1]
            self.v_vector = vc1
            # next initial condition
            x0 = x[1]

        plt.figure()

        plt.subplot(211)  # create window plot with 2 rows and 2 columns
        plt.subplots_adjust(hspace=0.5)
        plt.plot(t, il1, 'r', label='Indutor')
        plt.title('Corrente no Indutor')
        plt.xlabel('t (s)')
        plt.ylabel('I (A)')
        plt.grid(True)

        plt.legend()

        plt.subplot(212)
        plt.plot(t, vc1, 'b', label='Capacitor')
        plt.title('Tensão no Capacitor')
        plt.xlabel('t (s)')
        plt.ylabel('V (V)')
        plt.grid(True)

        plt.legend()
    
    def calculate_mean_rms(self, time, setpoint_min, setpoint_max):
        vector = list()

        i = 0
        if setpoint_max > setpoint_min:
            for a in time:
                if a >= setpoint_min and a <= setpoint_max:
                    vector.append(self.v_vector[i])
                i += 1
            print(len(vector))
        
        return [self.calc_average(vector), self.calc_rms(vector)]
    

if __name__ == '__main__':
    pass
