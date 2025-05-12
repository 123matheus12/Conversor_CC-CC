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

        result = result / len(List)

        return result

    def calc_rms(self, List):

        result = 0.0

        List = List[0:int(len(List) / 2)]

        for n in (List):
            result = result + n ** 2

        result = result / len(List)
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
        vl1 = np.empty_like(t)
        vc1 = np.empty_like(t)
        ic1 = np.empty_like(t)
        is1 = np.empty_like(t)
        vs1 = np.empty_like(t)
        vd1 = np.empty_like(t)
        id1 = np.empty_like(t)
        vr1 = np.empty_like(t)
        ir1 = np.empty_like(t)
        vg1 = np.empty_like(t)
        ig1 = np.empty_like(t)
        pwm_vec = self.pwm(t, P, D)

        # record initial conditions
        il1[0] = x0[0]
        vc1[0] = x0[1]
        ir1[0] = x0[1] / R
        # solve ODE
        for i in range(1, len(t)):
            # span for next time step
            tspan = [t[i - 1], t[i]]
            # solve for next step
            if pwm_vec[i] == 1:
                x = odeint(self.buck_model_closed, x0, tspan, args=(u1[i], Vg, R, Rl, Rds, L, C))
                is1[i] = x[1][0]
                vs1[i] = 0
                vl1[i] = Vg - x[1][1]
                vd1[i] = -Vg
                id1[i] = 0
                ig1[i] = x[1][0]
            else:
                x = odeint(self.buck_model_open, x0, tspan, args=(u1[i], Vg, Vd, R, Rl, Rd, L, C))
                is1[i] = 0
                vs1[i] = -(Vg - Vd)
                vl1[i] = x[1][1] - Vd
                vd1[i] = Vd
                id1[i] = x[1][0]
                ig1[i] = 0
            # store solution for plotting
            il1[i] = x[1][0]
            vc1[i] = x[1][1]
            vr1[i] = x[1][1]
            ir1[i] = x[1][1] / R
            ic1[i] = x[1][0] - ir1[i]
            vg1[i] = Vg

            self.v_vector = vc1
            # next initial condition
            x0 = x[1]

        plt.figure()

        plt.subplot(321)  # create window plot with 4 rows and 1 columns
        plt.subplots_adjust(hspace=0.5)
        plt.plot(t, ig1, 'r', label='I Fonte')
        plt.plot(t, vg1, 'b', label='V Fonte')
        plt.title('Fonte')
        plt.xlabel('t (s)')
        plt.grid(True)

        plt.legend()

        plt.subplot(322)  # create window plot with 4 rows and 1 columns
        plt.subplots_adjust(hspace=0.5)
        plt.plot(t, il1, 'r', label='I Indutor')
        plt.plot(t, vl1, 'b', label='V Indutor')
        plt.title('Indutor')
        plt.xlabel('t (s)')
        plt.grid(True)

        plt.legend()

        plt.subplot(323)
        plt.plot(t, vc1, 'b', label='V Capacitor')
        plt.plot(t, ic1, 'r', label='I Capacitor')
        plt.title('Capacitor')
        plt.xlabel('t (s)')
        plt.grid(True)

        plt.legend()

        plt.subplot(324)
        plt.plot(t, vd1, 'b', label='V Diodo')
        plt.plot(t, id1, 'r', label='I Diodo')
        plt.title('Diodo')
        plt.xlabel('t (s)')
        plt.grid(True)

        plt.legend()

        plt.subplot(325)
        plt.plot(t, vr1, 'b', label='V Carga')
        plt.plot(t, ir1, 'r', label='I Carga')
        plt.title('Carga')
        plt.xlabel('t (s)')
        plt.grid(True)

        plt.legend()

        plt.subplot(326)
        plt.plot(t, vs1, 'b', label='V Chave')
        plt.plot(t, is1, 'r', label='I Chave')
        plt.title('Chave')
        plt.xlabel('t (s)')
        plt.grid(True)

        plt.legend()

    def Solve_Dif_equations_boost(self, t, Vg, Vd, P, D, R, Rl, Rds, Rd, L, C):
        # initial condition
        x0 = [0, 0]

        # inputs
        u1 = Vg * np.ones(len(t))

        # store solution
        il1 = np.empty_like(t)
        vl1 = np.empty_like(t)
        vc1 = np.empty_like(t)
        ic1 = np.empty_like(t)
        is1 = np.empty_like(t)
        vs1 = np.empty_like(t)
        vd1 = np.empty_like(t)
        id1 = np.empty_like(t)
        vr1 = np.empty_like(t)
        ir1 = np.empty_like(t)
        vg1 = np.empty_like(t)
        ig1 = np.empty_like(t)
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
                vl1[i] = Vg
                ic1[i] = x[1][1] / R
                vs1[i] = 0
                is1[i] = x[1][0]
                vd1[i] = - x[1][0]
                id1[i] = 0
            else:
                x = odeint(self.boost_model_open, x0, tspan, args=(u1[i], Vg, Vd, R, Rl, Rd, L, C))
                vl1[i] = Vg - Vd - x[1][1]
                ic1[i] = x[1][0] - x[1][1] / R
                vs1[i] = x[1][1]
                is1[i] = 0
                vd1[i] = 0
                id1[i] = x[1][0]
            # store solution for plotting
            il1[i] = x[1][0]
            vc1[i] = x[1][1]
            vr1[i] = x[1][1]
            ir1[i] = x[1][1] / R
            ig1[i] = x[1][0]
            vg1[i] = Vg
            self.v_vector = vc1
            # next initial condition
            x0 = x[1]

        plt.figure()

        plt.subplot(321)  # create window plot with 4 rows and 1 columns
        plt.subplots_adjust(hspace=0.5)
        plt.plot(t, ig1, 'r', label='I Fonte')
        plt.plot(t, vg1, 'b', label='V Fonte')
        plt.title('Fonte')
        plt.xlabel('t (s)')
        plt.grid(True)

        plt.legend()

        plt.subplot(322)  # create window plot with 4 rows and 1 columns
        plt.subplots_adjust(hspace=0.5)
        plt.plot(t, il1, 'r', label='I Indutor')
        plt.plot(t, vl1, 'b', label='V Indutor')
        plt.title('Indutor')
        plt.xlabel('t (s)')
        plt.grid(True)

        plt.legend()

        plt.subplot(323)
        plt.plot(t, vc1, 'b', label='V Capacitor')
        plt.plot(t, ic1, 'r', label='I Capacitor')
        plt.title('Capacitor')
        plt.xlabel('t (s)')
        plt.grid(True)

        plt.legend()

        plt.subplot(324)
        plt.plot(t, vd1, 'b', label='V Diodo')
        plt.plot(t, id1, 'r', label='I Diodo')
        plt.title('Diodo')
        plt.xlabel('t (s)')
        plt.grid(True)

        plt.legend()

        plt.subplot(325)
        plt.plot(t, vr1, 'b', label='V Carga')
        plt.plot(t, ir1, 'r', label='I Carga')
        plt.title('Carga')
        plt.xlabel('t (s)')
        plt.grid(True)

        plt.legend()

        plt.subplot(326)
        plt.plot(t, vs1, 'b', label='V Chave')
        plt.plot(t, is1, 'r', label='I Chave')
        plt.title('Chave')
        plt.xlabel('t (s)')
        plt.grid(True)

        plt.legend()

    def Solve_Dif_equations_buck_boost(self, t, Vg, Vd, P, D, R, Rl, Rds, Rd, L, C):
        # initial condition
        x0 = [0, 0]

        # inputs
        u1 = Vg * np.ones(len(t))

        # store solution
        il1 = np.empty_like(t)
        vl1 = np.empty_like(t)
        vc1 = np.empty_like(t)
        ic1 = np.empty_like(t)
        is1 = np.empty_like(t)
        vs1 = np.empty_like(t)
        vd1 = np.empty_like(t)
        id1 = np.empty_like(t)
        vr1 = np.empty_like(t)
        ir1 = np.empty_like(t)
        vg1 = np.empty_like(t)
        ig1 = np.empty_like(t)
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
            vl1[i] = Rd + x[1][1]
            vc1[i] = x[1][1]
            vr1[i] = x[1][1]
            ir1[i] = x[1][1] / R
            self.v_vector = vc1
            # next initial condition
            x0 = x[1]

        plt.figure()

        plt.subplot(321)  # create window plot with 4 rows and 1 columns
        plt.subplots_adjust(hspace=0.5)
        plt.plot(t, ig1, 'r', label='I Fonte')
        plt.plot(t, vg1, 'b', label='V Fonte')
        plt.title('Fonte')
        plt.xlabel('t (s)')
        plt.grid(True)

        plt.legend()

        plt.subplot(322)  # create window plot with 4 rows and 1 columns
        plt.subplots_adjust(hspace=0.5)
        plt.plot(t, il1, 'r', label='I Indutor')
        plt.plot(t, vl1, 'b', label='V Indutor')
        plt.title('Indutor')
        plt.xlabel('t (s)')
        plt.grid(True)

        plt.legend()

        plt.subplot(323)
        plt.plot(t, vc1, 'b', label='V Capacitor')
        plt.plot(t, ic1, 'r', label='I Capacitor')
        plt.title('Capacitor')
        plt.xlabel('t (s)')
        plt.grid(True)

        plt.legend()

        plt.subplot(324)
        plt.plot(t, vd1, 'b', label='V Diodo')
        plt.plot(t, id1, 'r', label='I Diodo')
        plt.title('Diodo')
        plt.xlabel('t (s)')
        plt.grid(True)

        plt.legend()

        plt.subplot(325)
        plt.plot(t, vr1, 'b', label='V Carga')
        plt.plot(t, ir1, 'r', label='I Carga')
        plt.title('Carga')
        plt.xlabel('t (s)')
        plt.grid(True)

        plt.legend()

        plt.subplot(326)
        plt.plot(t, vs1, 'b', label='V Chave')
        plt.plot(t, is1, 'r', label='I Chave')
        plt.title('Chave')
        plt.xlabel('t (s)')
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

        return [self.calc_average(vector), self.calc_rms(vector)]


if __name__ == '__main__':
    pass
