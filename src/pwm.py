import numpy as np
import matplotlib.pyplot as plt


def pwm(t, D):
    pwm_vec = np.empty_like(t)

    p = D * 100
    Cycle = 0

    for i in range(0, len(pwm_vec)):

        if Cycle == 100:
            Cycle = 0

        if Cycle < p:
            pwm_vec[i] = 1
        else:
            pwm_vec[i] = 0

        Cycle += 1

    return pwm_vec

