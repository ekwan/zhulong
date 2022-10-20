# Helper functions for chemspeed.py

import math

"""
# Example for plateau_function:

timepoints_all = range(0,70,5) # 14 timepoints
values_all = [0.49, 13.57, 21.37, 26.58, 30.19, 32.82, 33.82, 34.56, 35.35, 35.15, 35.28, 35.37, 35.12, 34.84]

#plateau is at timepoint==45, which is timepoints_all[9]

history = {'times' : timepoints_all[:5], 'values': values_all[:5]}
plateau_function_1(history)

for i in range(len(timepoints_all)):
    history = {'times' : timepoints_all[:(i+1)], 'values': values_all[:(i+1)]}
    print("i = ", i, "time:", timepoints_all[i],"plateau: ", plateau_function_1(history))

Add to zhulong.py:
from chemspeed_utils import plateau_function_1
"""


def plateau_function_1(history):
    print("Plateau Algorithm No.1")
    u0 = 0.00734  # hyperparameter for plateau detection algorithm
    max_n_samples = 12  # additional termination criterion
    plateau_flag = False  # plateau indicator to be return
    timepoints = history['times']
    values = history['values']
    n = len(timepoints)  # the algorithm requires at least two measurements
    if n >= max_n_samples:
        plateau_flag = True
    else:
        if n >= 2:
            diff_value_t = values[-1] - values[-2]
            sigma_t = (values[-1] + values[-2]) * 0.5 * u0
            w = 1.96 * math.sqrt(2) * sigma_t
            if abs(diff_value_t) <= w:
                plateau_flag = True
            else:
                plateau_flag = False
    return plateau_flag
