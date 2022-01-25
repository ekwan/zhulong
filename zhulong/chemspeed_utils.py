# Helper functions for chemspeed.py

import math

"""
# Example for plateau_function: 

timepoints_all = range(0,70,5) # 14 timepoints
values_all = [0.49, 13.57, 21.37, 26.58, 30.19, 32.82, 33.82, 34.56, 35.35, 35.15, 35.28, 35.37, 35.12, 34.84]

#plateau is at timepoint==45, which is timepoints_all[9]

history = {'t' : timepoints_all[:5], 'v': values_all[:5]}
plateau_function(history)

for i in range(len(timepoints_all)):
    history = {'t' : timepoints_all[:(i+1)], 'v': values_all[:(i+1)]}
    print("i = ", i, "time:", timepoints_all[i],"plateau: ", plateau_function(history))

"""

def plateau_function(history):
    u0 = 0.00734 # hyperparameter for plateau detection algorithm
    plateau_flag = False # plateau indicator to be return
    # To-Do: get the list of timepoints and experimental outcome values from history
    #timepoints, values = a_function_of_history(history)
    timepoints, values = history['t'],history['v']
    n = len(timepoints) # the algorithm requires at least two measurements 
    if n >= 2:
        diff_value_t = values[-1]-values[-2]
        sigma_t = (values[-1]+values[-2])*0.5*u0
        w = 1.96*math.sqrt(2)*sigma_t
        if abs(diff_value_t) <= w:
            plateau_flag = True
        else:
            plateau_flag = False
    return plateau_flag
