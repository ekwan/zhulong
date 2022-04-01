import pandas as pd
import numpy as np
import pickle

import sys
sys.path.append('../zhulong/')

from optimizer import *
from optimizer_utils import *

with open('sim_model.pkl', 'rb') as f:
    cols_subset = pickle.load(f)
    reg = pickle.load(f)

def simulation_fun(X_new, cols_X_new, cols_subset, reg, sigma = 0):
    cols_subset_index = [i for i in range(len(cols_X_new)) if cols_X_new[i] in cols_subset]
    pred_z = reg.predict(X_new[:,cols_subset_index])
    pred_z = np.squeeze(pred_z)
    z_sim = 0.5*pred_z + \
        ((1-(X_new[:,7]-4) * (X_new[:,7]-4)*0.2) + \
        0.4*(1-X_new[:,1])*X_new[:,2] + \
        (-1)*X_new[:,5]*X_new[:,6]*0.25 + \
        (-1)*X_new[:,5]*X_new[:,7]*0.3)*0.2
    if sigma > 0:
        print(z_sim.shape)
        z_sim = z_sim + np.random.normal(0, sigma, z_sim.shape[0])[:,np.newaxis]
    f_sim = transform_z2f(z_sim)
    return f_sim

with open('parameter_bounds.pkl', 'rb') as f:
    parameter_bounds = pickle.load(f)

method = 'LM'# or'RS'
initMethod = 'batch8' # 'random' or 'batch8' for LM
n_total = 16

dat = pd.DataFrame()
f_best = 0
for r in range(n_total):
    new_para = optimization(dat, parameter_bounds,method, seed = 999, initMethod = initMethod)
    df_new_para = pd.DataFrame({k: [v] for k, v in new_para.items()})
    X_new, cols_X_new = preprocessing_lm(df_new_para,X_only = True)
    f_sim = simulation_fun(X_new, cols_X_new, cols_subset, reg)
    f_best = max(f_best, f_sim)
    df_new_para = df_new_para.assign(f = f_sim, f_best = f_best)
    dat = pd.concat([dat,df_new_para],ignore_index = True)
    print("Round: ", r+1, "f: ", f_sim,"f_best",f_best,"\n")

dat = dat.assign(Round = range(1,n_total+1))
dat = dat.assign(method = method)
dat.to_csv("sim_"+method+".csv")
