"""
Inputs:
    dat: A pandas data frame of available experimental results for modeling
    parameter_bounds: A dictionary contains the min/max/step-size of all experimental conditions to be explored, 
                      produced by parameter_space.get_bounds_dict(), where the 'parameter_space' is a ParameterSpace class object 
    method: Choose an optimization method from the list ['RS','LM','BO']
        'RS' = Random Search
        'LM' = Polynomial regression modeling
        'BO' = Bayesian optimization
Outputs:
    newPara: A dictionary contains the candidate experimental condition for next round

Example Usage: ../test/simulation.py
Example Usage from zhulong.py: ../test/opt_for_zhulong_DoNotRun.py

"""

import pandas as pd
import numpy as np
from numpy.linalg import inv
import math
from sklearn.linear_model import LinearRegression
from scipy.stats import norm, pearsonr
from optimizer_utils import *
import random
import os

def initial_oneSample(parameter_bounds):
    dict_listAll = {
        'Reagent': parameter_bounds['reagents']['names'],
        'Reagent_equiv': np.arange(parameter_bounds['reagents']['min_equivalents'],
                                    parameter_bounds['reagents']['max_equivalents']+parameter_bounds['reagents']['step_size']/10.0,
                                    step=parameter_bounds['reagents']['step_size']).tolist(),
        'Solvent': parameter_bounds['solvents']['names'], 
        'Additive': parameter_bounds['additives']['names'],
        'AdditiveLoading': np.arange(parameter_bounds['additives']['min_mole_percent'],
                                    parameter_bounds['additives']['max_mole_percent']+parameter_bounds['additives']['step_size']/10.0,
                                    step=parameter_bounds['additives']['step_size']).tolist(),
        'Temperature': np.arange(parameter_bounds['temperature']['min'],
                                    parameter_bounds['temperature']['max']+parameter_bounds['temperature']['step_size']/10.0,
                                    step=parameter_bounds['temperature']['step_size']).tolist(),
        'Stage':parameter_bounds['light stages']['names']
    }
    df_one = pd.DataFrame({k: random.sample(v,1) for k, v in dict_listAll.items()})
    return df_one 
                            

def initial_batch_LM(parameter_bounds):
    # balanced design for the first 8 points:
    saveCSVFileName = 'df_init.csv'
    if os.path.exists(saveCSVFileName):
        df_init = pd.read_csv(saveCSVFileName)
    else:
        dict_Reagent_Solvent = {'Reagent': parameter_bounds['reagents']['names'], 'Solvent': parameter_bounds['solvents']['names']}
        df_init = expand_grid(dict_Reagent_Solvent)
        df_init = pd.concat([df_init, df_init], axis = 0).reset_index(drop=True).assign(Additive=parameter_bounds['additives']['names'])
        Reagent_equiv_all = np.arange(parameter_bounds['reagents']['min_equivalents'],
                                        parameter_bounds['reagents']['max_equivalents']+parameter_bounds['reagents']['step_size']/10.0,
                                        step=parameter_bounds['reagents']['step_size'])
        df_init = df_init.assign(Reagent_equiv=np.random.permutation([Reagent_equiv_all[round(x)] for x in np.percentile(range(len(Reagent_equiv_all)), [0, 25, 25, 50, 50, 75, 75, 100])]))
        AdditiveLoading_all = np.arange(parameter_bounds['additives']['min_mole_percent'],
                                        parameter_bounds['additives']['max_mole_percent']+parameter_bounds['additives']['step_size']/10.0,
                                        step=parameter_bounds['additives']['step_size'])
        df_init = df_init.assign(AdditiveLoading=np.random.permutation([AdditiveLoading_all[round(x)] for x in np.percentile(range(len(AdditiveLoading_all)), [0, 25, 25, 50, 50, 75, 75, 100])]))
        Temperature_all = np.arange(parameter_bounds['temperature']['min'],
                                        parameter_bounds['temperature']['max']+parameter_bounds['temperature']['step_size']/10.0,
                                        step=parameter_bounds['temperature']['step_size'])
        df_init = df_init.assign(Temperature=np.random.permutation([Temperature_all[round(x)] for x in np.percentile(range(len(Temperature_all)), [0, 25, 25, 50, 50, 75, 75, 100])]))
        df_init = df_init.assign(Stage=random.sample([1,2,2,3,3,4,4,5],8))
        df_init.to_csv(saveCSVFileName, index = False)
    return df_init
    
def optFun_LM(dat, parameter_bounds, n_start = 8):
    if dat.shape[0] < n_start:
        #df_init = initial_batch_LM(parameter_bounds)
        #dict_select = df_init.iloc[dat.shape[0]].to_dict()
        dict_select = initial_oneSample(parameter_bounds).iloc[0].to_dict()
        return dict_select   
    X, cols_X, y = preprocessing_lm(dat,X_only = False)
    # decide whether to include the quadratic terms 
    cols_subset = list(range(8))
    for j in np.arange(8, 12):
        max_corr_j = np.max([abs(np.corrcoef(X[:,i], X[:,j])[0,1]) for i in np.arange(1,8)])
        if max_corr_j < 0.99:
            cols_subset.append(j)
    X = X[:,cols_subset]
    z = transform_f2z(y)
    reg = LinearRegression().fit(X[:,1:], z)
    df_para = parameter_df(parameter_bounds)
    X_all, cols_X_new = preprocessing_lm(df_para,X_only = True)
    X_all = X_all[:,cols_subset]
    pred_z = reg.predict(X_all[:,1:])
    pred_z = np.squeeze(pred_z)
    top_index = np.argwhere(pred_z == np.amax(pred_z)) #np.argmax(pred_z)
    np.random.shuffle(top_index)
    if top_index.shape[0] > 1:
        dist2old = []
        for i in range(top_index.shape[0]):
            min_dist2old_i = np.min([np.linalg.norm(X_all[top_index[i],:]-X[j,:]) for j in range(X.shape[0])])
            dist2old.append(min_dist2old_i) 
        select_index = top_index[np.argmax(dist2old)]
        dist2old_select = np.max(dist2old)
    else:
        select_index = top_index[0]
        dist2old_select = np.min([np.linalg.norm(X_all[select_index,:]-X[j,:]) for j in range(X.shape[0])])
    if dist2old_select>0:
        df_select = df_para.iloc[select_index].assign(target=['pred_f'],target_value = transform_z2f(pred_z[select_index])).reset_index(drop=True)
        dict_select = df_select.iloc[0].to_dict()
        return dict_select
    else:
        print("Didn't find new candidate by LM method. Use RS instead.")
        return optFun_RS(dat, parameter_bounds)

def optFun_RS(dat, parameter_bounds, n_start = 8, n_candidates = 100):
    if dat.shape[0] < n_start:
        #df_init = initial_batch_LM(parameter_bounds)
        #dict_select = df_init.iloc[dat.shape[0]].to_dict()
        dict_select = initial_oneSample(parameter_bounds).iloc[0].to_dict()
        return dict_select
    X, cols_X_1 = preprocessing_lm(dat,X_only = True)
    df_para = parameter_df(parameter_bounds)
    candidates = random.sample(range(df_para.shape[0]), k = n_candidates)
    df_para = df_para.iloc[candidates]
    X_all, cols_X_2 = preprocessing_lm(df_para,X_only = True)
    dist2old = []
    for i in range(len(candidates)):
        min_dist2old_i = np.min([np.linalg.norm(X_all[i,:]-X[j,:]) for j in range(X.shape[0])])
        dist2old.append(min_dist2old_i)
    select_index = np.argmax(dist2old)
    #print(max(dist2old))
    df_select = df_para.iloc[[select_index]].reset_index(drop=True)
    dict_select = df_select.iloc[0].to_dict()
    return dict_select

def optimization(dat, parameter_bounds, method = 'RS', n_start = 8):
    if method in ['RS','LM']:
        print("Round: ", dat.shape[0]+1, "\t Optimization method:", method, "\n")
    else:
        print("Round: ", dat.shape[0]+1, "\t Warning: Optimization method", method, "is not available! \n Use Random Search instead. \n")
        method = 'RS'
    #
    if method == 'LM':
        new_para = optFun_LM(dat, parameter_bounds, n_start)
    elif method == 'BO':
        # To-Do
        return 0
    else: 
        # RS
        new_para = optFun_RS(dat, parameter_bounds, n_start)
    #print(new_para)
    return new_para