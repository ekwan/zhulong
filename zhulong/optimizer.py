"""
Inputs:
    dat: A pandas data frame of available experimental results for modeling
    parameter_bounds: A dictionary contains the min/max/step-size of all experimental conditions to be explored, 
                      produced by parameter_space.get_bounds_dict(), where the 'parameter_space' is a ParameterSpace class object 
    method: Choose an optimization method from the list ['RS','LM','BO']
        'RS' = Random Search
        'LM' = Polynomial regression modeling
        'BO' = Bayesian optimization
    aquisition: 
        For BO: Choose an aquisition function from the list ['EI','PI','UCB']
            'EI' = Expected Improvement
            'PI' = Probability of Improvement
            'UCB' = GP Upper Confidence Bound 
        For LM: Only one option 'Pred'
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
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel
from scipy.stats import norm, pearsonr
from optimizer_utils import *
import random
import os
import pickle

def initial_oneSample(parameter_bounds, seed):
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
    random.seed(seed)
    df_one = pd.DataFrame({k: random.sample(v,1) for k, v in dict_listAll.items()})
    return df_one 
                            

def initial_batch_8(parameter_bounds, seed):
    # balanced design for the first 8 points, with D-optimality criterion 
    saveCSVFileName = os.path.join('temp_initial_batch','df_init_seed_'+str(seed)+'.csv')
    if not os.path.isdir('temp_initial_batch'):
        os.mkdir('temp_initial_batch')
    if os.path.exists(saveCSVFileName):
        df_init = pd.read_csv(saveCSVFileName)
    else:
        random.seed(seed)
        dict_Reagent_Solvent = {'Reagent': parameter_bounds['reagents']['names'], 'Solvent': parameter_bounds['solvents']['names']}
        df_init_factors = expand_grid(dict_Reagent_Solvent)
        df_init_factors = df_init_factors.sample(frac = 1)
        Additive_all = np.random.permutation(parameter_bounds['additives']['names'])
        df_init_factors = pd.concat([df_init_factors, df_init_factors], axis = 0).reset_index(drop=True).assign(Additive=Additive_all)
        det_list = []
        df_init_list = []
        for i in range(1000):
            df_init_i = pd.DataFrame()
            Reagent_equiv_all = np.arange(parameter_bounds['reagents']['min_equivalents'],
                                            parameter_bounds['reagents']['max_equivalents']+parameter_bounds['reagents']['step_size']/10.0,
                                            step=parameter_bounds['reagents']['step_size'])
            df_init_i = df_init_i.assign(Reagent_equiv=np.random.permutation([Reagent_equiv_all[round(x)] for x in np.percentile(range(len(Reagent_equiv_all)), [0, 25, 25, 50, 50, 75, 75, 100])]))
            AdditiveLoading_all = np.arange(parameter_bounds['additives']['min_mole_percent'],
                                            parameter_bounds['additives']['max_mole_percent']+parameter_bounds['additives']['step_size']/10.0,
                                            step=parameter_bounds['additives']['step_size'])
            df_init_i = df_init_i.assign(AdditiveLoading=np.random.permutation([AdditiveLoading_all[round(x)] for x in np.percentile(range(len(AdditiveLoading_all)), [0, 25, 25, 50, 50, 75, 75, 100])]))
            Temperature_all = np.arange(parameter_bounds['temperature']['min'],
                                            parameter_bounds['temperature']['max']+parameter_bounds['temperature']['step_size']/10.0,
                                            step=parameter_bounds['temperature']['step_size'])
            df_init_i = df_init_i.assign(Temperature=np.random.permutation([Temperature_all[round(x)] for x in np.percentile(range(len(Temperature_all)), [0, 25, 25, 50, 50, 75, 75, 100])]))
            df_init_i = df_init_i.assign(Stage=random.sample([1,2,2,3,3,4,4,5],8))
            df_init_i = df_init_i.assign(Reagent_DBDMH = np.array(df_init_factors['Reagent'] == 'DBDMH', dtype = 'int8'))
            df_init_i = df_init_i.assign(Solvent_DMC = np.array(df_init_factors['Solvent'] == 'DMC', dtype = 'int8'))
            X = df_init_i.to_numpy()
            for j in range(X.shape[1]):
                X[:,j] = (X[:,j]-np.mean(X[:,j]))/np.std(X[:,j])
            df_init_list.append(df_init_i)
            det_list.append(np.linalg.det(np.matmul(np.transpose(X),X)))
        ind_max = det_list.index(max(det_list))
        df_init = pd.concat([df_init_factors, df_init_list[ind_max]], axis = 1).reset_index(drop=True)
        df_init.drop(columns=['Reagent_DBDMH','Solvent_DMC'], inplace = True)
        df_init.to_csv(saveCSVFileName, index = False)
    return df_init
    
def optFun_LM(dat, parameter_bounds, seed, initMethod='random', n_start=8, model_dir='models', acquisition='Pred'):
    if (initMethod=='random') & (dat.shape[0] < n_start):
        print("Random search for initial points.")
        return optFun_RS(dat, parameter_bounds, seed)
    if (initMethod=='batch8') & (dat.shape[0] < 8):
        print("Batch of 8 diverse initial points with seed:", seed)
        df_init = initial_batch_8(parameter_bounds, seed)
        dict_select = df_init.iloc[dat.shape[0]].to_dict()
        return dict_select 
    random.seed(seed+5*dat.shape[0]+321)
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
    saveModelPath = os.path.join(model_dir, 'LM_'+initMethod+'_'+acquisition+'_'+str(seed)+'_Rd_'+str(dat.shape[0]+1)+'.pkl')
    with open(saveModelPath,'wb') as f:
        pickle.dump(reg,f)
        pickle.dump(X,f)
        pickle.dump(z,f)
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
        df_select = df_para.iloc[select_index].assign(target=[acquisition],
                                                      target_value = pred_z[select_index], 
                                                      target_f=transform_z2f(pred_z[select_index])).reset_index(drop=True)
        dict_select = df_select.iloc[0].to_dict()
        return dict_select
    else:
        print("Didn't find new candidate by LM method. Use RS instead.")
        return optFun_RS(dat, parameter_bounds, seed+3*dat.shape[0]+116)

def optFun_BO(dat, parameter_bounds, seed, initMethod='random', n_start=8, model_dir='models', acquisition='EI'):
    if (initMethod=='random') & (dat.shape[0] < n_start):
        print("Random search for initial points.")
        return optFun_RS(dat, parameter_bounds, seed)
    if (initMethod=='batch8') & (dat.shape[0] < 8):
        print("Batch of 8 diverse initial points with seed:", seed)
        df_init = initial_batch_8(parameter_bounds, seed)
        dict_select = df_init.iloc[dat.shape[0]].to_dict()
        return dict_select 
    random.seed(seed+6*dat.shape[0]+45)
    X, cols_X, y = preprocessing_GP(dat,X_only = False)
    z = transform_f2z(y)
    kernel = 1.0 * RBF(length_scale=1e-1, length_scale_bounds=(1e-2, 1e3)) + WhiteKernel(noise_level=1e-2, noise_level_bounds=(1e-10, 1e1))
    gpr = GaussianProcessRegressor(kernel=kernel, alpha=0.0, n_restarts_optimizer = 10, normalize_y = True, random_state = 123)
    gpr.fit(X, z)
    saveModelPath = os.path.join(model_dir, 'BO_'+initMethod+'_'+acquisition+'_'+str(seed)+'_Rd_'+str(dat.shape[0]+1)+'.pkl')
    with open(saveModelPath,'wb') as f:
        pickle.dump(gpr,f)
        pickle.dump(X,f)
        pickle.dump(z,f)
    df_para = parameter_df(parameter_bounds)
    X_all, cols_X_new = preprocessing_GP(df_para,X_only = True)
    pred_z, pred_sd_z = gpr.predict(X_all, return_std=True)
    pred_z = np.squeeze(pred_z)
    z_best = max(z)
    gamma_z = (z_best - pred_z)/pred_sd_z
    # Acquisition function
    if acquisition=='EI':
        # Expected Improvement (EI) 
        acq_values = pred_sd_z*(norm.pdf(gamma_z)-gamma_z*(1-norm.cdf(gamma_z))) # EI
    elif acquisition=='PI':
        # Probability of Improvement (PI)
        acq_values = 1-norm.cdf(gamma_z) # PI
    else:
        # GP Upper Confidence Bound (UCB)
        alpha = 0.8 # use 80% prediction interval
        acq_values = pred_z+norm.ppf(1-(1-alpha)/2.0)*pred_sd_z # upr_z
    top_index = np.argwhere(acq_values == np.amax(acq_values)) #np.argmax(acq_values)
    if top_index.shape[0] > 1:
        np.random.shuffle(top_index)
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
        df_select = df_para.iloc[select_index].assign(target=[acquisition],target_value = acq_values[select_index]).reset_index(drop=True)
        if acquisition=='EI':
            df_select=df_select.assign(target_f = transform_z2f(z_best+acq_values[select_index]))
        elif acquisition=='PI':
            df_select=df_select.assign(target_f = 100*acq_values[select_index])
        else:
            df_select=df_select.assign(target_f = transform_z2f(acq_values[select_index]))
        dict_select = df_select.iloc[0].to_dict()
        return dict_select
    else:
        print("Didn't find new candidate by BO method. Use RS instead.")
        return optFun_RS(dat, parameter_bounds, seed+2*dat.shape[0]+885)
def optFun_RS(dat, parameter_bounds, seed, n_start = 1, n_candidates = 100):
    if dat.shape[0] < n_start:
        #df_init = initial_batch_LM(parameter_bounds)
        #dict_select = df_init.iloc[dat.shape[0]].to_dict()
        dict_select = initial_oneSample(parameter_bounds, seed).iloc[0].to_dict()
        return dict_select
    X, cols_X_1 = preprocessing_lm(dat,X_only = True)
    df_para = parameter_df(parameter_bounds)
    random.seed(seed+17*dat.shape[0]+99)
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

def optimization(dat, parameter_bounds, method = 'BO', seed = 0, initMethod = 'random', n_start=8, acquisition='EI'):
    assert(initMethod in ['random','batch8'])
    if method in ['RS','LM','BO']:
        print("Round: ", dat.shape[0]+1, "\t Optimization method:", method, "\t Aquisition function: ", acquisition)
    else:
        print("Round: ", dat.shape[0]+1, "\t Warning: Optimization method", method, "is not available! \n Use Random Search instead. ")
        method = 'RS'
    # create folder to save trained models at each iteration
    saveModelFolder='models'
    if not os.path.exists(saveModelFolder):
        os.makedirs(saveModelFolder)
    if method == 'LM':
        assert(acquisition in ['Pred'])
        new_para = optFun_LM(dat, parameter_bounds, seed, initMethod, n_start, saveModelFolder, acquisition)
    elif method == 'BO':
        assert(acquisition in ['EI','PI','UCB'])
        new_para = optFun_BO(dat, parameter_bounds, seed, initMethod, n_start, saveModelFolder, acquisition)
    else: 
        # RS
        new_para = optFun_RS(dat, parameter_bounds, seed)
    #print(new_para)
    return new_para