import numpy as np
import pandas as pd
from itertools import product
import math
import pickle


"""
Transformations between f and z
    f is experimental outcome, a percentage between 0 and 100
    z is transformed variable for modeling， whose domain is real line
    Example:
        f = np.array([10, 20, 40, 60, 80, 90])
        z=transform_f2z(f)
        transform_z2f(z)
"""

def transform_f2z(f):
    if isinstance(f, np.ndarray) and (f.shape[0] > 1):
        z = [transform_f2z(ff) for ff in f]
        z = np.array(z)[:, np.newaxis]
    else:
        f = f / 100.0
        z = math.log(f / (1 - f))
    return z


def transform_z2f(z):
    if isinstance(z, np.ndarray) and (z.shape[0] > 1):
        f = [transform_z2f(zz) for zz in z]
        f = np.array(f)[:, np.newaxis]
    else:
        f = 100.0 / (1 + math.exp(-z))
    return f


"""
Full factorial design:
    Input: A dictionary which lists all levels of each experimental condition
    Output: A pandas data frame consists of all possible combinations of levels for all factors (outer product)
    Example:
        dict = {'A': [0,1,2], 'B': ['b1', 'b2']}
        expand_grid(dict)
"""


def expand_grid(dictionary):
    return pd.DataFrame([row for row in product(*dictionary.values())],
                        columns=dictionary.keys())


"""
Given a dictionary of parameter_bounds, generate a data frame of all possible combinations.
Example:

with open('parameter_bounds.pkl', 'rb') as f:
    parameter_bounds = pickle.load(f)

df_para = parameter_df(parameter_bounds) # shape: (3292800, 7)
"""


def parameter_df(parameter_bounds):
    dict_listAll = {
        'Reagent': parameter_bounds['reagents']['names'],
        'Reagent_equiv': np.arange(
            parameter_bounds['reagents']['min_equivalents'],
            parameter_bounds['reagents']['max_equivalents'] +
            parameter_bounds['reagents']['step_size'] /
            10.0,
            step=parameter_bounds['reagents']['step_size']),
        'Solvent': parameter_bounds['solvents']['names'],
        'Additive': parameter_bounds['additives']['names'],
        'AdditiveLoading': np.arange(
            parameter_bounds['additives']['min_mole_percent'],
            parameter_bounds['additives']['max_mole_percent'] +
            parameter_bounds['additives']['step_size'] /
            10.0,
            step=parameter_bounds['additives']['step_size']),
        'Temperature': np.arange(
            parameter_bounds['temperature']['min'],
            parameter_bounds['temperature']['max'] +
            parameter_bounds['temperature']['step_size'] /
            10.0,
            step=parameter_bounds['temperature']['step_size']),
        'Stage': parameter_bounds['light stages']['names']}
    output = expand_grid(dict_listAll)
    return (output)


"""
Preprocessing of experimental data for modeling
"""


def preprocessing_lm(dat, X_only=True):
    # Prepare data for polynomial regression
    # Output X has 12 columns:
    #   First column is ones for intercept
    #   2 factors: Reagent, Solvent, each with two levels
    #   1 numeric linear: AdditiveLoading
    #   4 numeric quadratic: Additive_pKa, Temperature, Stage, Reagent_equiv
    # add numerical descriptor for Additive if missing
    if 'Additive_pKa' not in dat:
        # Load numerical descriptor for the Additive factor
        df_Additive_desc = pd.read_csv("../data/df_Additive_desc.csv")
        # To match the hard-coded names in zhulong:
        df_Additive_desc.rename(
            columns={
                'Additive': 'Additive_oldName'},
            inplace=True)
        df_Additive_desc = pd.concat(
            [
                pd.DataFrame(
                    {
                        'Additive': [
                            'hydrochloric acid',
                            'sulfuric acid',
                            'picolinic acid',
                            'phenylphosphonic acid',
                            'phosphoric acid',
                            'lactic acid',
                            'acetic acid',
                            'water']}),
                df_Additive_desc],
            axis=1,
            ignore_index=False)
        dat = pd.merge(dat, df_Additive_desc, how='left', on='Additive')
    # only keep Additive with Additive_pKa descriptor
    dat = dat[~dat['Additive_pKa'].isnull()]
    # scale some numerical variables to avoid large numbers
    #   Temperature
    Temperature_scaled = dat.loc[:, "Temperature"] / 40.0
    dat.loc[:, "Temperature"] = Temperature_scaled
    #   AdditiveLoading
    AdditiveLoading_scaled = dat.loc[:, "AdditiveLoading"] / 25.0
    dat.loc[:, "AdditiveLoading"] = AdditiveLoading_scaled
    #   Additive_pKa
    Additive_pKa_scaled = dat.loc[:, "Additive_pKa"] / 10.0
    dat.loc[:, "Additive_pKa"] = Additive_pKa_scaled
    #   Stage
    Stage_scaled = dat.loc[:, "Stage"] / 5.0
    dat.loc[:, "Stage"] = Stage_scaled
    # create one-hot encoding to replace factors
    """
    # need to revise if there are more levels
    dat = pd.get_dummies(dat, columns = ['Reagent','Solvent'],drop_first = True)
    colnames_Reagent = [s for s in list(dat) if s.startswith("Reagent_") and (s!='Reagent_equiv')]
    colnames_Solvent = [s for s in list(dat) if s.startswith("Solvent_")]
    """
    dat = dat.assign(
        Reagent_DBDMH=np.array(
            dat['Reagent'] == 'DBDMH',
            dtype='int8'))
    dat = dat.assign(
        Solvent_DMC=np.array(
            dat['Solvent'] == 'DMC',
            dtype='int8'))
    colnames_Reagent = ['Reagent_DBDMH']
    colnames_Solvent = ['Solvent_DMC']
    X = dat[colnames_Reagent + colnames_Solvent + ["AdditiveLoading",
                                                   "Additive_pKa", "Temperature", "Reagent_equiv", "Stage"]].to_numpy()
    # add intercept column and quadratic terms for Additive_pKa, Temperature,
    # Reagent_equiv and Stage to X
    Additive_pKa_quadratic = dat[["Additive_pKa"]].to_numpy()
    Additive_pKa_quadratic = Additive_pKa_quadratic * Additive_pKa_quadratic
    Temperature_quadratic = dat[["Temperature"]].to_numpy()
    Temperature_quadratic = Temperature_quadratic * Temperature_quadratic
    Reagent_equiv_quadratic = dat[["Reagent_equiv"]].to_numpy()
    Reagent_equiv_quadratic = Reagent_equiv_quadratic * Reagent_equiv_quadratic
    Stage_quadratic = dat[["Stage"]].to_numpy()
    Stage_quadratic = Stage_quadratic * Stage_quadratic
    X = np.column_stack(
        (np.ones(
            shape=(
                X.shape[0],
                1)),
         X,
         Additive_pKa_quadratic,
         Temperature_quadratic,
         Reagent_equiv_quadratic,
         Stage_quadratic))
    cols_X = ['ones'] + colnames_Reagent + colnames_Solvent + ["AdditiveLoading",
                                                               "Additive_pKa",
                                                               "Temperature",
                                                               "Reagent_equiv",
                                                               "Stage"] + ['Additive_pKa_2',
                                                                           'Temperature_2',
                                                                           'Reagent_equiv_2',
                                                                           'Stage_2']
    if X_only:
        return X, cols_X
    else:
        y = dat[["f"]].to_numpy()
        return X, cols_X, y


def preprocessing_GP(dat, X_only=True):
    # Prepare data for Gaussian Process regression
    # Output X has 7 columns:
    #   2 factors: Reagent, Solvent (each with two levels)
    #   5 numeric: AdditiveLoading, Additive_pKa, Temperature, Stage, Reagent_equiv
    # add numerical descriptor for Additive if missing
    if 'Additive_pKa' not in dat:
        # Load numerical descriptor for the Additive factor
        df_Additive_desc = pd.read_csv("../data/df_Additive_desc.csv")
        # To match the hard-coded names in zhulong:
        df_Additive_desc.rename(
            columns={
                'Additive': 'Additive_oldName'},
            inplace=True)
        df_Additive_desc = pd.concat(
            [
                pd.DataFrame(
                    {
                        'Additive': [
                            'hydrochloric acid',
                            'sulfuric acid',
                            'picolinic acid',
                            'phenylphosphonic acid',
                            'phosphoric acid',
                            'lactic acid',
                            'acetic acid',
                            'water']}),
                df_Additive_desc],
            axis=1,
            ignore_index=False)
        dat = pd.merge(dat, df_Additive_desc, how='left', on='Additive')
    # only keep Additive with Additive_pKa descriptor
    dat = dat[~dat['Additive_pKa'].isnull()]
    # scale some numerical variables to avoid large numbers
    #   Temperature
    Temperature_scaled = dat.loc[:, "Temperature"] / 40.0
    dat.loc[:, "Temperature"] = Temperature_scaled
    #   AdditiveLoading
    AdditiveLoading_scaled = dat.loc[:, "AdditiveLoading"] / 25.0
    dat.loc[:, "AdditiveLoading"] = AdditiveLoading_scaled
    #   Additive_pKa
    Additive_pKa_scaled = dat.loc[:, "Additive_pKa"] / 10.0
    dat.loc[:, "Additive_pKa"] = Additive_pKa_scaled
    #   Stage
    Stage_scaled = dat.loc[:, "Stage"] / 5.0
    dat.loc[:, "Stage"] = Stage_scaled
    # create one-hot encoding to replace factors
    dat = dat.assign(
        Reagent_DBDMH=np.array(
            dat['Reagent'] == 'DBDMH',
            dtype='int8'))
    dat = dat.assign(
        Solvent_DMC=np.array(
            dat['Solvent'] == 'DMC',
            dtype='int8'))
    colnames_Reagent = ['Reagent_DBDMH']
    colnames_Solvent = ['Solvent_DMC']
    X = dat[colnames_Reagent + colnames_Solvent + ["AdditiveLoading",
                                                   "Additive_pKa", "Temperature", "Reagent_equiv", "Stage"]].to_numpy()
    cols_X = colnames_Reagent + colnames_Solvent + \
        ["AdditiveLoading", "Additive_pKa", "Temperature", "Reagent_equiv", "Stage"]
    if X_only:
        return X, cols_X
    else:
        y = dat[["f"]].to_numpy()
        return X, cols_X, y
