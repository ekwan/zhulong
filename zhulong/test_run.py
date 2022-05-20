import json
from experiment import Reagent, ParameterSpace, Experiment
from chemstation import Peak, get_yield_function
from chemspeed import ChemSpeed
import pandas as pd
from optimizer import *
from optimizer_utils import *
from chemspeed_utils import plateau_function_1
import pickle

# define the starting material
starting_material = Reagent(name="starting material", abbreviation="SM", min_volume=80, max_volume=80, concentration=0.250)

# define the reagents
# assume we will pick one of these
# reagents should have the same min/max volumes and stock solution concentrations
min_volume = 80       # uL
max_volume = 120      # uL
concentration = 0.250 # M
NBS = Reagent(name="NBS", abbreviation="NBS", min_volume=min_volume, max_volume=max_volume, concentration=concentration)
DBDMH = Reagent(name="DBDMH", abbreviation="DBDMH", min_volume=min_volume, max_volume=max_volume, concentration=concentration)
reagents = [NBS, DBDMH]

# define the additives
# assume we will pick one of these
# additives should have the same min/max volumes and stock solution concentrations
min_volume = 2        # uL
max_volume = 50       # uL
concentration = 0.100 # M
hydrochloric_acid = Reagent(name="hydrochloric acid", abbreviation="HCl", min_volume=min_volume, max_volume=max_volume, concentration=concentration)
sulfuric_acid = Reagent(name="sulfuric acid", abbreviation="H2SO4", min_volume=min_volume, max_volume=max_volume, concentration=concentration)
picolinic_acid = Reagent(name="picolinic acid", abbreviation="Picolinic", min_volume=min_volume, max_volume=max_volume, concentration=concentration)
phenylphosphonic_acid = Reagent(name="phenylphosphonic acid", abbreviation="Phenylphosphonic", min_volume=min_volume, max_volume=max_volume, concentration=concentration)
phosphoric_acid = Reagent(name="phosphoric acid", abbreviation="Phosphoric", min_volume=min_volume, max_volume=max_volume, concentration=concentration)
lactic_acid = Reagent(name="lactic acid", abbreviation="Lactic", min_volume=min_volume, max_volume=max_volume, concentration=concentration)
acetic_acid = Reagent(name="acetic acid", abbreviation="Acetic", min_volume=min_volume, max_volume=max_volume, concentration=concentration)
water = Reagent(name="water", abbreviation="Water", min_volume=min_volume, max_volume=max_volume, concentration=concentration)
additives = [ hydrochloric_acid, sulfuric_acid, picolinic_acid, phenylphosphonic_acid,
              phosphoric_acid, lactic_acid, acetic_acid, water ]

# define the solvents
# assume we will pick one of these
solvents = [ "MeCN", "DMC" ]

# define the parameter space to optimize over
parameter_space = ParameterSpace(starting_material, reagents, solvents,
                                 additives, light_stages=5,
                                 total_volume=250,          # add solvent to make the final volume in each experiment 250 uL
                                 min_temperature=5,         # in C
                                 max_temperature=35,       # in C
                                 temperature_step_size=5   # in C
                                )

# for convenience, here is a dictionary containing the bounds
# bounds_dict = parameter_space.get_bounds_dict()
# bounds_dict_string = json.dumps(bounds_dict, indent=2)
# print(bounds_dict_string)

# define HPLC peaks
starting_material_peak = Peak(name="starting material", min_retention_time=0.92, max_retention_time=1.05)
product_peak = Peak(name="bromo product", min_retention_time=1.20, max_retention_time=1.35)
dibromo_peak = Peak(name="dibromo product", min_retention_time=1.42, max_retention_time=1.52)
#internal_standard_peak = Peak(name="internal standard", min_retention_time=1.40, max_retention_time=1.50)

# for testing, just report the raw integral of the product peak
# but could also convert to chemical yield via response factor
# in this run, we report LC area % of the product peak
yield_function = get_yield_function(product_peak, internal_standard_peak=None, response_factor=1.0)
#yield_function = get_yield_function(product_peak, internal_standard_peak, response_factor=10.0)

# represents the ChemSpeed robot
chemspeed = ChemSpeed(chemspeed_csv_filename="Z:\\TEST\\Closed Loop APO.csv",
                      chemstation_folder="C:\\Users\\Public\\Documents\\ChemStation\\2\\Data\\Chemspeed\\Chemspeed*\\",
                      yield_function=yield_function,
                      plateau_function=plateau_function_1,
                      overwrite_existing_chemspeed_csv=True,
                      ignore_existing_chemstation_folders=True,
                      polling_interval=1)

# Modification of last section "run the experiments" in zhulong.py:

parameter_bounds = parameter_space.get_bounds_dict()
method = 'LM'# or'RS'
initMethod = 'batch8' # 'random' or 'batch8' for LM
n_total = 48 # total iterations
all_experiments = []

dat = pd.DataFrame()
f_best = 0
for r in range(n_total):
    # get a set of new parameters for round r
    para_r = optimization(dat, parameter_bounds,method, seed = 1234, initMethod = initMethod)
    df_para_r = pd.DataFrame({k: [v] for k, v in para_r.items()})
    # define an experiment object for this round
    experiment_r = Experiment.create(
                         parameter_space,
                         solvent=para_r['Solvent'],
                         temperature=int(para_r['Temperature']),
                         starting_material_volume=80,
                         reagent=para_r['Reagent'],
                         reagent_equivalents=para_r['Reagent_equiv'],
                         additive=para_r['Additive'],
                         additive_mole_percent=para_r['AdditiveLoading'],
                         light_stage=para_r['Stage'])
    print("running experiment")
    print(experiment_r)
    chemspeed.run_experiment(experiment_r)
    history = experiment_r.history
    print("final result is:")
    print(f"{history['times']=}")
    print(f"{history['values']=}")
    all_experiments.append(experiment_r)
    with open("APO_timecourse.pkl","wb") as file:
        pickle.dump(all_experiments, file)
    f_r = history['values'][-1]
    f_best = max(f_best, f_r)
    df_para_r = df_para_r.assign(f = f_r, f_best = f_best)
    dat = pd.concat([dat,df_para_r],ignore_index = True)
    print("Round: ", r+1, "f: ", f_r,"f_best",f_best,"\n")

dat = dat.assign(Round = range(1,n_total+1))
dat = dat.assign(method = method)
dat.to_csv("dat.csv")