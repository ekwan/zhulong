# Zhulong

Autonomous optimization code.  This code provides code to run autonomous optimization with a ChemSpeed Robot and an Agilent UPLC system.

### An Autonomous Optimization

Here's how to set up and run autonomous optimization experiments.  You can find this code in `zhulong/zhulong_optimizer.py`.

Here are the required imports:

```
from experiment import Reagent, ParameterSpace, Experiment
from chemstation import Peak, get_yield_function
from chemspeed import ChemSpeed
```

First, define the starting material and reagents:

```
# define the starting material
starting_material = Reagent(
    name="starting material",
    abbreviation="SM",
    min_volume=80,
    max_volume=80,
    concentration=0.250)

# define the reagents
# assume we will pick one of these
# reagents should have the same min/max volumes and stock solution
# concentrations
min_volume = 80       # uL
max_volume = 120      # uL
concentration = 0.250  # M
NBS = Reagent(
    name="NBS",
    abbreviation="NBS",
    min_volume=min_volume,
    max_volume=max_volume,
    concentration=concentration)
DBDMH = Reagent(
    name="DBDMH",
    abbreviation="DBDMH",
    min_volume=min_volume,
    max_volume=max_volume,
    concentration=concentration)
reagents = [NBS, DBDMH]
```

Starting materials must have a fixed volume.  Reagents can have any volume in their allowed range; makeup solvent will be computed later to ensure every experiment has a constant volume.

Now, let's define the additives and solvents:

```
# define the additives
# assume we will pick one of these
# additives should have the same min/max volumes and stock solution
# concentrations
min_volume = 2        # uL
max_volume = 50       # uL
concentration = 0.100  # M
hydrochloric_acid = Reagent(
    name="hydrochloric acid",
    abbreviation="HCl",
    min_volume=min_volume,
    max_volume=max_volume,
    concentration=concentration)
sulfuric_acid = Reagent(
    name="sulfuric acid",
    abbreviation="H2SO4",
    min_volume=min_volume,
    max_volume=max_volume,
    concentration=concentration)
picolinic_acid = Reagent(
    name="picolinic acid",
    abbreviation="Picolinic",
    min_volume=min_volume,
    max_volume=max_volume,
    concentration=concentration)
phenylphosphonic_acid = Reagent(
    name="phenylphosphonic acid",
    abbreviation="Phenylphosphonic",
    min_volume=min_volume,
    max_volume=max_volume,
    concentration=concentration)
phosphoric_acid = Reagent(
    name="phosphoric acid",
    abbreviation="Phosphoric",
    min_volume=min_volume,
    max_volume=max_volume,
    concentration=concentration)
lactic_acid = Reagent(
    name="lactic acid",
    abbreviation="Lactic",
    min_volume=min_volume,
    max_volume=max_volume,
    concentration=concentration)
acetic_acid = Reagent(
    name="acetic acid",
    abbreviation="Acetic",
    min_volume=min_volume,
    max_volume=max_volume,
    concentration=concentration)
water = Reagent(
    name="water",
    abbreviation="Water",
    min_volume=min_volume,
    max_volume=max_volume,
    concentration=concentration)
additives = [
    hydrochloric_acid,
    sulfuric_acid,
    picolinic_acid,
    phenylphosphonic_acid,
    phosphoric_acid,
    lactic_acid,
    acetic_acid,
    water]

# define the solvents
# assume we will pick one of these
solvents = ["MeCN", "DMC"]
```

Notice that the starting materials, reagents, and additives are all of type `Reagent` but solvents are of type `str`.

Now, we'll define the parameter space:

```
# define the parameter space to optimize over
parameter_space = ParameterSpace(starting_material, reagents, solvents,
                                 additives, light_stages=5,
                                 total_volume=350,
                                 # add solvent to make the final volume in each
                                 # experiment 250 uL
                                 min_temperature=5,         # in C
                                 max_temperature=35,       # in C
                                 temperature_step_size=5   # in C
                                 )

# for convenience, here is a dictionary containing the bounds
# bounds_dict = parameter_space.get_bounds_dict()
# bounds_dict_string = json.dumps(bounds_dict, indent=2)
# print(bounds_dict_string)
```

Now, we'll define the UPLC peaks of interest:

```
# define UPLC peaks
starting_material_peak = Peak(
    name="starting material",
    min_retention_time=0.92,
    max_retention_time=1.07)
product_peak = Peak(
    name="bromo product",
    min_retention_time=1.20,
    max_retention_time=1.35)
dibromo_peak = Peak(
    name="dibromo product",
    min_retention_time=1.42,
    max_retention_time=1.52)
# internal_standard_peak = Peak(
#     name="internal standard",
#     min_retention_time=1.00,
#     max_retention_time=1.10)
```

Then, we need to calculate the chemical yield:

```
# convert the ratio of the product peak over the internal standard peak to the chemical yield via the yield function
# in this project, we simply report UPLC area % of the product peak
yield_function = get_yield_function(
    product_peak,
    internal_standard_peak=None,
    response_factor=1.0)
```

Next, we create an object to represent the ChemSpeed robot and define the yield and plateau functions:

```
chemspeed = ChemSpeed(
    chemspeed_csv_filename="Z:\\TEST\\Closed Loop APO.csv",
    chemstation_folder="C:\\Users\\Public\\Documents\\ChemStation\\2\\Data\\Chemspeed\\Chemspeed*\\",
    yield_function=yield_function,
    plateau_function=plateau_function_1,
    overwrite_existing_chemspeed_csv=True,
    ignore_existing_chemstation_folders=True,
    polling_interval=1)
```

Dispense volumes will calculated and written to `Closed Loop APO.csv` for importing into the Chemspeed robot software once the experimental parameters are defined by the optimizer.

Now we define the optimization method, initialization method, and acquisition function:

```
# get the parameter bounds and define the optimization method
# a unique acquisition function can be defined for each optimization round
parameter_bounds = parameter_space.get_bounds_dict()
method = 'BO'  # or 'RS' or 'LM'
initMethod = 'batch8'  # 'random' or 'batch8' for LM and BO
acquisition = 'PI'  # 'EI', 'PI' or 'UCB' for BO, and 'Pred' for LM
acquisition_list = [
    'PI',
    'UCB',
    'EI',
    'PI',
    'UCB',
    'EI',
    'PI',
    'UCB',
    'EI',
    'PI',
    'UCB',
    'EI',
    'PI',
    'UCB',
    'EI',
    'PI',
    'UCB',
    'EI',
    'PI',
    'UCB',
    'EI',
    'PI',
    'UCB',
    'EI',
    'PI',
    'UCB',
    'EI',
    'PI',
    'UCB',
    'EI',
    'PI',
    'UCB',
    'EI',
    'PI',
    'UCB',
    'EI',
    'PI',
    'UCB',
    'EI',
    'PI',
    'UCB',
    'EI',
    'PI',
    'UCB',
    'EI',
    'PI',
    'UCB',
    'EI']
n_total = 48  # total iterations
all_experiments = []
```

Then we get a new set of parameters for each experimental round:
```
dat = pd.DataFrame()
f_best = 0
for r in range(n_total):
    # get a set of new parameters for round r
    # para_r = optimization(dat, parameter_bounds,method, seed = 1234, initMethod = initMethod)
    para_r = optimization(
        dat,
        parameter_bounds,
        method,
        seed=1234,
        initMethod=initMethod,
        acquisition=acquisition_list[r])
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
    # print the experiment for this round
    print("running experiment")
    print(experiment_r)
    chemspeed.run_experiment(experiment_r)
    history = experiment_r.history
    print("final result is:")
    print(f"{history['times']=}")
    print(f"{history['values']=}")
    # save the experiment for this round
    all_experiments.append(experiment_r)
    with open("APO_timecourse.pkl", "wb") as file:
        pickle.dump(all_experiments, file)
    f_r = history['values'][-1]
    f_best = max(f_best, f_r)
    df_para_r = df_para_r.assign(
        f=f_r, f_best=f_best, Round=(
            r + 1), method=method)
    dat = pd.concat([dat, df_para_r], ignore_index=True)
    print("Round: ", r + 1, "f: ", f_r, "f_best", f_best, "\n")
    dat.to_csv("dat.csv")
```

The output for each experimental round is:

```
Round:  1 	 Optimization method: BO 	 Aquisition function:  PI
Batch of 8 diverse initial points with seed: 1234
running experiment
solv=DMC (164 uL) temp=20°C SM = 80 uL reagent=DBDMH (vol=80 uL, 1.00 equiv) additive=H2SO4 (vol=26 uL, 13 mol%) light=2
```  

The code checks the `chemstation_folder` for new `.D` folders every `polling_interval` seconds.  If `ignore_existing_chemstation_folders` is set to `True`, then all the `.D` folders that exist when the `ChemSpeed` object is initialized are ignored.

Once the plateau algorithm determines that a reaction plateau is reached, the final result is reported.

The output for each experimental round is:   

```
found folder C:\Users\Public\Documents\ChemStation\2\Data\Chemspeed\Chemspeed 2022-05-12 16-00-50 - 10 Rounds - Copy\0312650-0725-001.D
parsing C:\Users\Public\Documents\ChemStation\2\Data\Chemspeed\Chemspeed 2022-05-12 16-00-50 - 10 Rounds - Copy\0312650-0725-001.D

Found C:\Users\Public\Documents\ChemStation\2\Data\Chemspeed\Chemspeed 2022-05-12 16-00-50 - 10 Rounds - Copy\0312650-0725-001.D/Report.TXT...parsing...
yield is 67.4741 (recorded at 1666358105.6932242 s since the epoch)
Plateau Algorithm No.1
plateaued=False

found folder C:\Users\Public\Documents\ChemStation\2\Data\Chemspeed\Chemspeed 2022-05-12 16-00-50 - 10 Rounds - Copy\0312650-0725-002.D
parsing C:\Users\Public\Documents\ChemStation\2\Data\Chemspeed\Chemspeed 2022-05-12 16-00-50 - 10 Rounds - Copy\0312650-0725-002.D

Found C:\Users\Public\Documents\ChemStation\2\Data\Chemspeed\Chemspeed 2022-05-12 16-00-50 - 10 Rounds - Copy\0312650-0725-002.D/Report.TXT...parsing...
yield is 81.8419 (recorded at 1666358106.1884985 s since the epoch)
Plateau Algorithm No.1
plateaued=False

found folder C:\Users\Public\Documents\ChemStation\2\Data\Chemspeed\Chemspeed 2022-05-12 16-00-50 - 10 Rounds - Copy\0312650-0725-003.D
parsing C:\Users\Public\Documents\ChemStation\2\Data\Chemspeed\Chemspeed 2022-05-12 16-00-50 - 10 Rounds - Copy\0312650-0725-003.D

Found C:\Users\Public\Documents\ChemStation\2\Data\Chemspeed\Chemspeed 2022-05-12 16-00-50 - 10 Rounds - Copy\0312650-0725-003.D/Report.TXT...parsing...
yield is 81.9339 (recorded at 1666358106.6911287 s since the epoch)
Plateau Algorithm No.1
plateaued=True
final result is:
history['times']=[1666358105.6932242, 1666358106.1884985, 1666358106.6911287]
history['values']=[67.4741, 81.8419, 81.9339]
Round:  1 f:  81.9339 f_best 81.9339 
```  

### Authors

Melodie Christensen, Eugene Kwan, Yuting Xu
