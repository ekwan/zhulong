import json
from experiment import Reagent, ParameterSpace, Experiment
from chemstation import Peak, get_yield_function
from chemspeed import ChemSpeed

# define the starting material
starting_material = Reagent(name="starting material", abbreviation="SM", min_volume=100, max_volume=100, concentration=0.630)

# define the reagents
# assume we will pick one of these
# reagents should have the same min/max volumes and stock solution concentrations
min_volume = 73       # uL
max_volume = 100      # uL
concentration = 0.945 # M
NBS = Reagent(name="NBS", abbreviation="NBS", min_volume=min_volume, max_volume=max_volume, concentration=concentration)
DBDMH = Reagent(name="DBDMH", abbreviation="DBDMH", min_volume=min_volume, max_volume=max_volume, concentration=concentration)
reagents = [NBS, DBDMH]

# define the additives
# assume we will pick one of these
# additives should have the same min/max volumes and stock solution concentrations
min_volume = 2        # uL
max_volume = 50       # uL
concentration = 0.315 # M
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
                                 min_temperature=-30,         # in C
                                 max_temperature=40,       # in C
                                 temperature_step_size=5   # in C
                                )

# for convenience, here is a dictionary containing the bounds
bounds_dict = parameter_space.get_bounds_dict()
bounds_dict_string = json.dumps(bounds_dict, indent=2)
print(bounds_dict_string)

# define an experiment using volumes
experiment1 = Experiment(parameter_space,                   # defines the parameter space this experiment is helping to explore
                         solvent="DMC",                     # stocks will automatically be added in the correct solvent
                         temperature=25,                    # in Celsius
                         starting_material_volume=100,      # in uL
                         reagent="NBS",                     # one of the reagents in the ParameterSpace (also accepted:
                                                            # list index as zero-indexed integer or Reagent object)
                         reagent_volume=75,                 # in uL
                                                            # (also accepted: zero-indexed list index)
                         additive="lactic acid",            # one of the additives in the ParameterSpace
                                                            # (also accepted: zero-indexed list index or Reagent object)
                         additive_volume=40,                # in uL
                         light_stage=5)                     # int: 1-5

# print out what would go into the chemspeed csv file
#headings, values = experiment.get_row(sampling=1)
#for h,v in zip(headings,values):
#    print(f"{h} : {v}")
#print()

# define an experiment using equivalents and mol%
experiment2 = Experiment.create(
                         parameter_space,                   # defines the parameter space this experiment is helping to explore
                         solvent="MeCN",                    # stocks will automatically be added in the correct solvent
                         temperature=25,                    # in Celsius
                         starting_material_volume=100,      # in uL
                         reagent="DBDMH",                   # one of the reagents in the ParameterSpace (also accepted:
                                                            # list index as zero-indexed integer or Reagent object)
                         reagent_equivalents=1.40,          # in equivalents
                                                            # (also accepted: zero-indexed list index)
                         additive="lactic acid",            # one of the additives in the ParameterSpace
                                                            # (also accepted: zero-indexed list index or Reagent object)
                         additive_mole_percent=5,           # in mol%
                         light_stage=2)                     # int: 1-5

# define HPLC peaks
starting_material_peak = Peak(name="starting material", min_retention_time=1.28, max_retention_time=1.48)
product_peak = Peak(name="bromo product", min_retention_time=1.53, max_retention_time=1.73)
internal_standard_peak = Peak(name="internal standard", min_retention_time=2.05, max_retention_time=2.25)

# for testing, just report the raw integral of the product peak
# but could also convert to chemical yield via response factor
yield_function = get_yield_function(product_peak, internal_standard_peak=None, response_factor=1.0)
#yield_function = get_yield_function(product_peak, internal_standard_peak, response_factor=10.0)

# represents the ChemSpeed robot
chemspeed = ChemSpeed(chemspeed_csv_filename="temp.csv",
                      chemstation_folder="../data/E-Z 2021-09-30 17-18-47/",
                      yield_function=yield_function,
                      plateau_function=None,
                      overwrite_existing_chemspeed_csv=True,
                      ignore_existing_chemstation_folders=False,
                      polling_interval=1)

# run the experiments
experiments = [experiment1, experiment2]
for e in experiments:
    print("running experiment")
    print(e)
    chemspeed.run_experiment(e)
    print()
    history = e.history
    print("final result is:")
    print(f"{history['times']=}")
    print(f"{history['values']=}")
    print()
    print("json looks like:")
    # to get experiment parameters in dictionary form
    # print(e.get_dict())
    print(e.get_json_string())
    print("--------")
