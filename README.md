# Zhulong

Autonomous optimization code.  This code provides code to run autonomous optimization with a ChemSpeed Robot and Agilent HPLC system.

### A Trivial Optimization

Here's how to run two points.  You can find this code in `zhulong/zhulong.py`.

Here are the required imports:

```
from experiment import Reagent, ParameterSpace, Experiment
from chemstation import Peak, get_yield_function
from chemspeed import ChemSpeed
```

First, define the starting material and some reagents:

```
# define the starting material
starting_material = Reagent(name="starting material", abbreviation="SM", min_volume=100, max_volume=100)

# define the reagents
# assume we will pick one of these
NBS = Reagent(name="NBS", abbreviation="NBS", min_volume=73, max_volume=100)
DBDMH = Reagent(name="DBDMH", abbreviation="DBDMH", min_volume=73, max_volume=100)
reagents = [NBS, DBDMH]
```

Starting materials must have a fixed volume.  Reagents can have any volume in their allowed range; makeup solvent will be computed later to ensure every experiment has a constant volume.

Now, let's define the additives and solvents:

```
# define the additives
# assume we will pick one of these
hydrochloric_acid = Reagent(name="hydrochloric acid", abbreviation="HCl", min_volume=2, max_volume=50)
sulfuric_acid = Reagent(name="sulfuric acid", abbreviation="H2SO4", min_volume=2, max_volume=50)
picolinic_acid = Reagent(name="picolinic acid", abbreviation="Picolinic", min_volume=2, max_volume=50)
phenylphosphonic_acid = Reagent(name="phenylphosphonic acid", abbreviation="Phenylphosphonic", min_volume=2, max_volume=50)
phosphoric_acid = Reagent(name="phosphoric acid", abbreviation="Phosphoric", min_volume=2, max_volume=50)
lactic_acid = Reagent(name="lactic acid", abbreviation="Lactic", min_volume=2, max_volume=50)
acetic_acid = Reagent(name="acetic acid", abbreviation="Acetic", min_volume=2, max_volume=50)
water = Reagent(name="water", abbreviation="Water", min_volume=2, max_volume=50)
additives = [ hydrochloric_acid, sulfuric_acid, picolinic_acid, phenylphosphonic_acid,
              phosphoric_acid, lactic_acid, acetic_acid, water ]

# define the solvents
# assume we will pick one of these
solvents = [ "MeCN", "DMC" ]
```

Notice that the starting materials, reagents, and additives are all of type `Reagent` but solvents are of type `str`.

Now, we'll define the parameter space:

```
parameter_space = ParameterSpace(starting_material, reagents, solvents,
                                    additives, light_stages=5,
                                    total_volume=200)   # add solvent to make the final volume in each experiment 200 uL
```

Here's how to define an experiment in the parameter space:

```
experiment1 = Experiment(parameter_space,                   # defines the parameter space this experiment is helping to explore
                         solvent="DMC",                     # stocks will automatically be added in the correct solvent
                         temperature=25,                    # in Celsius
                         starting_material_volume=50,       # in uL
                         reagent="NBS",                     # one of the reagents in the ParameterSpace (also accepted:
                                                            # list index as zero-indexed integer or Reagent object)
                         reagent_volume=30,                 # in uL
                                                            # (also accepted: zero-indexed list index)
                         additive="lactic acid",            # one of the additives in the ParameterSpace
                                                            # (also accepted: zero-indexed list index or Reagent object)
                         additive_volume=40,                # in uL
                         light_stage=5)                     # int: 1-5
```

These values will be automatically converted into the correct format in the ChemSpeed csv file.  The order of columns is hardcoded into `experiment.get_row()`.  This function produces the csv column names (`headings`) and the values for the row (`values`):

```
headings, values = experiment.get_row(sampling=1)
for h,v in zip(headings,values):
    print(f"{h} : {v}")
```

Before we can run the experiment, we'll need to tell *zhulong* how to compute the yield.  First, we'll need to tell it where the HPLC peaks are:

```
starting_material_peak = Peak(name="starting material", min_retention_time=1.28, max_retention_time=1.48)
product_peak = Peak(name="bromo product", min_retention_time=1.53, max_retention_time=1.73)
internal_standard_peak = Peak(name="internal standard", min_retention_time=2.05, max_retention_time=2.25)
```

Then, we need to define a function that converts a ChemStation report into a yield.  In this test, we'll just take the "yield" as the raw product integral and multiply by the response factor (1.0) here.  In production, we should measure the response factor and calculate the yield as response_factor * product/IS:

```
# for testing, just report the raw integral of the product peak
# but could also convert to chemical yield via response factor
yield_function = get_yield_function(product_peak, internal_standard_peak=None, response_factor=1.0)
#yield_function = get_yield_function(product_peak, internal_standard_peak, response_factor=10.0)
```

Next, create an object to represent the ChemSpeed robot:

```
chemspeed = ChemSpeed(chemspeed_csv_filename="temp.csv",
                      chemstation_folder="../data/E-Z 2021-09-30 17-18-47",
                      yield_function=yield_function,
                      overwrite_existing_chemspeed_csv=True,
                      ignore_existing_chemstation_folders=False,
                      polling_interval=1)
```

The volumes will be written to `temp.csv`

Finally, to run the experiment:

```
chemspeed.run_experiment(experiment1)
```

*Zhulong* will check `chemstation_folder` for new `.D` folders every `polling_interval` seconds.  If `ignore_existing_chemstation_folders` is set to `True`, then all the `.D` folders that exist when the `ChemSpeed` object is initialized will be ignored.

Note that `plateau_function` is not set here, which means it takes on its default value of `None`.  When there is no plateau function, a plateau is assumed to have occurred after one sample and a `0` will be written to the `Sampling` column in the csv file after the first new `.D` folder is parsed.  Optionally, you can provide a function that takes a list of yields and returns True if we have reached a plateau.  If a plateau has not been reached, then sampling will continue (i.e., by leaving the csv file with `Sampling` set to 1).

You can access previous results in `chemspeed.history`, which is a `dict` in which the keys are `Experiment` objects and the values are lists of yields (as `float`).

(No code for optimizers is provided.  Currently, some debugging message will be printed when `run_experiment()` is executed.)

### Authors

Melodie Christensen, Eugene Kwan, Yuting Xu
