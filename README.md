# Zhulong

Autonomous optimization code.  This code provides code to run autonomous optimization with a ChemSpeed Robot and Agilent HPLC system.

### A Basic Optimization

Here's how to run the simplest kind of optimization: three random points from a two-parameter space.  This code is available at the bottom of `optimizer.py`.
First, define the parameter space as one continuous and one categorical parameter:

```
param1 = ContinuousParameter(name="temperature", min_value=-10.0, max_value=100.0)
param2 = CategoricalParameter(name="HCl", allowed_levels={"low":1.0, "medium":2.0, "high":3.0})
parameter_space = [ param1, param2 ]
```

Next, setup the compound to be monitored by HPLC:

```
starting_material = Compound(name="starting material", min_retention_time=1.28, max_retention_time=1.48)
product = Compound(name="bromo product", min_retention_time=1.53, max_retention_time=1.73)
internal_standard = Compound(name="internal standard", min_retention_time=2.05, max_retention_time=2.25)
```

Create a `DummyOptimizer` which just samples three random points.  There are already three example HPLC files in the repo, so we can trick the program into thinking that three new files are available.  (Files are parsed from oldest to newest, ignoring files that have already been seen.)

```
optimizer = DummyOptimizer(name="dummy optimizer", chemspeed_csv_filename=chemspeed_csv_filename,
                           chemstation_folder="/Users/kwaneu/research/zhulong/data/E-Z 2021-09-30 17-18-47",
                           parameter_space=parameter_space)
```

We'll also need to register the compounds with the optimizer:

```
optimizer.starting_material = starting_material
optimizer.product = product
optimizer.internal_standard = internal_standard
```

Finally, run the optimization:

```
optimizer.run()
```

Here's the output:

```
starting opt
==============
parameter_values=[56.22086187592684, 'medium']
new row ['experiment-1', 56.22086187592684, 2.0, 1]
wait
parsing /Users/kwaneu/research/zhulong/data/E-Z 2021-09-30 17-18-47/001-1-0312650-0713-SM.D/
plateaued=True
objective_function_value=0.0
new row ['experiment-1', 56.22086187592684, 2.0, 0]
existing
     identifier  temperature  HCl  sampling
0  experiment-1    56.220862  2.0         1
-
==============
parameter_values=[56.98726334005265, 'medium']
new row ['experiment-2', 56.98726334005265, 2.0, 1]
existing
     identifier  temperature  HCl  sampling
0  experiment-1    56.220862  2.0         0
-
wait
parsing /Users/kwaneu/research/zhulong/data/E-Z 2021-09-30 17-18-47/002-2-0312650-0713-BrPR.D/
plateaued=True
objective_function_value=1734.63611
new row ['experiment-2', 56.98726334005265, 2.0, 0]
existing
     identifier  temperature  HCl  sampling
0  experiment-1    56.220862  2.0         0
1  experiment-2    56.987263  2.0         1
-
==============
parameter_values=[6.52613148585219, 'low']
new row ['experiment-3', 6.52613148585219, 1.0, 1]
existing
     identifier  temperature  HCl  sampling
0  experiment-1    56.220862  2.0         0
1  experiment-2    56.987263  2.0         0
-
wait
parsing /Users/kwaneu/research/zhulong/data/E-Z 2021-09-30 17-18-47/003-3-0312650-0713-IS.D/
plateaued=True
objective_function_value=0.0
new row ['experiment-3', 6.52613148585219, 1.0, 0]
existing
     identifier  temperature  HCl  sampling
0  experiment-1    56.220862  2.0         0
1  experiment-2    56.987263  2.0         0
2  experiment-3     6.526131  1.0         1
-
finished opt
```

### Optimizer Class

Here are some more details on how to implement an `Optimizer`:

- Extend the `optimizer.Optimizer` class.
- `experiment_identifier` (`str`): every Experiment is assigned a unique string like "experiment-1".  Sampled points from the same experiment correspond to the same experiment identifier.
- Register `chemstation.Compound` instances as needed with the `Optimizer` so that the `compute_objective_function_value` method can access them.  These are needed to determine which integrals are useful for calculating yields or selectivities.

Here are useful fields in Optimizer:

- `parameter_space` (`list(Parameter)`): the space to optimizer over
- `chemstation_history` (`dict`): `experiment_identifier` : `list(Report)` (all the integration reports for this experiment)
- `parameter_values_history` (`dict`): `experiment_identifier` : `list(float/str)` (stores the experimental parameters)
- `objective_function_history` (`dict`): `experiment_identifier` : `float` (all the objective function values so far)

Four abstract methods need implementation:

- **`next_point(self)**
    Creates the next point in the parameter space as a list.  List members should be floats for `ContinuousParameter`s and strings for `CategoricalParameter`s.
    The list should be ordered the same as `optimizer.parameter_space`.
- **`compute_objective_function_value(self, experiment_identifier)`**
    Returns the objective function value for a particular set of experimental parameters.  It might be helpful to use `self.chemstation_history` to compute this.
- **`is_finished(self)`**
    Returns True if the optimization is finished.
- **`check_for_plateau(self, experiment_identifier)`**
    Returns True if we should stop sampling for this set of experimental parameters.  It might be helpful to use `self.chemstation_history` to compute this.

See `DummyOptimizer` in `optimizer.py` for an example.

### Base Classes

Here's some code that demonstrates some of the base classes in `zhulong`.  Let's define one continuous and one categorical parameter:

```
temperature_template = ContinuousParameter(name="temperature", min_value=-10.0, max_value=100.0)
additive_HCl_template = CategoricalParameter(name="HCl", allowed_levels={"low":1.0, "medium":2.0, "high":3.0})
```

Thus, temperature can be between -10 and 100 and HCl can be low (1.0), medium (2.0), or high (3.0).  The numbers will be written into the ChemSpeed `.csv` file directly, so whatever units are appropriate there are appropriate here.

Note that when the value of the HCl parameter is "low" the value 1.0 will be written to the file.

Let's come up with a few experiments to run:

```
conditions_list = [ [ 0.0, "low" ], [ 20.0, "low"], [ 0.0, "high" ], [ 20.0, "high" ] ]
parameter_templates = [ temperature_template, additive_HCl_template ]
experiments = []
for i, conditions in enumerate(conditions_list):
    identifier = f"experiment-{i}"
    parameters = [ p.copy() for p in parameter_templates ]
    for p,c in zip(parameters,conditions):
        p.value = c
    experiment = Experiment(identifier=identifier, parameters=parameters, sampling=True)
    experiments.append(experiment)
```

This tries a 2x2 factorial grid of 0.0 and 20.0 (temperature) and low and high (HCl).  Each experiment is assigned a corresponding `Experiment` object.  Note that every experiment is assigned a name ("experiment-1", "experiment-2", etc.).  Also note that each parameter is deep copied from a template, so that mutating the parameter only changes its value for its assigned Experiment.  However, there should be no reason to change Parameters once created.  Fresh ones should be made for new Experiments.

To write the Experiments to the `.csv` file:

```
# write to file
for experiment in experiments:
    experiment.update_csv(filename)
```

To update an experiment, use the `experiment.update_csv(filename)` method.  This might be used to stop sampling:

```
# update an experiment
experiments[1].sampling=False
experiments[1].update_csv(filename)
```

### Authors

Melodie Christensen, Eugene Kwan, Yuting Xu
