# Zhulong

Autonomous optimization code.  This code provides code to run autonomous optimization with a ChemSpeed Robot and Agilent HPLC system.

### Example Code

Here's some code that demonstrates some of the base classes in `zhulong`.  Let's define one continuous and one categorical parameter:

```
temperature_template = ContinuousParameter(name="temperature", min_value=-10.0, max_value=100.0)
additive_HCl_template = CategoricalParameter(name="HCl", allowed_levels={"low":1.0, "medium":2.0, "high":3.0})
```

Thus, temperature can be between -10 and 100 and HCl can be low (1.0), medium (2.0), or high (3.0).  The numbers will be written into the ChemSpeed `.csv` file directly, so whatever units are appropriate there are appropriate here.

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
