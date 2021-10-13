import pandas as pd
from pandas import DataFrame, Series
import numpy as np
import os
from glob import glob
from typing import Any, List, Callable
from dataclasses import dataclass
from copy import deepcopy
from experiment import Experiment, Parameter, ContinuousParameter, CategoricalParameter
from chemstation import Report, Compound
from time import sleep
import random

# Represents an optimizer that runs on the ChemSpeed platform.
@dataclass
class Optimizer():
    name : str                                 # name of the optimizer
    chemspeed_csv_filename : str               # where to send instructions to ChemSpeed (full path)
    chemstation_folder : str                   # watch this folder for new ChemStation files (/path/to/folder, no trailing slash)
    parameter_space : List[Parameter]          # the space of experimental parameters

    # auto-initialize fields:
    #  chemstation_parsed_folders : List[str]  # which ChemStation folders have been parsed already
    #  chemstation_history : dict              # { experiment identifier (str) : List[Report] }
    #  parameter_values_history : dict         # { experiment identifier (str) : List[float/str] }
    #  objective_function_history : dict       # { experiment identifier (str) : objective function value (float) }
    #  experiment_counter : int                # how many Experiments have been generated so far
    def __post_init__(self):
        self.chemstation_parsed_folders = []
        self.chemstation_history = {}
        self.parameter_values_history = {}
        self.objective_function_history = {}
        self.experiment_counter = 0

    # make a new Experiment object
    # values : the parameter values for the Experiment (list of float/str)
    # sampling : whether to continue sampling or not (bool)
    def _generate_new_experiment(self, values, sampling):
        assert len(values) == len(self.parameter_space), f"expected {len(self.parameter_space)} values, but got {len(values)} instead"
        parameters = [ p.copy() for p in self.parameter_space ]
        for p,v in zip(parameters, values):
            p.value = v
        self.experiment_counter += 1
        identifier = f"experiment-{self.experiment_counter}"
        experiment = Experiment(identifier=identifier, parameters=parameters, sampling=sampling)
        return experiment

    # run the specified experiment
    # values: the function values to try
    # polling interval: how often to check that the experiment is finished in seconds
    def _run_experiment(self, values, polling_interval=1):
        assert isinstance(values, list), f"expected list but got f{type(experiment)}"
        assert polling_interval >= 1, f"polling interval must be greater than 1 second"
        experiment = self._generate_new_experiment(values=values, sampling=True)
        self.parameter_values_history[experiment.identifier] = values
        experiment.update_csv(self.chemspeed_csv_filename)

        # loop until the experiment has plateaued
        while True:
            print("wait")
            # wait polling_interval seconds before checking status again
            sleep(polling_interval)

            # see if there are new results from this run
            # if there is more than one match, the oldest file that has not been
            # parsed yet will be used
            directories = glob(f"{self.chemstation_folder}/*.D/")
            new_directory_found = False
            for directory in sorted(directories, key=os.path.getmtime):
                if directory in self.chemstation_parsed_folders:
                    continue

                # found an unparsed directory, so parse it and mark it as parsed
                self.chemstation_parsed_folders.append(directory)
                new_directory_found = True
                break

            if not new_directory_found:
                continue

            # parse and store the data from ChemStation
            print(f"parsing {directory}")
            report = Report(directory)
            if experiment.identifier not in self.chemstation_history:
                self.chemstation_history[experiment.identifier] = []
            reports = self.chemstation_history[experiment.identifier]
            reports.append(report)

            # see if the experiment has plateaued
            plateaued = self.check_for_plateau(experiment.identifier)
            print(f"{plateaued=}")
            if plateaued:
                # compute the objective function value and store it
                objective_function_value = self.compute_objective_function_value(experiment.identifier)
                print(f"{objective_function_value=}")
                self.objective_function_history[experiment.identifier] = objective_function_value

                # tell the ChemSpeed to stop sampling
                experiment.sampling = False
                experiment.update_csv(self.chemspeed_csv_filename)
                break

    # start an optimization run
    def run(self):
        print("starting opt")
        while not self.is_finished():
            print("==============")
            parameter_values = self.next_point()
            print(f"{parameter_values=}")
            self._run_experiment(parameter_values)
        print("finished opt")

    ### abstract methods ###

    # get the next set of values to try
    # values can be floats (for continuous parameters) or strings (for categorical parameters)
    # returns list of values
    def next_point(self):
        pass

    # compute the objective function value
    # look at self.chemstation_history[experiment_identifier],
    # which stores all the ChemStation Reports for this experiment
    # experiment_identifier: str
    # returns: float
    def compute_objective_function_value(self, experiment_identifier):
        pass

    # determine whether this optimization is finished
    # returns: boolean
    def is_finished(self):
        pass

    # check if the specified experiment has plateaued
    # experiment_identifier (str): which experiment to check
    # returns: True if plateaued
    def check_for_plateau(self, experiment_identifier):
        pass

# Runs some random points for optimization.
class DummyOptimizer(Optimizer):
    # generate a random point within the parameter space
    def next_point(self):
        values = []
        for parameter in self.parameter_space:
            if isinstance(parameter, ContinuousParameter):
                min_value = parameter.min_value
                max_value = parameter.max_value
                random_value = random.uniform(min_value, max_value)
                values.append(random_value)
            elif isinstance(parameter, CategoricalParameter):
                allowed_levels = list(parameter.allowed_levels.keys())
                random_level = random.choice(allowed_levels)
                values.append(random_level)
            else:
                raise ValueError(f"unexpected parameter type: {type(parameter)}")
        return values

    # return a dummy objective function value
    def compute_objective_function_value(self, experiment_identifier):
        assert isinstance(experiment_identifier, str), f"expected str for experiment_identifier but got {type(experiment_identifier)} instead"
        reports = self.chemstation_history[experiment_identifier]  # all the reports for this id
        report = reports[-1] # look at last report only
        starting_material_integral = report.get_integration(self.starting_material)
        product_integral = report.get_integration(self.product)
        internal_standard_integral = report.get_integration(self.internal_standard)
        objective_function_value = product_integral / (internal_standard_integral+1)
        return objective_function_value

    # stop after three points
    def is_finished(self):
        return self.experiment_counter >= 3

    # assume "plateau" after one point
    def check_for_plateau(self, experiment_identifier):
        return True

# for testing
if __name__ == "__main__":
    # setup parameter space to run the optimization over
    param1 = ContinuousParameter(name="temperature", min_value=-10.0, max_value=100.0)
    param2 = CategoricalParameter(name="HCl", allowed_levels={"low":1.0, "medium":2.0, "high":3.0})
    parameter_space = [ param1, param2 ]

    # setup compounds to monitor by HPLC
    starting_material = Compound(name="starting material", min_retention_time=1.28, max_retention_time=1.48)
    product = Compound(name="bromo product", min_retention_time=1.53, max_retention_time=1.73)
    internal_standard = Compound(name="internal standard", min_retention_time=2.05, max_retention_time=2.25)

    # delete csv file if exists
    chemspeed_csv_filename = "temp.csv"
    if os.path.isfile(chemspeed_csv_filename):
        os.remove(chemspeed_csv_filename)

    # create the optimizer
    optimizer = DummyOptimizer(name="dummy optimizer", chemspeed_csv_filename=chemspeed_csv_filename,
                               chemstation_folder="/Users/kwaneu/research/zhulong/data/E-Z 2021-09-30 17-18-47",
                               parameter_space=parameter_space)

    # register compounds with optimizer
    optimizer.starting_material = starting_material
    optimizer.product = product
    optimizer.internal_standard = internal_standard

    # run the optimization
    optimizer.run()
