import pandas as pd
from pandas import DataFrame, Series
import numpy as np
import os
from glob import glob
from typing import Any, List
from dataclasses import dataclass
from copy import deepcopy

# Represents an experimental parameter.
class Parameter():
    def __init__(self, name):
        self.name = name

    @property
    def value(self) -> Any:
        return self._value

    @value.setter
    def value(self, value: Any):
        assert self._is_valid(value), "trying to set invalid value"
        self._value = value

    def _is_valid(self, value):
        return True

    def copy(self):
        return deepcopy(self)

# Represents a continuous experimental parameter.
class ContinuousParameter(Parameter):
    def __init__(self, name : str, min_value : float, max_value : float):
        super().__init__(name)
        self.min_value = min_value
        self.max_value = max_value

    def __repr__(self):
        report_string = f"ContinuousParameter({self.name=}"
        if hasattr(self, "value"):
            report_string += f", {self.value=}"
        report_string += f", min={self.min_value}, max={self.max_value})"
        return report_string

    def _is_valid(self, value):
        is_valid = self.min_value <= value <= self.max_value
        assert is_valid, f"{value=} out of range {self.min_value} to {self.max_value}"
        return True

# Represents a categorical experimental parameter.
class CategoricalParameter(Parameter):
    # allowed_levels = { additive_level_name : float }
    def __init__(self, name : str, allowed_levels : dict):
        super().__init__(name)
        all_levels = set()
        for additive_level_name, additive_level in allowed_levels.items():
            assert isinstance(additive_level_name, str) and len(additive_level_name) > 0, f"invalid: {additive_level_name=}"
            assert isinstance(additive_level, float) and additive_level >= 0.0, f"invalid: {additive_level=}"
            assert additive_level not in all_levels, f"duplicate level name for {name} at level {additive_level}"
            all_levels.add(additive_level)
        self.allowed_levels = allowed_levels

    def __repr__(self):
        report_string = f"CategoricalParameter({self.name}"
        if hasattr(self, "value"):
            report_string += f"={self.value}"
        report_string += ", allowed_levels={"
        for i, (additive_level_name, additive_value) in enumerate(self.allowed_levels.items()):
            if i > 0:
                report_string += ", "
            report_string += f"{additive_level_name}={additive_value:.2f}"
        report_string += "})"
        return report_string

    # checks if the value is one of the allowed level names
    def _is_valid(self, value):
        assert isinstance(value, str), f"expected str but got {type(value)} for CategoricalParameter value"
        is_valid = value in self.allowed_levels
        assert is_valid, f"{value=} does not correspond to one of the level names in CategoricalParameter({self.name})"
        return is_valid

# Represents a single experiment
@dataclass
class Experiment():
    identifier : str
    parameters : List[Parameter]
    sampling : bool  # whether to continue sampling (True) or stop (False)

    def copy(self):
        return deepcopy(self)

    # Updates the given csv with the specified values for this experiment
    # If the entry already exists, it is updated; if it does not, it is appended.
    # Expected columns are ["identifier", "param1", "param2", ..., "paramN", "sampling"]
    def update_csv(self, filename):
        # make columns for csv
        columns = ["identifier"]
        columns.extend([ p.name for p in self.parameters ])
        columns.append("sampling")

        # create the new row
        new_row = [ self.identifier ]
        for p in self.parameters:
            # continuous parameter, so write the value directly
            if isinstance(p, ContinuousParameter):
                new_row.append(p.value)

            # categorical parameter, so write the numerical value
            elif isinstance(p, CategoricalParameter):
                numerical_value = p.allowed_levels[p.value]
                new_row.append(numerical_value)

            # other parameter types must be supported explicitly
            else:
                raise ValueError("unexpected parameter type: {type(parameter)}")

        sampling_value = 1 if self.sampling else 0  # write 1 to continue sampling and 0 to stop
        new_row.append(sampling_value)
        print("new row", new_row)

        # check if csv file already exists
        if os.path.isfile(filename):
            # file exists, so read it
            df = pd.read_csv(filename)
            print("existing")
            print(df)
            print("-")

            # ensure columns are as expected
            for c,c2 in zip(columns, df.columns):
                if c != c2:
                    print()
                    print(f"{c} does not match {c2}")
                    print()
                    lines = read_text_file(filename)
                    print(lines[0:2])
                    raise ValueError(f"Unexpected columns in {filename}.")
        else:
            # create the df
            df = DataFrame([new_row], columns=columns)
            df.to_csv(filename, index=False)
            return

        query_df = df.query(f"identifier == '{self.identifier}'")
        if len(query_df) == 0:
            new_row = Series(new_row, index = df.columns)
            df = df.append(new_row, ignore_index=True)
        elif len(query_df) == 1:
            df[df.identifier == self.identifier] = new_row
        else:
            print()
            print(query_df)
            raise ValueError(f"{self.identifier=} found multiple times in {filename}")
        df.to_csv(filename, index=False)

# get contents of text file for debugging
def read_text_file(filename):
    with open(filename, "r") as f:
        lines = f.readlines()
    return "".join(lines)

# for testing
if __name__ == "__main__":
    param1_template = ContinuousParameter(name="temperature", min_value=-10.0, max_value=100.0)
    print(param1_template)
    param1_template.value=0.0
    print(param1_template)
    try:
        param1_template.value=-20.0
    except Exception as e:
        print(e)
    print()
    param2_template = CategoricalParameter(name="HCl", allowed_levels={"low":1.0, "medium":2.0, "high":3.0})
    print(param2_template)
    param2_template.value="low"
    print(param2_template)
    try:
        param2_template.value="adsfadsf"
    except Exception as e:
        print(e)
    print()

    # come up with a few experiments to try
    conditions_list = [ [ 0.0, "low" ], [ 20.0, "low"], [ 0.0, "high" ], [ 20.0, "high" ] ]
    parameter_templates = [ param1_template, param2_template ]
    experiments = []
    for i, conditions in enumerate(conditions_list):
        identifier = f"experiment-{i}"
        parameters = [ p.copy() for p in parameter_templates ]
        for p,c in zip(parameters,conditions):
            p.value = c
        experiment = Experiment(identifier=identifier, parameters=parameters, sampling=True)
        experiments.append(experiment)

    # delete file if exists
    filename = "temp.csv"
    if os.path.isfile(filename):
        os.remove(filename)

    # write to file
    for experiment in experiments:
        experiment.update_csv(filename)
        print("final csv")
        print(read_text_file(filename))
        print("---")

    # update an experiment
    print("update")
    experiments[1].sampling=False
    experiments[1].update_csv(filename)
    print(read_text_file(filename))

