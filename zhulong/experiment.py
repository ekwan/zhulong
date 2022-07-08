import json
from chemstation import Peak, AgilentReport

# represents a reagent, either neat or in solution
class Reagent():
    # name (str) : name of reagent
    # abbreviation (str) : abbreviation for reagent
    # min_volume (int) : minimum volume in uL
    # max_volume (int) : maximum volume in uL
    # concentration (float) : in mol/L (optional)
    def __init__(self, name, abbreviation, min_volume=None, max_volume=None, concentration=None):
        assert isinstance(name, str)
        assert len(name) > 0
        self.name = name
        assert isinstance(abbreviation, str)
        assert len(abbreviation) > 0
        self.abbreviation = abbreviation
        assert isinstance(min_volume, int)
        assert isinstance(max_volume, int)
        assert min_volume > 0
        assert max_volume >= min_volume
        self.min_volume = min_volume
        self.max_volume = max_volume
        if concentration is not None:
            assert isinstance(concentration, float)
            assert concentration > 0
        self.concentration = concentration

# represents the space of parameters to explore in an optimization
class ParameterSpace():
    def __init__(self, starting_material,
                       reagents,
                       solvents,
                       additives,
                       light_stages,
                       total_volume,
                       min_temperature,
                       max_temperature,
                       temperature_step_size):
        assert isinstance(starting_material, Reagent)
        assert starting_material.min_volume == starting_material.max_volume
        self.starting_material = starting_material

        assert isinstance(reagents, list)
        assert len(reagents) > 0
        for r in reagents:
            assert isinstance(r, Reagent)
        assert len(reagents) == len(set(reagents)), "duplicate reagent?"
        self.reagents = reagents

        assert isinstance(solvents, list)
        assert len(solvents) > 0
        for s in solvents:
            assert isinstance(s, str)
            assert len(s) > 0
        assert len(solvents) == len(set(solvents)), "duplicate solvent?"
        self.solvents = solvents

        assert isinstance(additives, list)
        assert len(additives) > 0
        for a in additives:
            assert isinstance(a, Reagent)
        assert len(additives) == len(set(additives)), "duplicate additive?"
        self.additives = additives

        assert isinstance(light_stages, int)
        assert light_stages > 0
        self.light_stages = light_stages

        assert isinstance(total_volume, int)
        assert total_volume > 0
        self.total_volume = total_volume

        assert isinstance(min_temperature, (int,float))
        assert isinstance(max_temperature, (int,float))
        assert isinstance(temperature_step_size, (int,float))
        assert temperature_step_size > 0
        assert max_temperature >= min_temperature
        self.min_temperature = int(min_temperature)
        self.max_temperature = int(max_temperature)
        self.temperature_step_size = int(temperature_step_size)

        self.experiment_count = 0

    def get_bounds_dict(self):
        bounds_dict = {}

        reagent = self.reagents[0]
        min_equivalents, max_equivalents, step_size = get_bounds(self.starting_material, reagent)
        bounds_dict["reagents"] = {
            #"name" : reagent.name,
            "names" : [ reagent.name for reagent in self.reagents ],# YX: make parameter_space more complete???
            "min_equivalents" : min_equivalents,
            "max_equivalents" : max_equivalents,
            "step_size" : step_size,
        }

        # assumes all additives have the same volume limits
        additive = self.additives[0]
        min_equivalents, max_equivalents, step_size = get_bounds(self.starting_material, additive)
        bounds_dict["additives"] = {
            "names" : [ additive.name for additive in self.additives ],
            "min_mole_percent" : min_equivalents*100,
            "max_mole_percent" : max_equivalents*100,
            "step_size" : step_size*100,
        }

        bounds_dict["temperature"] = {
            "min" : self.min_temperature,
            "max" : self.max_temperature,
            "step_size" : self.temperature_step_size,
        }

        bounds_dict["solvents"] = {
            "names" : self.solvents,
        }

        bounds_dict["light stages"] = {
            "names" : list(range(1,self.light_stages+1)),
        }

        return bounds_dict

# represents a single experiment inside a ParameterSpace
class Experiment():
    # parameter_space : the parameter space this experiment should correspond to
    # solvent (str) : as a string or list index for parameter_space.solvents
    # temperature : int, in celsius
    # starting_material_volume : int, in uL
    # reagent (Reagent) : as a string, list index for parameter_space.reagents, or Reagent
    # reagent_volume : int, in uL
    # additive (Reagent) : None (because additive is optional),
    #                      string, list index for parameter_space.additives, or Reagent
    # additive_volume : int, in uL
    # light_stage : int

    # additional fields that will be generated:
    # identifier : int, unique number for this experiment
    # starting_material (Reagent) : points to parameter_space
    # solvent_volume : int, in uL, makeup solvent will be added to satisfy parameter_space.total_volume
    def __init__(self, parameter_space,
                       solvent, temperature,
                       starting_material_volume,
                       reagent, reagent_volume,
                       additive, additive_volume,
                       light_stage, factory=False):
        assert isinstance(parameter_space, ParameterSpace)
        self.parameter_space = parameter_space

        if isinstance(solvent, str):
            self.solvent = None
            for s in parameter_space.solvents:
                if solvent == s:
                    self.solvent = s
                    break
            assert self.solvent is not None, f"couldn't find solvent {solvent}"
        elif isinstance(solvent, int):
            assert 0 <= solvent <= len(parameter_space.solvents)
            self.solvent = parameter_space.solvents[solvent]
        else:
            raise ValueError(f"unexpected solvent type: {type(solvent)}")

        assert isinstance(temperature, int)
        self.temperature = temperature

        self.parameter_space.experiment_count += 1
        self.identifier = self.parameter_space.experiment_count

        self.starting_material = self.parameter_space.starting_material
        assert isinstance(starting_material_volume, int)
        assert starting_material_volume == parameter_space.starting_material.min_volume, "check starting material volume"
        self.starting_material_volume = starting_material_volume

        if isinstance(reagent, str):
            self.reagent = None
            for r in parameter_space.reagents:
                if reagent == r.name:
                    self.reagent = r
                    break
            assert self.reagent is not None, f"couldn't find reagent {reagent}"
        elif isinstance(reagent, int):
            assert 0 <= reagent < len(parameter_space.reagents)
            self.reagent = parameter_space.reagents[reagent]
        elif isinstance(reagent, Reagent):
            assert reagent in parameter_space.reagents
            self.reagent = reagent
        else:
            raise ValueError(f"unexpected reagent type: {type(reagent)}")

        assert isinstance(reagent_volume, int)
        if not factory:
            assert self.reagent.min_volume <= reagent_volume <= self.reagent.max_volume,\
               f"reagent volume of {reagent_volume} is outside the range {self.reagent.min_volume}-{self.reagent.max_volume}"
        self.reagent_volume = reagent_volume

        if additive is None:
            self.additive = None
            self.additive_volume = 0
        else:
            if isinstance(additive, str):
                self.additive = None
                for r in parameter_space.additives:
                    if additive == r.name:
                        self.additive = r
                        break
                assert self.additive is not None, f"couldn't find additive {additive}"
            elif isinstance(additive, int):
                assert 0 <= additive < len(parameter_space.additives)
                self.additive = parameter_space.additives[additive]
            elif isinstance(additive, Reagent):
                assert additive in parameter_space.additives
                self.additive = additive
            else:
                raise ValueError("unexpected additive")

            assert isinstance(additive_volume, int)
            if not factory:
                assert self.additive.min_volume <= additive_volume <= self.additive.max_volume
            self.additive_volume = additive_volume

        total_volume = starting_material_volume + reagent_volume + additive_volume
        self.solvent_volume = self.parameter_space.total_volume - total_volume
        assert self.solvent_volume >= 0, f"total volume for this experiment is {total_volume}, which exceeds the max of {self.parameter_space.total_volume}"

        if isinstance(light_stage, int):
            assert 1 <= light_stage <= parameter_space.light_stages
            self.light_stage = light_stage
        else:
            raise ValueError("unexpected light_stage")

    def __str__(self):
        return_string  = f"solv={self.solvent} ({self.solvent_volume} uL) temp={self.temperature}°C "
        return_string += f"SM = {self.starting_material_volume} uL "
        return_string += f"reagent={self.reagent.abbreviation} (vol={self.reagent_volume} uL"
        if hasattr(self, "reagent_equivalents"):
            return_string += f", {self.reagent_equivalents:.2f} equiv"
        return_string += f") "
        return_string += f"additive={self.additive.abbreviation} (vol={self.additive_volume} uL"
        if hasattr(self, "additive_mole_percent"):
            return_string += f", {self.additive_mole_percent:.0f} mol%"
        return_string += f") "
        return_string += f"light={self.light_stage}"
        return return_string

    def get_dict(self):
        json_dict = {
            "solvent" : self.solvent,
            "temperature" : self.temperature,
            "reagent" : self.reagent.name,
            "reagent_volume" : self.reagent_volume,
            "reagent_equivalents" : self.reagent_equivalents if hasattr(self, "reagent_equivalents") else "None",
            "additive" : self.additive.name,
            "additive_volume" : self.additive_volume,
            "additive_mole_percent" : self.additive_mole_percent if hasattr(self, "additive_mole_percent") else "None",
            "light_stage" : self.light_stage,
            "history_times" : str(self.history["times"]),
            "history_values" : str(self.history["values"]),
        }
        return json_dict

    def get_json_string(self):
        return json.dumps(self.get_dict(), indent=2)

    # convenience factory method to create experiments using molar equivalents and mole percent instead of volumes
    @staticmethod
    def create(parameter_space, solvent, temperature, starting_material_volume, reagent,
               reagent_equivalents, additive, additive_mole_percent, light_stage):
        # make an experiment with minimum volumes for reagent and additive then adjust the volumes
        # based on the requested equivalents
        experiment = Experiment(parameter_space, solvent, temperature, starting_material_volume,
                                reagent, 1, additive, 1, light_stage, factory=True)

        # check that we have enough information
        assert parameter_space.starting_material.concentration is not None, "must specify starting material concentration"
        assert experiment.reagent.concentration is not None, "must specify reagent concentration"
        assert experiment.additive.concentration is not None, "must specify the additive concentration"

        # sanity checks
        assert reagent_equivalents >= 0
        assert additive_mole_percent >= 0

        # uL * mol / L = umol
        starting_material_concentration = parameter_space.starting_material.concentration
        starting_material_micromoles = starting_material_volume * starting_material_concentration

        # calculate reagent volume
        reagent_concentration = experiment.reagent.concentration
        reagent_volume_uL = int(starting_material_micromoles * reagent_equivalents / reagent_concentration)
        if reagent_volume_uL < experiment.reagent.min_volume or reagent_volume_uL > experiment.reagent.max_volume:
            raise ValueError("requested reagent equivalents translates to a volume that is outside the allowed range")
        experiment.reagent_volume = reagent_volume_uL
        experiment.reagent_equivalents = reagent_equivalents

        # calculate additive volume
        additive_concentration = experiment.additive.concentration
        additive_equivalents = additive_mole_percent / 100
        additive_volume_uL = int(starting_material_micromoles * additive_equivalents / additive_concentration)
        if additive_volume_uL < experiment.additive.min_volume or additive_volume_uL > experiment.additive.max_volume:
            raise ValueError("requested additive mol% translates to a volume that is outside the allowed range")
        experiment.additive_volume = additive_volume_uL
        experiment.additive_mole_percent = additive_mole_percent

        total_volume = starting_material_volume + experiment.reagent_volume + experiment.additive_volume
        experiment.solvent_volume = parameter_space.total_volume - total_volume
        assert experiment.solvent_volume >= 0, f"total volume for this experiment is {total_volume}, which exceeds the max of {parameter_space.total_volume}"

        # return result
        return experiment

    # get new row for chemspeed csv file
    # the order of columns is defined here
    #
    # sampling (int) : whether to keep sampling this well
    # returns: column_names (list of strings),
    #          values (list of integers)
    def get_row(self, sampling=1):
        assert sampling in [0,1]

        # initialize lists
        column_names = ["Experiment_ID", "Rxn_temp"]
        values = [self.identifier, self.temperature]
        solvents = self.parameter_space.solvents

        # starting material
        column_names.extend([f"{self.starting_material.abbreviation}_{solvent}" for solvent in solvents])
        values.extend([ 0 for solvent in solvents ])
        solvent_index = 2+solvents.index(self.solvent)
        values[solvent_index] = self.starting_material_volume
        assert len(column_names) == len(values), f"len(column_names)={len(column_names)} but len(values)={len(values)}"

        # reagent
        reagents = self.parameter_space.reagents
        for reagent in reagents:
            for solvent in solvents:
                name = f"{reagent.abbreviation}_{solvent}"
                column_names.append(name)
                if reagent == self.reagent and solvent == self.solvent:
                    value = self.reagent_volume
                else:
                    value = 0
                values.append(value)
        assert len(column_names) == len(values), f"len(column_names)={len(column_names)} but len(values)={len(values)}"

        # additive
        additives = self.parameter_space.additives
        for additive in additives:
            for solvent in solvents:
                name = f"{additive.abbreviation}_{solvent}"
                column_names.append(name)
                if additive == self.additive and solvent == self.solvent:
                    value = self.additive_volume
                else:
                    value = 0
                values.append(value)
        assert len(column_names) == len(values), f"len(column_names)={len(column_names)} but len(values)={len(values)}"

        # makeup solvent
        for solvent in solvents:
            column_names.append(solvent)
            if solvent == self.solvent:
                value = self.solvent_volume
            else:
                value = 0
            values.append(value)
        assert len(column_names) == len(values), f"len(column_names)={len(column_names)} but len(values)={len(values)}"

        # light stage
        column_names.append("Stage")
        values.append(self.light_stage)

        # sampling
        column_names.append("Sampling")
        values.append(sampling)

        # return result 
        assert len(column_names) == len(values), f"len(column_names)={len(column_names)} but len(values)={len(values)}"
        return column_names, values

def volume_to_equivalents(starting_material, reagent, reagent_volume):
    assert isinstance(starting_material, Reagent)
    assert starting_material.name == "starting material"
    assert isinstance(reagent, Reagent)
    assert hasattr(starting_material, "concentration")
    assert hasattr(reagent, "concentration")
    starting_material_moles = starting_material.concentration * starting_material.min_volume
    reagent_moles = reagent.concentration * reagent_volume
    return reagent_moles/starting_material_moles

def get_bounds(starting_material, reagent):
    min_volume, max_volume = reagent.min_volume, reagent.max_volume
    min_equivalents = volume_to_equivalents(starting_material, reagent, min_volume)
    max_equivalents = volume_to_equivalents(starting_material, reagent, max_volume)
    step_size = volume_to_equivalents(starting_material, reagent, 1)
    return min_equivalents, max_equivalents, step_size
