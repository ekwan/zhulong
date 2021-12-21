from chemstation import Peak, AgilentReport

# represents a reagent, either neat or in solution
class Reagent():
    # name (str) : name of reagent
    # abbreviation (str) : abbreviation for reagent
    # min_volume (int) : minimum volume in uL
    # max_volume (int) : maximum volume in uL
    def __init__(self, name, abbreviation, min_volume=None, max_volume=None):
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

# represents the space of parameters to explore in an optimization
class ParameterSpace():
    def __init__(self, starting_material, reagents, solvents, additives, light_stages, total_volume):
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

        assert isinstance(light_stages,int)
        assert light_stages > 0
        self.light_stages = light_stages

        assert isinstance(total_volume, int)
        assert total_volume > 0
        self.total_volume = total_volume

        self.experiment_count = 0

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
                       light_stage):
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
        assert starting_material_volume > 0
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
        assert reagent_volume > 0
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
            assert additive_volume > 0
            self.additive_volume = additive_volume

        self.solvent_volume = self.parameter_space.total_volume - starting_material_volume - reagent_volume - additive_volume
        assert self.solvent_volume >= 0

        if isinstance(light_stage, int):
            assert 1 <= light_stage <= parameter_space.light_stages
            self.light_stage = light_stage
        else:
            raise ValueError("unexpected light_stage")

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

