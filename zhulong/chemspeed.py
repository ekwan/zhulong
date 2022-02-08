import os
from glob import glob
from time import sleep
import pandas as pd
from pandas import DataFrame

from experiment import Experiment

# represents the ChemSpeed robot during an optimization run
# chemspeed_csv_filename (str) : where to write the volumes to
# chemstation_folder (str): what folder should be monitored for new .D folders
# overwrite_existing_chemspeed_csv (bool) : whether to overwrite chemspeed_csv_filename if it already exists
# ignore_existing_chemstation_folders (bool) : whether to ignore any existing
#                                              chemstation .D folders and only parse new ones
# yield_function (function) : parses the specified .D folder and returns a yield
# plateau_function (function) : receives a list of yields and returns True if we are plateaued and should stop sampling
#                               (if None, only one sample will be drawn) 
# polling_interval (int) : how often in seconds to check the chemstation folder for .D files
#
# additional fields:
# chemstation_parsed_folders (list): names of .D directories that have already been parsed
class ChemSpeed():
    def __init__(self, chemspeed_csv_filename,
                       chemstation_folder,
                       yield_function,
                       plateau_function=None,
                       overwrite_existing_chemspeed_csv=False,
                       ignore_existing_chemstation_folders=True,
                       polling_interval=1):
        assert isinstance(chemspeed_csv_filename, str)
        assert len(chemspeed_csv_filename) > 0
        assert isinstance(ignore_existing_chemstation_folders, bool)
        if ignore_existing_chemstation_folders and os.path.exists(chemspeed_csv_filename):
            raise ValueError("chemstation csv already exists (set overwrite_existing_chemspeed_csv=True to overwrite)")
        if os.path.exists(chemspeed_csv_filename):
            os.remove(chemspeed_csv_filename)
        self.chemspeed_csv_filename = chemspeed_csv_filename
        self.ignore_existing_chemstation_folders = ignore_existing_chemstation_folders

        assert isinstance(chemstation_folder, str)
        assert os.path.isdir(chemstation_folder), "chemstation folder does not exist"
        self.chemstation_folder = chemstation_folder

        assert isinstance(ignore_existing_chemstation_folders, bool)
        self.ignore_existing_chemstation_folders = ignore_existing_chemstation_folders
        self.chemstation_parsed_folders = []
        if ignore_existing_chemstation_folders:
            for directory in glob(f"{chemstation_folder}/*.D"):
                if os.path.isdir(directory):
                    self.chemstation_parsed_folders.append(directory)

        assert callable(yield_function)
        self.yield_function = yield_function

        if plateau_function is not None:
            assert callable(plateau_function)
        else:
            def plateau_function(history):
                return True
        self.plateau_function = plateau_function

        assert isinstance(polling_interval, int)
        assert polling_interval > 0
        self.polling_interval = polling_interval


    # runs one experiment and blocks until it is finished
    # experiment (Experiment) : the experiment to run
    def run_experiment(self, experiment):
        assert isinstance(experiment, Experiment)

        # write row with sampling turned on
        headings, values = experiment.get_row(sampling=1)
        self.append_to_chemspeed_csv(headings, values)

        # loop until plateau_function reports that we should stop
        experiment.history = {
                "times" : [],   # in epoch seconds
                "values" : [],  # objective function values (yields)
            }
        times = experiment.history["times"]
        values = experiment.history["values"]
        while True:
            # wait
            print("waiting")
            sleep(self.polling_interval)
            print("checking for new results")

            # see if there are new results from this run
            # if there is more than one match, the oldest file that has not been
            # parsed yet will be used
            directories = glob(f"{self.chemstation_folder}/*.D")
            new_directory_found = False
            for directory in sorted(directories, key=os.path.getmtime):
                if directory in self.chemstation_parsed_folders or not os.path.isdir(directory):
                    continue

                # found an unparsed directory, so parse it and mark it as parsed
                self.chemstation_parsed_folders.append(directory)
                new_directory_found = True
                print(f"found folder {directory}")
                folder_time = os.path.getmtime(directory)
                break
            if not new_directory_found:
                print("didn't find any new results")
                continue

            # parse and store the data from ChemStation
            print(f"parsing {directory}")
            chemical_yield = self.yield_function(directory)
            times.append(folder_time)
            values.append(chemical_yield)
            print(f"yield is {chemical_yield} (recorded at {folder_time} s since the epoch)")

            # see if the experiment has plateaued
            plateaued = self.plateau_function(experiment.history)
            print(f"{plateaued=}")
            if plateaued:
                # write row with sampling turned off
                headings, values = experiment.get_row(sampling=0)
                self.append_to_chemspeed_csv(headings, values)
                return


    # appends to csv file, creating the csv file if it does not already exist
    def append_to_chemspeed_csv(self, headings, values):
        if not hasattr(self, "_chemspeed_df"):
            # create new DataFrame and store it for future experiments
            df = DataFrame([values], columns=headings)
            self._chemspeed_df = df
        else:
            # note: not checking for duplicate experiment IDs
            df = self._chemspeed_df
            this_id = int(values[0])
            last_id = int(df.iloc[-1,0])
            #print(this_id, type(this_id), last_id, type(last_id), this_id==last_id)
            if this_id == last_id:
                # replace last row
                df.loc[len(df)-1] = values
            else:
                # append new row
                df.loc[len(df)] = values
            #print(df)

            assert list(df.columns) == headings
        df.to_csv(self.chemspeed_csv_filename, index=False)
