import pandas as pd
from pandas import DataFrame, Series
import numpy as np
import os
from dataclasses import dataclass
from glob import glob

# represents a ChemStation integration report
class Report():
    # fields:
    #
    # filename - which report file was parsed to generate the report
    # df       - DataFrame with peak centers and areas
    def __init__(self, folder):
        # look for the report file
        filename = f"{folder}/Report.TXT"
        if not os.path.exists(filename):
            raise ValueError(f"Error: can't find report file {filename}!")
        self.filename = filename

        # read the report file into the lines array
        with open(filename, encoding="utf-16") as f:
            lines = f.readlines()

        # keep track of where we are in the file
        i = 0
        in_block = False

        # store retention times and areas in this list
        rows = []

        # read file line by line
        while i < len(lines):
            line = lines[i]

            # search for the block that has the integrations
            if line.startswith("Signal 1: DAD1"):
                in_block = True
                i += 5
                continue
            elif line.startswith("Totals : "):
                break

            # if we're in the integration block, parse the
            # fixed-width text
            # sample:
            #
            #   Peak RetTime Type  Width     Area      Height     Area  
            #   [min]        [min]   [mAU*s]     [mAU]        %
            #            1         2         3         4         5         6
            #  0123456789012345678901234567890123456789012345678901234567890
            # "----|-------|----|-------|----------|----------|--------|"
            # "   1   0.389 VB    0.0281   34.52074   14.83516   0.9311"
            if in_block:
                try:
                    retention_time = float(line[5:12])
                    area = float(line[26:36])
                except:
                    raise ValueError(f"Error parsing integrals in {filename} for this line:\n{line}")
                row = [retention_time, area]
                rows.append(row)

            # go to next line
            i += 1

        # didn't find the block
        if len(rows) == 0:
            raise ValueError(f"Didn't find integrations in {filename}!")

        # create DataFrame
        df = DataFrame(rows)
        df.columns = ["retention_time", "area"]
        self.df = df

    # return the integration for the specified compound
    #
    # if multiple peaks are found within the retention time window for the compound,
    # return the largest area
    # if no match is found, return None
    def get_integration(self, compound):
        assert compound.min_retention_time <= compound.max_retention_time, f"check retention times for compound {compound.name}"
        query_df = self.df.query(f"{compound.min_retention_time} <= retention_time <= {compound.max_retention_time}")
        if len(query_df) == 0:
            return None
        return query_df.area.max()

@dataclass
# represents a compound of interest
class Compound:
    name: str
    min_retention_time: float
    max_retention_time: float

# for testing
if __name__ == '__main__':
    starting_material = Compound("starting material", 1.28, 1.48)
    product = Compound("bromo product", 1.53, 1.73)
    internal_standard = Compound("internal standard", 2.05, 2.25)

    directory = "/Users/kwaneu/research/zhulong/data/E-Z 2021-09-30 17-18-47"
    for sub_directory in sorted(glob(f"{directory}/*.D")):
        print(sub_directory)
        report = Report(sub_directory)
        print(report.df)
        print(f"starting material: {report.get_integration(starting_material)}")
        print(f"product:           {report.get_integration(product)}")
        print(f"internal standard: {report.get_integration(internal_standard)}")
        print()
