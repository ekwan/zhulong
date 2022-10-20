import pandas as pd
from pandas import DataFrame, Series
import numpy as np
import os
from dataclasses import dataclass
from glob import glob
from typing import List
import time

POLLING_INTERVAL = 5

# represents a ChemStation integration report


class AgilentReport():
    # fields:
    #
    # filename - which report file was parsed to generate the report
    # df       - DataFrame with peak centers and areas
    def __init__(self, folder):
        # look for the report file
        filename = f"{folder}/Report.TXT"

        while not os.path.exists(filename):
            print(f"Looking for {filename}...", end="\r")
            time.sleep(POLLING_INTERVAL)
        time.sleep(POLLING_INTERVAL)
        print(f"\nFound {filename}...parsing...")
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
            #         [min]        [min]   [mAU*s]     [mAU]        %
            #            1         2         3         4         5         6
            #  0123456789012345678901234567890123456789012345678901234567890
            # "----|-------|----|-------|----------|----------|--------|"
            # "   1   0.389 VB    0.0281   34.52074   14.83516   0.9311"
            if in_block:
                try:
                    retention_time = float(line[5:12])
                    area = float(line[48:56])
                except BaseException:
                    raise ValueError(
                        f"Error parsing integrals in {filename} for this line:\n{line}")
                row = [retention_time, area]
                rows.append(row)

            # go to next line
        # setup compounds to monitor by HPLC
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
    # if no match is found, return 0
    def get_integration(self, compound):
        assert compound.min_retention_time <= compound.max_retention_time, f"check retention times for compound {compound.name}"
        query_df = self.df.query(
            f"{compound.min_retention_time} <= retention_time <= {compound.max_retention_time}")
        if len(query_df) == 0:
            return 0.0
        return query_df.area.max()


@dataclass
# represents a compound of interest
class Peak:
    name: str
    min_retention_time: float
    max_retention_time: float

# create a method that computes yield from HPLC peaks
# peak_of_interest: product peak
# internal_standard_peak: if None, return peak_of_interest integral * response_factor
# if Peak, return response_factor * peak integral / IS integral


def get_yield_function(
        peak_of_interest,
        internal_standard_peak=None,
        response_factor=1.0):
    assert isinstance(peak_of_interest, Peak)
    if internal_standard_peak is not None:
        assert isinstance(internal_standard_peak, Peak)
    assert isinstance(response_factor, (int, float))
    assert response_factor > 0

    def yield_function(dot_D_folder):
        report = AgilentReport(dot_D_folder)
        peak_area = report.get_integration(peak_of_interest)

        if internal_standard_peak is not None:
            internal_standard_area = report.get_integration(
                internal_standard_peak)
            if internal_standard_area == 0:
                print(
                    f"Warning, internal standard not detected in {dot_D_folder}!")
                return 0.0
            chemical_yield = response_factor * peak_area / internal_standard_area
            return chemical_yield
        else:
            return response_factor * peak_area
    return yield_function


# for testing
if __name__ == '__main__':
    starting_material_peak = Peak(
        name="starting material",
        min_retention_time=1.28,
        max_retention_time=1.48)
    product_peak = Peak(
        name="bromo product",
        min_retention_time=1.53,
        max_retention_time=1.73)
    internal_standard_peak = Peak(
        name="internal standard",
        min_retention_time=2.05,
        max_retention_time=2.25)
    yield_function = get_yield_function(product_peak, response_factor=1.0)
    chemical_yield = yield_function(
        "../data/E-Z 2021-09-30 17-18-47/002-2-0312650-0713-BrPR.D")
    print(chemical_yield)
