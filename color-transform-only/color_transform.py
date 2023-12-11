from argparse import ArgumentParser
from typing import Tuple
import re

import pandas as pd
import numpy as np
from sklearn import linear_model
from scipy.spatial import KDTree

DYE_BASES = ["G", "C", "A", "T"]
DICTIONARY_SPOTS = DYE_BASES + ["BG"]

FILENAME_FORMAT = re.compile(r'(\d+)_(\d+)_(C\d+).tif')
def parse_image_filename(filename: str) -> Tuple[int, int]:
    parsed = FILENAME_FORMAT.search(filename)
    # TODO: Convert to int?
    cycle = int(parsed.group(3).lstrip("C"))
    wavelength = int(parsed.group(2))
    return (cycle, wavelength)

def select_cycle(spots: pd.DataFrame, cycle: int) -> pd.DataFrame:
    # Manually select the wavelengths to ensure the order is consistent
    return spots[[(cycle, 645), (cycle, 590), (cycle, 525), (cycle, 445)]]

def deduplicate_spots(spots: pd.DataFrame, min_distance_between_rois: int) -> pd.DataFrame:
    tree = KDTree(spots["position"])
    spot_indices_to_remove = set()
    for first_index, second_index in tree.query_pairs(min_distance_between_rois):
        first = spots.iloc[first_index]
        second = spots.iloc[second_index]

        # First, ensure we keep the dictionary spots
        if first.name in DICTIONARY_SPOTS:
            spot_indices_to_remove.add(second_index)
        elif second.name in DICTIONARY_SPOTS:
            spot_indices_to_remove.add(first_index)
        else:
            # Otherwise, keep the brighter spot across all wavelengths for the first cycle
            first_signal = sum(select_cycle(first, 1))
            second_signal = sum(select_cycle(second, 1))
            if first_signal > second_signal:
                spot_indices_to_remove.add(second_index)
            else:
                # `else`: this case is selected if there is a tie
                spot_indices_to_remove.add(first_index)

    return spots.drop(spots.iloc[list(spot_indices_to_remove)].index)

def calculate_transformation(spots: pd.DataFrame) -> linear_model.LinearRegression:
    # TODO: Automatically find reasonable dictionary spots instead
    # TODO: Why do we care about background/BG? We never use it
    dictionary_spots = spots[spots.index.isin(DICTIONARY_SPOTS)]
    if len(dictionary_spots) != len(DICTIONARY_SPOTS):
        print(dictionary_spots)
        raise Exception(f"Dictionary missing from spots")

    dye_spot_to_index_map = {dye_spot: i for i, dye_spot in enumerate(DYE_BASES)}

    # Set up a linear regression to determine each dye's contribution to each channel.
    # X is the input data from the dictionary spots for the first cycle only, where the bases are known.
    X = select_cycle(dictionary_spots, 1)
    # Y is the identity matrix we are trying to transform into.
    Y = np.zeros((len(dictionary_spots), len(DYE_BASES)))
    for i, base_spotname in enumerate(dictionary_spots.index):
        if base_spotname != "BG":
            Y[i, dye_spot_to_index_map[base_spotname]] = 1

    transformation = linear_model.LinearRegression()
    return transformation.fit(X, Y)

def apply_transformation(spots: pd.DataFrame, transformation: linear_model.LinearRegression) -> pd.DataFrame:
    cycles = set([col[0] for col in spots.columns if isinstance(col[0], int)])
    transformed = spots.drop(columns=cycles)

    def transformation_impl(x):
        # Temporarily reshape to satisfy sklearn
        return transformation.predict(x.to_numpy().reshape(1, -1)).reshape(-1)

    for cycle in cycles:
        cycle_spots = select_cycle(spots, cycle)
        applied = cycle_spots.apply(transformation_impl, axis="columns", result_type="expand")
        applied.columns = pd.MultiIndex.from_tuples(zip([cycle]*len(DYE_BASES), DYE_BASES))

        transformed = pd.concat([transformed, applied], axis="columns")

    return transformed

def convert_to_color_transformed_spots(transformed: pd.DataFrame, output_path: str, reindex: bool = True) -> None:
    unique_spots = transformed.index.unique()
    spot_name_to_index = {name: index for index, name in enumerate(unique_spots)}
    # TODO: Make this a CSV instead
    print(spot_name_to_index)

    cycles = set([col[0] for col in transformed.columns if isinstance(col[0], int)])
    # Take only the cycle data...
    color_transformed_spots = transformed[list(cycles)]
    # ... make a column for cycles...
    color_transformed_spots = color_transformed_spots.stack(0, future_stack=True)
    # ... rename the columns...
    color_transformed_spots.index.rename(["spot_name", "cycle"], inplace=True)
    color_transformed_spots = color_transformed_spots.reset_index()
    # ... create the spot_index column in the right spot...
    color_transformed_spots["spot_index"] = color_transformed_spots["spot_name"].apply(spot_name_to_index.get)
    color_transformed_spots = color_transformed_spots.reindex(["spot_index", "spot_name", "cycle", "G", "C", "A", "T"], axis="columns")
    # ... and rename if requested.
    if reindex:
        color_transformed_spots["spot_name"] = color_transformed_spots["spot_index"]

    print(color_transformed_spots)
    color_transformed_spots.to_csv(output_path, index=False)

parser = ArgumentParser()
parser.add_argument("spots_path")
parser.add_argument("-o", default="color_transformed_spots.csv")
parser.add_argument("-r", type=int, default=4, help="Minimum distance between ROIs")

if __name__ == "__main__":
    args = parser.parse_args()

    # Load and format data
    spots = pd.read_csv(args.spots_path)
    spots.columns = pd.MultiIndex.from_tuples([("spot", "id"), ("position", "x"), ("position", "y"), *map(parse_image_filename, spots.columns[3:])])
    spots = spots.set_index(("spot", "id")).astype(np.uint32)

    # Run multiple times to ensure only one detected spot remains for each cluster
    if args.r != 0:
        for _ in range(3):
            spots = deduplicate_spots(spots, args.r)

    spots.to_csv("deduplicated_spots.csv")

    transformation = calculate_transformation(spots)
    transformed_spots = apply_transformation(spots, transformation)
    transformed_spots.to_csv("transformed_spots.csv")

    convert_to_color_transformed_spots(transformed_spots, args.o)
