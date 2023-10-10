import csv
from argparse import ArgumentParser
from pathlib import Path
from typing import Dict, List, Tuple

from roifile import ImagejRoi
from scipy.spatial import KDTree

DATA_ROOT = Path() / "data"
OUTPUT_ROOT = Path() / "output"

def convert_and_join_rois(data_name: str, min_distance_between_rois: int):
    data_path = DATA_ROOT / data_name
    if not data_path.is_dir():
        raise FileNotFoundError(data_path)

    intensities: Dict[Tuple[int, int], int] = {}

    for data_file in data_path.glob("*.csv"):
        print("Loading ROIs from", data_file)
        with data_file.open() as csv_file:
            csv_reader = csv.DictReader(csv_file)
            for row in csv_reader:
                # Yes, we're throwing away the channel this came from.
                # This is okay -- we should have started from normalized images so intensities can be compared across channels.
                intensities[(int(row["X"]), int(row["Y"]))] = int(row["Max"])

    points = list(intensities.keys())
    tree = KDTree(points)

    print(f"Pruning ROIs that are less than {min_distance_between_rois} pixels away from each other")
    # Need to copy this since we will be modifying the dictionary
    for point, intensity in list(intensities.items()):
        neighbor_indices = tree.query_ball_point([point], r=min_distance_between_rois)
        for neighbor_index in neighbor_indices[0]:
            neighbor_point = points[neighbor_index]
            if point == neighbor_point:
                continue

            # Found a pair that's too close. Keep the larger one...
            neighbor_intensity = intensities.get(neighbor_point)
            if not neighbor_intensity:
                # ... unless it was already deleted...
                continue
            if intensity < neighbor_intensity:
                # ... by removing the smaller one from the `intensities` dictionary. We have to keep the `points` list intact because we are still relying on its indices.
                print(f"Remove {points.index(point)} Keep {neighbor_index}")
                intensities.pop(point, None)
            else:
                print(f"Remove {neighbor_index} Keep {points.index(point)}")
                intensities.pop(neighbor_point, None)

    output_file = OUTPUT_ROOT / f"{data_name}.zip"
    print("Saving ROIs to", output_file)
    for i, point in enumerate(intensities.keys()):
        roi = ImagejRoi.frompoints([point], name=str(i))
        roi.tofile(output_file, name=str(i))

parser = ArgumentParser()
parser.add_argument("name")
parser.add_argument("-r", type=int, default=2, help="Minimum distance between ROIs")

if __name__ == "__main__":
    args = parser.parse_args()
    convert_and_join_rois(args.name, args.r)
