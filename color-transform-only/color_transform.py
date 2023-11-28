from argparse import ArgumentParser
import math
import os
import re
from typing import Tuple
import pandas as pd
import numpy as np
from sklearn import linear_model
import cv2 as cv
import matplotlib.pyplot as plt
import matplotlib.colors
from scipy import ndimage
from scipy.spatial import KDTree
import roifile
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# TODO: Refactor users to use intensity values directly from the CSV
def OLD_get_cycle_files():
    pass

DYE_BASES = ["G", "C", "A", "T"]
DICTIONARY_SPOTS = DYE_BASES + ["BG"]

OLIGO_SEQUENCES = {
    '355': 'ACGTGACTAGTGCATCACGTGACTAGTGCATC',
    '357': 'ATGCAGTCGACGTACTATGCAGTCGACGTACT',
    '358': 'CGTATCGACTATGCAGCGTATCGACTATGCAG',
    '360': 'GACTCGATGCTCAGTAGACTCGATGCTCAGTA',
    '364': 'TCAGTACGATGACTGCTCAGTACGATGACTGC',
    '370': 'ACGTGACTAGTGCATCACGTGACTAGTGCATC',
    '372': 'ATGCAGTCGACGTACTATGCAGTCGACGTACT',
    '373': 'CGTATCGACTATGCAGCGTATCGACTATGCAG',
    '375': 'GACTCGATGCTCAGTAGACTCGATGCTCAGTA',
    '377': 'GTCAGCTACGACTGATGTCAGCTACGACTGAT',
    '379': 'TCAGTACGATGACTGCTCAGTACGATGACTGC',
    '574': 'GGGGGGGGGGTAAGAA',
    '575': 'AAAAAAAAAATAAGAA',
    '576': 'CCCCCCCCCCTAAGAA',
    '577': 'TTTTTTTTTTTAAGAA',
    '632': 'AAATGCAGTCGACGTACTATGCAGTC',
    '633': 'CCCGTATCGACTATGCAGCGTATCGA',
    '634': 'GGGACTCGATGCTCAGTAGACTCGAT',
    '635': 'TTTCAGTACGATGACTGCTCAGTACG',
    '648': 'TTTGCATTAAGAAATTAAAAAAGCTAAAAAAAAAA',
    '649': 'AAAGCATTAAGAAATTAAAAAAGCTAAAAAAAAAA',
    '650': 'GGGGCATTAAGAAATTAAAAAAGCTAAAAAAAAAA',
    '651': 'CCCGCATTAAGAAATTAAAAAAGCTAAAAAAAAAA',
    '657': 'GGGCATCTCGTATGCC',
    '662': 'ACTGATCTCGTATGCC',
    '663': 'GCTGATCTCGTATGCC',
    '664': 'CAGCATCTCGTATGCC',
    '665': 'TCTGATCTCGTATGCC',
}

DEFAULT_BASE_COLOR_MAP = {
    'A': 'green',
    'C': 'blue',
    'T': 'red',
    'G': 'black',
    'BG': 'lightgrey'
}


def OLD_get_roi_mask(shape, rois):
    pass



def plot_bars(
        df: pd.DataFrame,
        title: str
) -> go.Figure:

    spot_indizes = df.spot_index.unique()

    '''
    # Create figure with secondary x-axis
#    fig = go.Figure(layout=layout)
    layout = go.Layout(
        title="Basecalls, Spot " + spot_name,
        xaxis=XAxis(
            title="Cycles"
        ),
        xaxis2=XAxis(
            title="ACGT",
            overlaying='x',
            side='top',
        ),
        yaxis=dict(
            title="Y values"
        ),
    )
    '''

    cols = 4
    fig = make_subplots(
        rows=math.ceil(len(spot_indizes)/cols), cols=cols
    )

    for i, spot_index in enumerate(spot_indizes):

        r = (i // cols)+1
        c = (i % cols)+1

        df_spot = df.loc[(df['spot_index'] == spot_index)]
        spot_name = df_spot.spot_name.unique()[0]
        print(f"spot: {i}, idx: {spot_index}, name: {spot_name}, row={r}, col={c}")

        # Add traces
        for base_spot_name in DYE_BASES:
            fig.add_trace(
                # Scatter, Bar
                go.Bar(
                    x=df_spot['cycle'],
                    y=df_spot[base_spot_name],
                    name=base_spot_name,
                    marker_color=DEFAULT_BASE_COLOR_MAP[base_spot_name],
                    legendgroup=base_spot_name, showlegend=(i == 0)
                ),
                row=r, col=c
            )

        fig.add_trace(
            # Scatter, Bar
            go.Scatter(
                x=df_spot['cycle'],
                y=df_spot['G']/1000000+1,
                text=df_spot[DYE_BASES].idxmax(axis=1),  # column with the highest value
                marker_color="black",
                mode="text",
                textposition="top center",
                textfont_size=26,
                showlegend=False
            ),
            row=r, col=c
        )

        fig.update_xaxes(
            title_text=str(spot_index) + '  ' + spot_name + "  (" + OLIGO_SEQUENCES.get(spot_name, "")[:16] + ")",
            title_font={"size": 24},
            row=r, col=c
        )

        fig.update_yaxes(
            range=[-0.2, 1.2],
            row=r, col=c
        )

    fig.update_layout(
        height=3000,
        width=3000,
        title_text=title
    )

    fig.update_layout(
        legend=dict(title_font_family="Times New Roman",
                    font=dict(size=40)
                    )
    )

    return fig


def OLD_calculate_and_apply_transformation(
    OLD_spot_data_filename: str,
    OLD_roizipfilepath: str,
    OLD_input_directory_path: str,
    output_directory_path: str
) -> pd.DataFrame:
    
    channel_names = ['M445', 'M525', 'M590', 'M645']

    df = pd.read_csv(OLD_spot_data_filename)
    print(df)

    '''
        spot  pixel_i  timestamp_ms  cycle  R365  ...   G590   B590   R645   G645   B645
    0   BG          0         14241      1  4592  ...   5688   4736   5712   4984   4080
    1   BG          1         14241      1  4496  ...   5344   4784   6880   5136   4272
    70  SC         10         14241      1  5824  ...  65520  52528  65520  65520  35360
    71  SC         11         14241      1  4720  ...  65520  40736  65520  65520  29776
    '''

#    df = df[(df.spot != "B0")] # TODO
    # only look at these spots
#    df = df[df.spot.isin(["D1", "D2", "D3", "D4", "L0"])]

    df = df[df.spot_name.isin(DYE_BASES+['BG'])]

    n_features = len(channel_names)
    n_targets = len(DYE_BASES)

    unique_spot_names = df['spot_name'].unique()
    print("spots:", unique_spot_names)
    # make sure all DYE_BASES are in the spots_ist

    print("DYE_BASES:", DYE_BASES)

    dye_spot_to_index_map = {dye_spot: idx for idx, dye_spot in enumerate(DYE_BASES)}
    print(dye_spot_to_index_map)

    # intensity data
    X = df[channel_names].to_numpy()

    mono = 'M445' in df.columns

    # camera offset correction
    if not mono:
        offset = 4096
        BG_threshold = 25500
        SC_threshold = 64000*4
        X[X < offset] = offset
        X -= offset

    print("Datamatrix dimension:", X.shape)

    '''
    Generate Y
    1 0 0 0
    1 0 0 0
    0 1 0 0
    0 0 1 0
    0 0 1 0
    0 0 0 1
    0 0 0 1
    '''
    # (n_targets, n_features)
    # Y is a 2D matrix
    Y = np.zeros((len(df), n_targets))
    for i, base_spotname in enumerate(df['spot_name']):
        if base_spotname != "BG":
            Y[i, dye_spot_to_index_map[base_spotname]] = 1
#        else:
            #map to zero

    print(Y)

    model_ols = linear_model.LinearRegression()
    reg = model_ols.fit(X, Y)
    print(type(reg))
    print(reg.coef_)
    coef = model_ols.coef_
    intercept = model_ols.intercept_
    print('coef= ', coef)
    print('intercept= ', intercept)
#    np.set_printoptions(suppress=True)
    np.set_printoptions(threshold=np.inf)
    np.set_printoptions(linewidth=np.inf)
    print('coef= \n', coef*10000)
    print('intercept= ', intercept)


    # apply matrix to each cycle

    df_files = OLD_get_cycle_files(OLD_input_directory_path)
    print(df_files)

    rois = roifile.ImagejRoi.fromfile(OLD_roizipfilepath)
    print(rois)

    rows_list = []

    nb_cycles = max(df_files['cycle'])
    print("cycles:", nb_cycles)
#    nb_cycles = 8

    for cycle in range(1, nb_cycles+1):
        image_map = {}

        print("Apply transformation matrix on:")
        cyclefilenames = (df_files[df_files['cycle'] == cycle]).tail(5)
        print("cycle:", cycle, cyclefilenames.to_string())
        cycle_timestamp = cyclefilenames.iloc[0]['timestamp']

        for i, cyclefilename in cyclefilenames.iterrows():
            filenamepath = cyclefilename['filenamepath']
            wavelength = cyclefilename['wavelength']
            image = cv.imread(filenamepath, cv.IMREAD_UNCHANGED)  # 16bit data
            print("WL:", filenamepath, wavelength, image.shape)

            if mono:
                # no offset correction required, microscope camera is true 16 bit
                image_map['M' + str(wavelength)] = image
            else:
                image = image[:, :, ::-1]  # BGR to RGB
                # safe offset correction
                image[image < offset] = offset
                image -= offset

                image_map['R' + str(wavelength)] = image[:, :, 0]
                image_map['G' + str(wavelength)] = image[:, :, 1]
                image_map['B' + str(wavelength)] = image[:, :, 2]



        channels = [image_map[channel_name] for channel_name in channel_names]

        A = np.stack(channels, axis=2)
        (dim0, dim1, dim2) = A.shape
        print(f"Matrix A shape: {A.shape}")
        assert (dim2 == n_features)

        # apply transformation to each pixel, reshape temporarily
        a = reg.predict(A.reshape(dim0*dim1, n_features))
        # reshape back
        a = a.reshape(dim0, dim1, n_targets)

        assert (a.shape[0] == dim0)
        assert (a.shape[1] == dim1)
        assert (a.shape[2] == n_targets)
        print("Transformation applied, shape", a.shape, type(a))

        oligo_mask, nb_labels = OLD_get_roi_mask((dim0, dim1), rois)

#        nb_labels = len(rois)
        label_ids = np.arange(1, nb_labels + 1)  # range(1, nb_labels + 1)
        mean_list = []
        for i in range(n_targets):
            mean_list.append(ndimage.labeled_comprehension(a[:, :, i], oligo_mask, label_ids, np.mean, float, 0))

            img = a[:, :, i]
            cv.imwrite(os.path.join(output_directory_path, f"C{cycle:03d}_{DYE_BASES[i]}_{cycle_timestamp:09d}_gray.png"), (img+1)*100)
#            cv.imwrite(os.path.join(output_directory_path, f"C{cycle:03d}_{DYE_BASES[i]}_{cycle_timestamp:09d}_gray.tif"), img)


        for j, roi in enumerate(rois):
            if __debug__:
                print(roi.name, roi.top, roi.bottom, roi.left, roi.right, roi.roitype, roi.subtype, roi.options, roi.version, roi.props, roi.position)

            dict_entry = {
                'spot_index': j+1,
                'spot_name': roi.name,
                'cycle': cycle,
            }
            # base vector coefficients
            for i, base_spot_name in enumerate(DYE_BASES):
                dict_entry[base_spot_name] = mean_list[i][j]
            rows_list.append(dict_entry)


        # debug subplots
        if cycle == 1:

            fig, axs = plt.subplots(1, 5)

            for i in range(n_targets):
                img = a[:, :, i]
                print(f"min:  {img.min()}  , max: {img.max()}")

                cax_01 = axs[i].imshow(img, cmap='gray')
                fig.colorbar(cax_01, ax=axs[i])
                #        axs[i].xaxis.set_major_formatter(plt.NullFormatter())
                #        axs[i].yaxis.set_major_formatter(plt.NullFormatter())

                #    plt.show()

            mask = np.zeros(A.shape[:2], dtype=np.uint16)
            print(f"mask: {type(mask)}, {mask.dtype}, {mask.shape}")

            counts = [0] * n_targets
            for r in range(a.shape[0]):
                for c in range(a.shape[1]):
                    pixel = a[r, c]
                    label = np.array(pixel).argmax()
                    #            print(sum(A[r,c]), label)
                    #            if BG_threshold < sum(A[r,c]) and sum(A[r,c]) < SC_threshold:
                    mask[r, c] = label
                    #            else:
                    #                mask[r, c] = 5 # TODO

                    counts[label] += 1
            print("counts per base:", counts)

            colors = [
                DEFAULT_BASE_COLOR_MAP['BG'],
                DEFAULT_BASE_COLOR_MAP['G'],
                DEFAULT_BASE_COLOR_MAP['C'],
                DEFAULT_BASE_COLOR_MAP['A'],
                DEFAULT_BASE_COLOR_MAP['T'],
                'yellow', 'orange', 'magenta'
            ]
            scale = [0, 1, 2, 3, 4, 5, 6, 250]
            cmap = matplotlib.colors.ListedColormap(colors)
            norm = matplotlib.colors.BoundaryNorm(scale, len(colors))
            axs[4].imshow(mask, aspect="auto", cmap=cmap, norm=norm)

            #    axs[5].imshow(mask, aspect="auto", cmap=cmap, norm=norm, extent=[0, 400, 0, 300])

            #    x = np.random.normal(170, 10, 250)
            #    axs[5].hist(x)

            #    x = range(300)
            #    axs[5].plot(x, x, '--', linewidth=5, color='firebrick')

            #    plt.imshow(mask, aspect="auto", cmap=cmap, norm=norm)
            plt.show()

    # create final dataframe
    df_out = pd.DataFrame(rows_list)
    df_out.sort_values(by=['spot_index', 'spot_name', 'cycle'], inplace=True)
#           print(f"Writing {outputfilename}")
    df_out.to_csv(os.path.join(output_directory_path, "color_transformed_spots.csv"), index=False)
    print(df_out.to_string(index=False))

    return df_out

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

parser = ArgumentParser()
parser.add_argument("spots_path")
parser.add_argument("-r", type=int, default=2, help="Minimum distance between ROIs")

if __name__ == "__main__":
    args = parser.parse_args()

    # Load and format data
    spots = pd.read_csv(args.spots_path)
    spots.columns = pd.MultiIndex.from_tuples([("spot", "id"), ("position", "x"), ("position", "y"), *map(parse_image_filename, spots.columns[3:])])
    spots = spots.set_index(("spot", "id")).astype(np.uint32)

    # Run multiple times to ensure only one detected spot remains for each cluster
    for _ in range(3):
        spots = deduplicate_spots(spots, args.r)

    transformation = calculate_transformation(spots)
    print(transformation.coef_)
    print(transformation.intercept_)
