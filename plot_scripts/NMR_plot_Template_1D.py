from pathlib import Path
import matplotlib.pyplot as plt
from nmr_setup import setup_nmr_project, build_plot_dict

# Template_1D for Plot_func V1.4 (1D spectra: 1H, 13C, 19F, ...)
# Parameters that are not listed here get their default from nmr_setup.py

# ==========================================
# 1. AUTO-ROUTING & SYSTEM SETUP
# ==========================================
PROJECT_NAME, paths = setup_nmr_project(__file__)
import Plot_func as pf

# %% =======================================
# 2. IMPORT ALL SPECTRA
# ==========================================

# directory of your NMR files
base_dir = Path("/User/Documents/NMR_Spectra")

# Data Format (folder, name, scaling factor (1 = unchanged), color, label.csv(optional))
# The label CSV needs the columns "Residue" and "position_1"
data = [
    ("250905_Protein_1/1",      "Protein 1",        1,      "black",          "Protein_1_1H.csv"),
    ("250905_Protein_2/1",      "Protein 2",        1,      "tab:blue"),
    ("250905_Protein_3/1",      "Protein 3",        2,      "tab:orange"),
    ("250905_Protein_4/1",      "Protein 4",        1,      [0.839, 0.153, 0.157]),
]

files, file_names, cont, colors, csv_files = pf.parse_plot_data(data, base_dir)
# colors = pf.colormap(files, maps='viridis')  # use a colormap

# %% =======================================
# 3. SETUP PARAMETERS
# ==========================================

dpi = 300               # dpi of the figure
save_png = True         # saves png files
save_svg = False        # saves svg files (slower)

alpha = 0.7             # transparency in overlay_1d
legend = True           # visibility of the legend in overlay plots

# NORMALISATION AND STACKING
normalize = None        # None, "max" (highest peak in xlim), "region" (norm_region) or "ns" (NS*RG)
norm_region = None      # e.g. (0.1, -0.1) for normalize = "region" (reference signal)
stack_offset = 1.0      # vertical distance in stack_1d (1 = height of the highest peak)

# GRIDPLOT PARAMETERS
title_y = 0.88          # height of the panel title in grid_plot_1d (fraction of the panel)
first_x = True          # hides the x-tick label at the high-ppm end in grid panels

# PEAKLABEL PARAMETERS
peak_labels = False     # toggle visibility of peak labels from CSV files
label_fontsize = 4      # size of the text
connector_width = 0.5   # linewidth of the connector
label_rotation_1d = 90  # rotation of the labels
label_offset_1d = 0.04  # distance label - peak (fraction of the intensity range)
label_window_1d = 0.02  # ppm window in which the peak maximum is searched


# %% SPECTRA TYPE PARAMETERS
# -- 1D styles only need the x axis; "ylim_1d" (intensity range) is optional --

SPECTRA_STYLES = {
    "1H": {
        "xlim": (10.5, -0.5),
        "x_label": "$^1$H [ppm]",
        "xticks": 1,
        "xminorticks": 0.2,
        "line_width": 0.6,
        "xsize": 6,
        "ysize": 3,
    },

    "1H_amide": {
        "xlim": (10.5, 6.0),
        "x_label": "$^1$H [ppm]",
        "xticks": 0.5,
        "xminorticks": 0.1,
        "line_width": 0.6,
        "xsize": 6,
        "ysize": 3,
    },

    "1H_methyl": {
        "xlim": (1.5, -0.5),
        "x_label": "$^1$H [ppm]",
        "xticks": 0.5,
        "xminorticks": 0.1,
        "line_width": 0.6,
        "xsize": 6,
        "ysize": 3,
    },

    "13C": {
        "xlim": (200, 0),
        "x_label": "$^{13}$C [ppm]",
        "xticks": 20,
        "xminorticks": 5,
        "line_width": 0.5,
        "xsize": 6,
        "ysize": 3,
    },
}

# %% ADDITIONAL SETUP PARAMETERS

# LINEWIDTH AND FONT SIZE
plt.rcParams["axes.linewidth"] = 1
plt.rcParams["svg.fonttype"] = "none"  # set to "path" if you want the text as curves

tick_width = 1
plt.rcParams["xtick.major.width"] = tick_width
plt.rcParams["xtick.minor.width"] = tick_width

axis_labels = True
titles = True
plt.rcParams["font.size"] = 10
plt.rcParams["axes.titlesize"] = 12
plt.rcParams["axes.labelsize"] = 10
plt.rcParams["xtick.labelsize"] = 10
plt.rcParams["legend.fontsize"] = 8
plt.rcParams["figure.titlesize"] = 20

# FONT FAMILY
# -- Choose ARIAL or TIMES--
Arial = True
if Arial:
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans"]
    plt.rcParams["mathtext.fontset"] = "dejavusans"
else:
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = ["Times New Roman", "DejaVu Serif"]
    plt.rcParams["mathtext.fontset"] = "stix"

plt.rcParams["xtick.direction"] = "out"
plt.rcParams["axes.facecolor"] = "None"

# %% LOAD DATA INTO THE MEMORY
print("Loading NMR data into memory...")
dic_all, data_all = pf.read_data(files)

# %% LOAD AND TRANSFER INFO PARAMETERS TO THE FUNCTIONS
p = build_plot_dict(**locals())

# %% SPECTRUM MENU
# this gives you the number of the spectrum for the overlay plot
pf.spectrum_menu(p)

# %% =======================================
# 4. PLOT YOUR 1D SPECTRA
# ==========================================

# -- set your spectra style --
pf.set_style(p, "1H")

# -- plot all 1D spectra --
if True:
    pf.plot_everything_1d(p)

# %% OVERLAYS
# -- you need to define a "name" for the overlay plot, the last number will be on top --
if False:
    pf.overlay_1d(p, [0, 1], name="over1")
    pf.overlay_1d(p, [1, 0], name="over2")

    # -- you can change the style for different plots --
    pf.set_style(p, "1H_amide")
    pf.overlay_1d(p, [0, 1], name="over1")
    pf.set_style(p, "1H")

# %% STACKED PLOTS
# -- the first number is at the bottom, the last one at the top --
# -- names: "text" (next to the spectra), "legend" or None --
if False:
    pf.stack_1d(p, [0, 1, 2, 3], name="stack")
    pf.stack_1d(p, [0, 1, 2, 3], name="stack_close", offset=0.5)          # closer together
    pf.stack_1d(p, [0, 1, 2, 3], name="stack_3D", shift=0.1, names="legend")  # pseudo-3D look

# %% GRID
# -- indices=None plots all 1D spectra; sharey=True uses one intensity scale --
if False:
    pf.grid_plot_1d(p, row=4, col=1)
    pf.grid_plot_1d(p, indices=[0, 1], row=1, col=2, name="grid_pair", sharey=True)

# %% NORMALISATION
# -- compare intensities: by the highest peak, by a reference signal or by NS*RG --
if False:
    p["normalize"] = "max"
    pf.overlay_1d(p, [0, 1], name="over_norm_max")

    p["normalize"] = "region"
    p["norm_region"] = (0.1, -0.1)
    pf.overlay_1d(p, [0, 1], name="over_norm_ref")

    p["normalize"] = "ns"
    pf.overlay_1d(p, [0, 1], name="over_norm_ns")

    p["normalize"] = None
