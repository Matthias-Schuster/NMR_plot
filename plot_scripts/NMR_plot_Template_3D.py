from pathlib import Path
import matplotlib.pyplot as plt
from nmr_setup import setup_nmr_project, build_plot_dict

# Template_3D for Plot_func V1.4 (2D strips from 3D spectra: HNCO, HNCACB, ...)
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

# Data Format (folder, name, contour, color (second color for negative peaks))
data = [
    ("250905_Protein_4/701",    "HNCO",             6e8,    ["tab:blue", "tab:orange"]),
    ("250905_Protein_4/702",    "HNcaCO",           6e8,    ["tab:green", "tab:red"]),
    ("250905_Protein_4/703",    "HNCACB",           4e8,    ["black", "tab:red"]),
    ("250905_Protein_4/704",    "CBCAcoNH",         4e8,    ["tab:blue", "tab:orange"]),
]

files, file_names, cont, colors, csv_files = pf.parse_plot_data(data, base_dir)
# colors = pf.colormap(files, maps='viridis')  # use a colormap

# %% =======================================
# 3. SETUP PARAMETERS
# ==========================================

dpi = 300               # dpi of the figure
save_png = True         # saves png files
save_svg = False        # saves svg files (makes it much slower)

lines = 25              # number of contour lines
factor = 1.2            # distance between contour lines

labelpad_x = -8         # distance of the x label to the axis in individual strips
labelpad_y = 9          # distance of the y label to the axis in individual strips

alpha = 0.7             # transparency for overlay
first_x = False         # hides the x-tick label at the high-ppm end (left)
first_y = True          # hides the y-tick label at the low-ppm end (top)
legend = True           # visibility of the legend in overlay strips

# STRIP GRID PARAMETERS
title_y = 0.83          # distance of the title to the axis in the strip grid
grid_x = 0.03           # distance of the x label to the axis in the strip grid
grid_y = 0.045          # distance of the y label to the axis in the strip grid


# %% SPECTRA TYPE PARAMETERS
# -- xlim is replaced by x_slice +/- strip_width/2 when an x_slice is given --

SPECTRA_STYLES = {
    "HNCO": {
        "xlim": (8.0, 7.7),
        "ylim": (135.9, 125.8),
        "x_label": "$^1$H [ppm]",
        "y_label": "$^{13}$C\n[ppm]",
        "y_label_grid": "$^{13}$C [ppm]",
        "xticks": 0.1,
        "xminorticks": 0.025,
        "yticks": 2,
        "yminorticks": 1,
        "line_width": 1.5,
        "xsize": 2.5,
        "ysize": 5.0,
    },

    "HNCACB": {
        "xlim": (8.0, 7.7),
        "ylim": (72, 12),
        "x_label": "$^1$H [ppm]",
        "y_label": "$^{13}$C\n[ppm]",
        "y_label_grid": "$^{13}$C [ppm]",
        "xticks": 0.1,
        "xminorticks": 0.025,
        "yticks": 10,
        "yminorticks": 2,
        "line_width": 1.0,
        "xsize": 2.5,
        "ysize": 6.0,
    },
}

# %% ADDITIONAL SETUP PARAMETERS

# LINEWIDTH AND FONT SIZE
plt.rcParams["axes.linewidth"] = 1
plt.rcParams["svg.fonttype"] = "none"  # set to "path" if you want the text as curves

tick_width = 1
plt.rcParams["xtick.major.width"] = tick_width
plt.rcParams["ytick.major.width"] = tick_width
plt.rcParams["xtick.minor.width"] = tick_width
plt.rcParams["ytick.minor.width"] = tick_width

axis_labels = True
titles = True
plt.rcParams["font.size"] = 10
plt.rcParams["axes.titlesize"] = 12
plt.rcParams["axes.labelsize"] = 9
plt.rcParams["xtick.labelsize"] = 10
plt.rcParams["ytick.labelsize"] = 10
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
plt.rcParams["ytick.direction"] = "out"
plt.rcParams["axes.facecolor"] = "None"

# %% LOAD DATA INTO THE MEMORY
print("Loading NMR data into memory...")
dic_all, data_all = pf.read_data(files)

# %% LOAD AND TRANSFER INFO PARAMETERS TO THE FUNCTIONS
p = build_plot_dict(**locals())

# %% SPECTRUM MENU
# this gives you the number of the spectrum for the strip plots
pf.spectrum_menu(p)

# %% =======================================
# 4. PLOT 2D STRIPS FROM YOUR 3D SPECTRA
# ==========================================
# slice_axis: which dimension is fixed by z_slice
#   0: F1 (proc3s) -> F2-F3 plane
#   1: F2 (proc2s) -> F1-F3 plane (default, e.g. 15N plane of an HNCO / HNCACB)
#   2: F3 (procs)  -> F1-F2 plane

# -- set your spectra style --
pf.set_style(p, "HNCO")

# %% SINGLE STRIPS AND STRIP OVERLAYS
# -- z_slice: ppm of the fixed dimension (or slice_index for a point number) --
# -- x_slice: centre of the strip, strip_width: width in ppm, the last number will be on top --
if True:
    pf.plot_strip(p, 0, z_slice=120.75, x_slice=8.23, name="HNCO", strip_width=0.2)
    pf.plot_strip(p, [1, 0], z_slice=122.2, x_slice=7.85, name="HNCO and HNcaCO", strip_width=0.2)

# %% GRID OF STRIPS (e.g. sequential walk)
# -- each slice: a ppm value or a dict with z_slice / slice_index, x_slice and title --
if False:
    my_slices = [
        {"z_slice": 120.75, "x_slice": 8.23, "title": "Asn12"},
        {"z_slice": 121.38, "x_slice": 8.08, "title": "Glu13"},
        {"z_slice": 122.47, "x_slice": 7.85, "title": "Val14"},
    ]

    pf.strip_grid_plot(p, 0, slices=my_slices, row=1, col=3, strip_width=0.2, name="HNCO")
    pf.strip_grid_plot(
        p, [0, 1], slices=my_slices, row=1, col=3, strip_width=0.2, name="HNCO_HNcaCO_walk"
    )

# %% SEQUENTIAL WALK WITH HNCACB / CBCAcoNH
if False:
    pf.set_style(p, "HNCACB")

    walk = [
        {"z_slice": 120.75, "x_slice": 8.23, "title": "Asn12"},
        {"z_slice": 121.38, "x_slice": 8.08, "title": "Glu13"},
        {"z_slice": 122.47, "x_slice": 7.85, "title": "Val14"},
        {"z_slice": 119.90, "x_slice": 8.01, "title": "Lys15"},
    ]

    pf.strip_grid_plot(p, [2, 3], slices=walk, row=1, col=4, strip_width=0.2, name="HNCACB_walk")
