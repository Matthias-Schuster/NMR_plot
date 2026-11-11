from pathlib import Path
import matplotlib.pyplot as plt
from nmr_setup import setup_nmr_project, build_plot_dict

# Template_all for Plot_func V1.4 (1D, 2D and 3D spectra in one script)
# Specialised templates: NMR_plot_Template_1D.py, _2D.py, _3D.py

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

# Data Format (folder, name, contour, color (second color for negative peaks), label.csv(optional))
# 1D spectra can be added to the same list: for them "contour" is a scaling factor (1 = unchanged)
data = [
    ("250905_Protein_1/700",    "Protein 1",        7e7,    "black",          "Protein_1_pos.csv"),
    ("250905_Protein_2/700",    "Protein 2",        5e8,    "#000000"),
    ("250905_Protein_3/700",    "Protein 3",        6e8,    ["tab:blue", "tab:orange"]),
    ("250905_Protein_4/700",    "Protein 4",       12e8,    [0.839, 0.153, 0.157]),
    ("250905_Protein_4/701",    "Protein 4 HNCO",   6e8,    ["tab:blue", "tab:orange"]),
    ("250905_Protein_4/702",    "Protein 4 HNcaCO", 6e8,    ["tab:green", "tab:red"]),
    ("250905_Protein_1/1",      "Protein 1 1H",     1,      "black"),
    ("250905_Protein_2/1",      "Protein 2 1H",     1,      "tab:blue"),
]

files, file_names, cont, colors, csv_files = pf.parse_plot_data(data, base_dir)
# colors = pf.colormap(files, maps='viridis')  # use a colormap

# %% =======================================
# 3. SETUP PARAMETERS
# ==========================================

dpi = 300               # dpi of the figure
save_png = True         # saves png files
save_svg = False        # saves svg files (makes it much slower)

xsize = 5               # size of the figure in inch (or of one subplot in the grid-plot)
ysize = 4

lines = 25              # number of contour lines
factor = 1.2            # distance between contour lines

labelpad_x = -8         # distance of the x label to the axis in individual plots
labelpad_y = 9          # distance of the y label to the axis in individual plots

alpha = 0.7             # transparency for overlay
first_x = True          # hides the x-tick label at the high-ppm end (left)
first_y = True          # hides the y-tick label at the low-ppm end (top)
legend = True           # visibility of the legend in overlay plots

# GRIDPLOT PARAMETERS
title_y = 0.88          # distance of the title to the axis in the grid_plot
grid_x = 0.05           # distance of the x label to the axis in the grid_plot
grid_y = 0.06           # distance of the y label to the axis in the grid_plot

# PEAKLABEL PARAMETERS
peak_labels = False     # toggle visibility of peak labels from CSV files
label_fontsize = 4      # size of the text
connector_width = 0.5   # linewidth of the connector

show_ellipse = False    # show ellipses around peaks that will be avoided by the text
peak_ellipse_x = 0.06   # ellipse width in x ppm
peak_ellipse_y = 0.6    # ellipse height in y ppm
starting_angle = 42     # starting angle in degrees
random_seed = 42        # random seed for the initial placement
jitter = 0.3            # strength of the random jitter of the initial placement
cycles = 10             # number of peak and text avoiding cycles
iterations = 10         # iterations for each cycle

# 1D PARAMETERS
normalize = None        # None, "max" (highest peak in xlim), "region" (norm_region) or "ns" (NS*RG)
norm_region = None      # e.g. (0.1, -0.1) for normalize = "region" (reference signal)
stack_offset = 1.0      # vertical distance in stack_1d (1 = height of the highest peak)
label_rotation_1d = 90  # rotation of the 1D peak labels

# PROJECTION PARAMETERS
proj_size = 0.2         # size of the projection panels relative to the 2D plot
proj_line_width = 0.5   # line width of the projections


# %% SPECTRA TYPE PARAMETERS

SPECTRA_STYLES = {
    "15N": {
        "xlim": (11.2, 6.2),
        "ylim": (134, 104),
        "x_label": "$^1$H [ppm]",
        "y_label": "$^{15}$N\n[ppm]",
        "y_label_grid": "$^{15}$N [ppm]",
        "xticks": 1,
        "xminorticks": 0.5,
        "yticks": 5,
        "yminorticks": 1,
        "line_width": 0.3,
    },

    "13C": {
        "xlim": (1.51, -0.58),
        "ylim": (28, 18.2),
        "x_label": "$^1$H [ppm]",
        "y_label": "$^{13}$C\n[ppm]",
        "y_label_grid": "$^{13}$C [ppm]",
        "xticks": 0.5,
        "xminorticks": 0.1,
        "yticks": 5,
        "yminorticks": 1,
        "line_width": 0.6,
    },

    "CON": {
        "xlim": (178, 168),
        "ylim": (140, 108),
        "x_label": "$^{13}$C [ppm]",
        "y_label": "$^{15}$N\n[ppm]",
        "y_label_grid": "$^{15}$N [ppm]",
        "xticks": 2,
        "xminorticks": 0.5,
        "yticks": 5,
        "yminorticks": 1,
        "line_width": 0.3,
    },

    "IDP": {
        "xlim": (8.9, 6.6),
        "ylim": (130, 110),
        "x_label": "$^1$H [ppm]",
        "y_label": "$^{15}$N\n[ppm]",
        "y_label_grid": "$^{15}$N [ppm]",
        "xticks": 0.5,
        "xminorticks": 0.25,
        "yticks": 5,
        "yminorticks": 1,
        "line_width": 0.3,
        "xsize": 4,
        "ysize": 7,
    },

    "ZOOM": {
        "xlim": (8.8, 8),
        "ylim": (118, 114),
        "x_label": "$^1$H [ppm]",
        "y_label": "$^{15}$N\n[ppm]",
        "y_label_grid": "$^{15}$N [ppm]",
        "xticks": 0.2,
        "xminorticks": 0.05,
        "yticks": 2,
        "yminorticks": 0.5,
        "line_width": 0.6,
    },

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
        "grid_x": 0.03,
        "grid_y": 0.045,
        "title_y": 0.83,
        "first_x": False,
    },

    # --- 1D styles (only the x axis is needed; "ylim_1d" is optional) ---
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
# this gives you the number of the spectrum for the overlay plot
pf.spectrum_menu(p)

# %% =======================================
# 4. PLOT YOUR 2D SPECTRA
# ==========================================

# -- set your spectra style --
pf.set_style(p, "15N")

# -- plot all 2D spectra --
if True:
    pf.plot_everything(p)

# %% PLOT INDIVIDUAL SPECTRA OR OVERLAYS
# -- you need to define a "name" for the overlay plot --
# -- you can plot all combinations, the last number will be on top --
if False:
    pf.overlay(p, [1], name="prot1")
    pf.overlay(p, [1, 2], name="over1")
    pf.overlay(p, [2, 1], name="over2")

    # -- you can change the style for different plots --
    pf.set_style(p, "ZOOM")
    pf.overlay(p, [1], name="prot1")
    pf.overlay(p, [1, 2], name="over1")
    pf.overlay(p, [2, 1], name="over2")
    pf.set_style(p, "15N")

# %% PLOT ALL SPECTRA IN A GRID
# -- in "grid_plot" you need to define a row and a col for the grid --
if False:
    pf.grid_plot(p, row=2, col=2)

# %% GRID: ONE SPECTRUM OVERLAID ON ALL OTHERS
# -- in "grid_plot_over" you need to define a spectrum to overlay all others --
if False:
    pf.grid_plot_over(p, over=0, row=2, col=2, reverse=False)

# %% GRID: CUSTOM OVERLAY GROUPS
# -- in "grid_plot_over_xp" you need to define overlay groups for each subplot --
if False:
    pf.grid_plot_over_xp(
        p,
        row=2,
        col=2,
        overlay_groups=[
            [1],
            [2],
            [2, 1],
            [1, 2, 3],
        ],
    )

# %% 2D SPECTRUM WITH 1D PROJECTIONS
# -- top / left: "skyline" (max. projection of the 2D), index of a 1D spectrum, or None --
if False:
    pf.plot_with_projections(p, 0)
    pf.plot_with_projections(p, 0, top=6, left="skyline", name="Protein 1 with 1H")
    pf.plot_with_projections(p, [1, 2], top="skyline", left=None, name="over1")

# %% EXTRACT 2D SLICES FROM A 3D SPECTRUM
if False:
    pf.set_style(p, "HNCO")

    pf.plot_strip(p, 4, z_slice=120.75, x_slice=8.23, name="HNCO", strip_width=0.2)
    pf.plot_strip(p, [5, 4], z_slice=122.2, x_slice=7.85, name="HNCO and HNcaCO", strip_width=0.2)

# %% GRID OF 3D STRIPS (e.g. sequential walk)
if False:
    pf.set_style(p, "HNCO")

    my_slices = [
        {"z_slice": 120.75, "x_slice": 8.23, "title": "Asn12"},
        {"z_slice": 121.38, "x_slice": 8.08, "title": "Glu13"},
        {"z_slice": 122.47, "x_slice": 7.85, "title": "Val14"},
    ]

    pf.strip_grid_plot(
        p,
        4,
        slices=my_slices,
        slice_axis=1,
        row=1,
        col=3,
        strip_width=0.2,
        name="HNCO",
    )

    pf.strip_grid_plot(
        p,
        [4, 5],
        slices=my_slices,
        slice_axis=1,
        row=1,
        col=3,
        strip_width=0.2,
        name="Sequential_Walk",
        reverse=False,
    )

# %% =======================================
# 5. PLOT YOUR 1D SPECTRA
# ==========================================
# -- plot all 1D spectra --
if False:
    pf.set_style(p, "1H")
    pf.plot_everything_1d(p)

# %% 1D OVERLAY, STACK AND GRID
# -- the last number will be on top (overlay) / at the top (stack) --
if False:
    pf.set_style(p, "1H_amide")
    pf.overlay_1d(p, [6, 7], name="1H_overlay")
    pf.stack_1d(p, [6, 7], name="1H_stack")                        # names next to the spectra
    pf.stack_1d(p, [6, 7], name="1H_stack_3D", shift=0.1, names="legend")
    pf.grid_plot_1d(p, row=2, col=1)                               # all 1D spectra
