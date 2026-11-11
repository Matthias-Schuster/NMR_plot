from pathlib import Path
import sys


# Output sub-folders inside results/<project>/
OUTPUT_FOLDERS = {
    "out_single": "single",
    "out_overlay": "overlay",
    "out_grid": "grid",
    "out_strip": "strips",
    "out_1d": "1D",
}

# Default values of all plot parameters.
# A variable with the same name in the plot script overrides the default.
DEFAULTS = {
    # Output
    "dpi": 300,
    "save_png": True,
    "save_svg": False,
    # Figure size (inch, or size of one panel in grid plots)
    "xsize": 5,
    "ysize": 4,
    # Contours
    "lines": 25,
    "factor": 1.2,
    "alpha": 0.7,
    # Axes, titles, legend
    "labelpad_x": -8,
    "labelpad_y": 9,
    "titles": True,
    "first_x": True,
    "first_y": True,
    "legend": True,
    # Grid plots
    "title_y": 0.88,
    "grid_x": 0.05,
    "grid_y": 0.06,
    # Peak labels
    "peak_labels": False,
    "label_fontsize": 4,
    "connector_width": 0.5,
    "show_ellipse": False,
    "peak_ellipse_x": 0.06,
    "peak_ellipse_y": 0.6,
    "starting_angle": 42,
    "random_seed": 42,
    "jitter": 0.3,
    "cycles": 10,
    "iterations": 10,
    # 1D spectra
    "normalize": None,
    "norm_region": None,
    "stack_offset": 1.0,
    "label_offset_1d": 0.04,
    "label_rotation_1d": 90,
    "label_window_1d": 0.02,
    "label_headroom_1d": 0.3,
    # 2D with projections
    "proj_size": 0.2,
    "proj_line_width": 0.5,
}

# Parameters that normally come from SPECTRA_STYLES via set_style()
STYLE_KEYS = (
    "xlim",
    "ylim",
    "x_label",
    "y_label",
    "y_label_grid",
    "xticks",
    "xminorticks",
    "yticks",
    "yminorticks",
    "line_width",
    "ylim_1d",
)

# Loaded data and the columns of the data list
DATA_KEYS = ("dic_all", "data_all", "files", "file_names", "cont", "colors")


def setup_nmr_project(
    script_file,
    prefix="NMR_plot_",
    results_folder="results",
    csv_folder="input_csv",
    create_output_dirs=True,
    verbose=True,
):
    """
    Automatic setup for NMR plotting scripts.

    It determines:
    - where the script is located
    - where the project root is
    - the project name from the script filename
    - where results should be saved
    - where project-specific CSV files are stored
    """

    # Location of the currently running plot script
    script_path = Path(script_file).resolve()
    script_dir = script_path.parent

    # Project root is one level above plot_scripts/
    root_dir = script_dir.parent

    # Make root directory modules importable
    root_dir_str = str(root_dir)
    if root_dir_str not in sys.path:
        sys.path.insert(0, root_dir_str)

    # Extract project name from filename
    script_name = script_path.stem
    project_name = script_name[len(prefix):] if script_name.startswith(prefix) else script_name

    # Define paths
    project_out = root_dir / results_folder / project_name
    paths = {
        "script_dir": script_dir,
        "root_dir": root_dir,
        "project_out": project_out,
        "csv_dir": root_dir / csv_folder / project_name,
    }
    for key, folder in OUTPUT_FOLDERS.items():
        paths[key] = project_out / folder

    # Create output folders
    if create_output_dirs:
        for key in (*OUTPUT_FOLDERS, "csv_dir"):
            paths[key].mkdir(parents=True, exist_ok=True)

    if verbose:
        print("\n----- NMR Project Setup -----")
        print(f"Project name: {project_name}")
        print(f"Script dir:   {paths['script_dir']}")
        print(f"Root dir:     {paths['root_dir']}")
        print(f"CSV dir:      {paths['csv_dir']}")
        print(f"Output dir:   {paths['project_out']}")
        print("-----------------------------\n")

    return project_name, paths


def build_plot_dict(**kwargs):
    """
    Builds the 'p' dictionary from the local variables of the main script
    (called as build_plot_dict(**locals())). Missing variables get their DEFAULTS.
    """
    paths = kwargs.get("paths") or {}

    p = {key: kwargs.get(key, default) for key, default in DEFAULTS.items()}
    p.update({key: kwargs.get(key) for key in STYLE_KEYS})
    p.update({key: kwargs.get(key) for key in DATA_KEYS})

    p["labels"] = kwargs.get("axis_labels", True)
    p["csv_files"] = kwargs.get("csv_files", [])
    p["spectra_type"] = kwargs.get("spectra_type")  # older scripts set it directly
    p["SPECTRA_STYLES"] = kwargs.get("SPECTRA_STYLES")

    p["csv_dir"] = paths.get("csv_dir", "")
    for key in OUTPUT_FOLDERS:
        p[key] = paths.get(key, "")

    return p
