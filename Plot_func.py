"""
Plot_func V1.4

Publication-quality plots of Bruker NMR data (pdata):
    - 2D spectra: single plots, overlays, grids, overlaid grids
    - 2D spectra with 1D projections (skyline or external 1D spectra)
    - 2D strips extracted from 3D spectra
    - 1D spectra: single plots, overlays, stacked plots, grids

Includes contour caching, negative contours (second colour), CSV peak labels
and switchable spectra styles (SPECTRA_STYLES / set_style).
"""

import contextlib
import os
import warnings
from pathlib import Path

import contourpy
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import nmrglue as ng
import numpy as np
import pandas as pd
from adjustText import adjust_text
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse
from matplotlib.transforms import Bbox


# Contour lines (in ppm coordinates) per spectrum. The key always contains a stable
# identifier of the spectrum (index, slice, ...), never only id() of a temporary array.
_CONTOUR_CACHE = {}

# Style parameters every plot type needs (normally set by set_style)
_KEYS_2D = ("xlim", "ylim", "xticks", "xminorticks", "yticks", "yminorticks", "line_width")
_KEYS_1D = ("xlim", "xticks", "xminorticks", "line_width")

_MISSING = object()


# =====================================================================
# --- GENERAL HELPERS ---
# =====================================================================


def clear_cache():
    """Empties the contour cache. Called automatically by read_data()."""
    _CONTOUR_CACHE.clear()


def _check_style(p, keys, what, grid=False, ylabel=True):
    """Raises a clear error if style parameters are missing (e.g. set_style forgotten)."""
    keys = list(keys)
    if p.get("labels", True):
        keys.append("x_label")
        if ylabel:
            keys.append("y_label_grid" if grid else "y_label")

    missing = [k for k in keys if p.get(k) is None]
    if missing:
        style = p.get("spectra_type")
        hint = (
            f"The style '{style}' does not define it."
            if style
            else "Did you forget pf.set_style(p, ...)?"
        )
        raise ValueError(f"{what}: missing parameter(s) {', '.join(missing)}. {hint}")


def _as_list(indices):
    if isinstance(indices, (int, np.integer)):
        return [int(indices)]
    return list(indices)


def _ndim(p, idx):
    data = p["data_all"][idx]
    if p["dic_all"][idx] is None or data is None:
        return None
    return data.ndim


def _valid(p, indices, ndim, verbose=True):
    """Returns only the indices that belong to loaded spectra with 'ndim' dimensions."""
    valid = []
    for idx in _as_list(indices):
        if _ndim(p, idx) == ndim:
            valid.append(idx)
        elif verbose:
            print(f"  -> Skipping [{idx}] {p['file_names'][idx]}: not a {ndim}D spectrum.")
    return valid


def _split_colors(color):
    """Returns (positive, negative) colour. 'negative' is None for a single colour."""
    if (
        isinstance(color, (list, tuple, np.ndarray))
        and len(color) == 2
        and isinstance(color[0], (str, list, tuple, np.ndarray))
    ):
        return color[0], color[1]
    return color, None


def _in_range(values, limits):
    lo, hi = sorted(limits)
    return (values >= lo) & (values <= hi)


def _crop(ppm, y, region, pad=2):
    """Cuts a trace to the plotted region (keeps SVG files small)."""
    inside = np.nonzero(_in_range(ppm, region))[0]
    if inside.size == 0:
        return ppm, y
    start = max(inside[0] - pad, 0)
    stop = min(inside[-1] + pad + 1, len(ppm))
    return ppm[start:stop], y[start:stop]


def _auto_range(traces, region, top=0.05, bottom=0.03):
    """Intensity range of the visible part of several (ppm, y) traces, with margins."""
    lo, hi = np.inf, -np.inf
    for ppm, y in traces:
        mask = _in_range(ppm, region)
        if mask.any():
            lo = min(lo, float(y[mask].min()))
            hi = max(hi, float(y[mask].max()))
    if not np.isfinite(lo):
        return None
    span = (hi - lo) or 1.0
    return lo - bottom * span, hi + top * span


@contextlib.contextmanager
def _quiet():
    """Silences prints and UserWarnings (adjustText is quite chatty)."""
    with open(os.devnull, "w") as fnull, contextlib.redirect_stdout(
        fnull
    ), warnings.catch_warnings():
        warnings.simplefilter("ignore", category=UserWarning)
        yield


def _zip_axes(axes, items):
    """Pairs grid axes with items; surplus axes get None."""
    if len(items) > len(axes):
        print(f"  -> Warning: {len(items)} panels requested, but the grid has only {len(axes)}.")
    for i, ax in enumerate(axes):
        yield ax, (items[i] if i < len(items) else None)


# =====================================================================
# --- DATA INPUT ---
# =====================================================================


def parse_plot_data(data, base_dir, default_sub_dir="pdata/1"):
    """
    Parses the input data list into individual Series for plotting.
    Handles optional columns like csv_file dynamically.
    Automatically appends the default pdata folder if none is specified in the path.
    """
    df = pd.DataFrame(data)

    # Rename the first 4 columns, and grab the 5th if it exists
    if len(df.columns) >= 5:
        df.columns = ["folder", "name", "contour", "color", "csv_file"] + list(df.columns[5:])
    else:
        df.columns = ["folder", "name", "contour", "color"]
        df["csv_file"] = None  # Create an empty column if no CSVs are provided

    # If the user already wrote "pdata" in the folder name, use it as-is.
    # Otherwise, append the default_sub_dir (pdata/1).
    df["full_path"] = df["folder"].apply(
        lambda x: str(base_dir / x) if "pdata" in str(x) else str(base_dir / x / default_sub_dir)
    )

    return (df["full_path"], df["name"], df["contour"], df["color"], df["csv_file"])


def read_data(files):
    """Reads all Bruker pdata folders. Missing or unreadable folders become None."""
    clear_cache()  # new data -> old contours are no longer valid
    data_all, dic_all = [], []
    for f in files:
        dic, data = None, None
        if not Path(f).exists():
            print(f"  -> Path missing: {f}")
        else:
            try:
                dic, data = ng.bruker.read_pdata(str(f))
            except Exception as e:
                print(f"  -> Could not read {f}: {e}")
        dic_all.append(dic)
        data_all.append(data)
    return dic_all, data_all


def _ppm_axis(size, offset_ppm, sw_hz, sf_mhz):
    sw_ppm = float(sw_hz) / float(sf_mhz)
    return float(offset_ppm) - (np.arange(size, dtype=float) * (sw_ppm / size))


def _proc_ppm(proc, size):
    return _ppm_axis(size, proc["OFFSET"], proc["SW_p"], proc["SF"])


def read_hsqc_bruker_pdata(dic, data):
    """ppm axes (x = direct dimension, y = indirect dimension) of a 2D spectrum."""
    x_ppm = _proc_ppm(dic["procs"], data.shape[1])
    y_ppm = _proc_ppm(dic["proc2s"], data.shape[0])
    return x_ppm, y_ppm


# =====================================================================
# --- SAVING AND FORMATTING ---
# =====================================================================


def result_folder(name):
    Path(name).mkdir(parents=True, exist_ok=True)


def save_and_clear(fig, folder, name, p):
    """Handles multi-format saving, showing, and memory cleanup."""

    if p.get("spectra_type"):
        name = f"{name}_{p['spectra_type']}"

    # Append the extension as text: with_suffix() would treat ".75ppm" in
    # "strip_120.75ppm" as an extension and replace it (-> files overwrite each other)
    save_base = Path(folder) / name

    if p.get("save_png", False):
        fig.savefig(
            save_base.parent / f"{save_base.name}.png",
            transparent=False,
            bbox_inches="tight",
            dpi=p.get("dpi", 300),
        )
    if p.get("save_svg", False):
        fig.savefig(save_base.parent / f"{save_base.name}.svg", transparent=True, bbox_inches="tight")

    if plt.get_backend().lower() != "agg":
        plt.show()

    plt.close(fig)


class _HideEdgeLabel(mticker.Formatter):
    """
    Wraps a tick formatter and blanks the label of the visible tick with the
    highest ('max') or lowest ('min') value. Works on the tick values, so it stays
    correct after redraws and limit changes (unlike hiding a label by list index).
    """

    def __init__(self, base, edge):
        self.base = base
        self.edge = edge

    def set_axis(self, axis):
        super().set_axis(axis)
        self.base.set_axis(axis)

    def __call__(self, x, pos=None):
        return self.base(x, pos)

    def get_offset(self):
        return self.base.get_offset()

    def format_data_short(self, value):
        return self.base.format_data_short(value)

    def format_ticks(self, values):
        labels = list(self.base.format_ticks(values))
        if self.axis is None or len(values) == 0:
            return labels
        lo, hi = sorted(self.axis.get_view_interval())
        tol = 1e-9 * max(abs(hi - lo), 1.0)
        visible = [i for i, v in enumerate(values) if lo - tol <= v <= hi + tol]
        if visible:
            pick = max if self.edge == "max" else min
            labels[pick(visible, key=lambda i: values[i])] = ""
        return labels


def _hide_edge_label(axis, edge):
    fmt = axis.get_major_formatter()
    if isinstance(fmt, _HideEdgeLabel):  # shared axes share their formatter
        fmt.edge = edge
    else:
        axis.set_major_formatter(_HideEdgeLabel(fmt, edge))


def _set_title(ax, p, title, is_grid=False, loc="center"):
    if p.get("titles", True) and title:
        ax.set_title(title, y=p["title_y"] if is_grid else 1.0, loc=loc)


def _format_x(ax, p, hide_first=False):
    ax.set_xlim(p["xlim"])
    ax.xaxis.set_major_locator(mticker.MultipleLocator(p["xticks"]))
    ax.xaxis.set_minor_locator(mticker.MultipleLocator(p["xminorticks"]))
    # Hides the label at the high-ppm end (left), where it clashes with the
    # x-axis label or with the neighbouring panel in grid plots
    if hide_first:
        _hide_edge_label(ax.xaxis, "max")


def apply_formatting(ax, p, title=None, is_grid=False, add_labels=False):
    """Handles titles, limits, ticks, and optionally adds axis labels (2D plots)."""
    _set_title(ax, p, title, is_grid)
    _format_x(ax, p, hide_first=p.get("first_x", False))

    ax.set_ylim(p["ylim"])
    ax.yaxis.set_major_locator(mticker.MultipleLocator(p["yticks"]))
    ax.yaxis.set_minor_locator(mticker.MultipleLocator(p["yminorticks"]))

    if not is_grid:
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")

    # Hides the label at the low-ppm end of y (top)
    if p.get("first_y", False):
        _hide_edge_label(ax.yaxis, "min")

    if add_labels and p["labels"]:
        ax.set_xlabel(p["x_label"], loc="left", labelpad=p["labelpad_x"], backgroundcolor="None")
        ax.set_ylabel(
            p["y_label"],
            loc="top",
            labelpad=p["labelpad_y"],
            rotation="horizontal",
            ma="center",
            backgroundcolor="None",
        )


def _new_figure(p):
    fig = plt.figure(dpi=p["dpi"], figsize=(p["xsize"], p["ysize"]))
    return fig, fig.add_subplot()


def _new_grid(p, row, col, sharex=True, sharey=True, ylabel=True, suplabels=True):
    fig, axes = plt.subplots(
        row,
        col,
        dpi=p["dpi"],
        figsize=(p["xsize"] * col, p["ysize"] * row),
        sharex=sharex,
        sharey=sharey,
        gridspec_kw={"wspace": 0, "hspace": 0},
    )
    if p["labels"] and suplabels:
        fig.supxlabel(p["x_label"], y=p["grid_x"])
        if ylabel:
            fig.supylabel(p["y_label_grid"], x=p["grid_y"])
    return fig, np.atleast_1d(axes).flatten()


def _add_legend(ax, handles, names, p):
    if handles and p.get("legend", True):
        leg = ax.legend(handles, names, loc="upper left", framealpha=0.8, handlelength=1.5)
        for line in leg.get_lines():
            line.set_linewidth(2.0)


# =====================================================================
# --- CONTOURS ---
# =====================================================================


def _contour_lines(generator, levels):
    lines = []
    for level in levels:
        lines.extend(seg for seg in generator.lines(level) if len(seg) > 1)
    return lines


def draw_contours(ax, dic, data, p, contour_start, color, cache_id=None):
    """
    Draws positive (solid) and, if two colours are given, negative (dashed) contours.

    The contour lines are cached under 'cache_id' (any hashable that identifies the
    spectrum, e.g. ("2d", idx, id(data))), so repeated plots are fast.
    cache_id=None disables caching. Returns a legend handle.
    """
    pos_color, neg_color = _split_colors(color)
    negative = neg_color is not None
    levels = contour_start * p["factor"] ** np.arange(p["lines"])

    key = None
    if cache_id is not None:
        key = (cache_id, data.shape, float(contour_start), p["factor"], p["lines"], negative)

    cached = _CONTOUR_CACHE.get(key) if key is not None else None
    if cached is None:
        x_ppm, y_ppm = read_hsqc_bruker_pdata(dic, data)
        generator = contourpy.contour_generator(
            x_ppm, y_ppm, np.asarray(data, dtype=float), line_type="Separate"
        )
        pos_paths = _contour_lines(generator, levels)
        neg_paths = _contour_lines(generator, np.sort(-levels)) if negative else []
        cached = (pos_paths, neg_paths)
        if key is not None:
            _CONTOUR_CACHE[key] = cached

    pos_paths, neg_paths = cached
    style = {"linewidths": p["line_width"], "alpha": p["alpha"]}

    ax.add_collection(LineCollection(pos_paths, colors=[pos_color], **style), autolim=False)
    if negative and neg_paths:
        ax.add_collection(
            LineCollection(neg_paths, colors=[neg_color], linestyles="dashed", **style),
            autolim=False,
        )

    return Line2D([0], [0], color=pos_color, lw=p["line_width"])


def _draw_2d(ax, p, indices):
    """Draws the contours of several 2D spectra (in this order) into one axes."""
    handles, names = [], []
    for idx in indices:
        data = p["data_all"][idx]
        handles.append(
            draw_contours(
                ax,
                p["dic_all"][idx],
                data,
                p,
                p["cont"][idx],
                p["colors"][idx],
                cache_id=("2d", idx, id(data)),
            )
        )
        names.append(p["file_names"][idx])
    return handles, names


# =====================================================================
# --- LABEL FUNCTIONS ---
# =====================================================================


def _has_csv(csv_file):
    return isinstance(csv_file, (str, Path)) and str(csv_file) != ""


def _read_label_csv(csv_file, p, columns):
    """Reads a label CSV from csv_dir. Returns None if disabled, missing or invalid."""
    if not p.get("peak_labels", True) or not _has_csv(csv_file):
        return None

    csv_path = Path(p.get("csv_dir") or "") / csv_file
    if not csv_path.exists():
        print(f"  -> Warning: CSV file not found at {csv_path.resolve()}")
        return None

    try:
        df = pd.read_csv(csv_path)
    except Exception as e:
        print(f"  -> Error reading {csv_file}: {e}")
        return None

    missing = [c for c in columns if c not in df.columns]
    if missing:
        print(f"  -> Warning: CSV file {csv_file} is missing the column(s) {', '.join(missing)}.")
        return None
    return df


def add_labels_from_csv(ax, csv_file, p):
    """
    Reads positions from a CSV, creates text objects on the axis,
    and returns the objects and coordinates to be adjusted later.
    """
    df = _read_label_csv(csv_file, p, ("Residue", "position_1", "position_2"))
    if df is None:
        return [], [], []

    xmin, xmax = sorted(ax.get_xlim())
    ymin, ymax = sorted(ax.get_ylim())
    df = df[
        (df["position_1"] >= xmin)
        & (df["position_1"] <= xmax)
        & (df["position_2"] >= ymin)
        & (df["position_2"] <= ymax)
    ]

    texts = [
        ax.text(
            row.position_1,
            row.position_2,
            str(row.Residue),
            fontsize=p.get("label_fontsize", 6),
            color="black",
            ha="center",
            va="center",
            clip_on=False,
        )
        for row in df.itertuples()
    ]
    return texts, df["position_1"].tolist(), df["position_2"].tolist()


def _collect_labels(ax, p, indices):
    """Creates the CSV labels of several spectra and pools them for adjust_all_labels."""
    texts, xs, ys = [], [], []
    for idx in indices:
        t, x, y = add_labels_from_csv(ax, p["csv_files"][idx], p)
        texts += t
        xs += x
        ys += y
    return texts, xs, ys


# --- LABEL ADJUSTMENT HELPER FUNCTIONS ---


def hide_texts_outside_axes(ax, texts, arrows=None):
    """Hide text labels and arrows whose bounding box lies outside the axes."""
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    ax_bbox = ax.get_window_extent(renderer)

    for i, text in enumerate(texts):
        text_bbox = text.get_window_extent(renderer)
        inside = (
            text_bbox.x0 >= ax_bbox.x0
            and text_bbox.x1 <= ax_bbox.x1
            and text_bbox.y0 >= ax_bbox.y0
            and text_bbox.y1 <= ax_bbox.y1
        )
        if not inside:
            text.set_visible(False)
            if arrows is not None and i < len(arrows) and arrows[i] is not None:
                arrows[i].set_visible(False)


def _ellipse_bboxes_from_data(ax, x, y, width, height, expand=1.0):
    """Compute display-space bounding boxes for ellipses directly from data coords."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    half_w, half_h = width / 2.0, height / 2.0
    corners_1 = np.column_stack([x - half_w, y - half_h])
    corners_2 = np.column_stack([x + half_w, y + half_h])

    disp_1 = ax.transData.transform(corners_1)
    disp_2 = ax.transData.transform(corners_2)

    x0 = np.minimum(disp_1[:, 0], disp_2[:, 0])
    x1 = np.maximum(disp_1[:, 0], disp_2[:, 0])
    y0 = np.minimum(disp_1[:, 1], disp_2[:, 1])
    y1 = np.maximum(disp_1[:, 1], disp_2[:, 1])

    if expand != 1.0:
        cx, cy = 0.5 * (x0 + x1), 0.5 * (y0 + y1)
        hw, hh = 0.5 * (x1 - x0) * expand, 0.5 * (y1 - y0) * expand
        x0, x1, y0, y1 = cx - hw, cx + hw, cy - hh, cy + hh

    return [Bbox.from_extents(a, b, c, d) for a, b, c, d in zip(x0, y0, x1, y1)]


def _bbox_list_to_arrays(bboxes):
    """Convert list of Bbox objects to numpy arrays for fast vectorized checks."""
    arr = np.array([[b.x0, b.y0, b.x1, b.y1] for b in bboxes], dtype=float)
    return arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3]


def push_texts_out_of_obstacle_bboxes_fast(
    fig, ax, texts, obstacle_bboxes, padding_px=5, max_iter=80
):
    """Vectorized hard constraint step. Pushes labels out of peak ellipses."""
    if not texts or not obstacle_bboxes:
        return

    obs_x0, obs_y0, obs_x1, obs_y1 = _bbox_list_to_arrays(obstacle_bboxes)
    # Text extents only need a renderer, not a full (slow) redraw of all contours
    renderer = fig.canvas.get_renderer()

    for _ in range(max_iter):
        moved_any = False

        for text in texts:
            if not text.get_visible():
                continue

            tb = text.get_window_extent(renderer)
            tx0, tx1 = min(tb.x0, tb.x1), max(tb.x0, tb.x1)
            ty0, ty1 = min(tb.y0, tb.y1), max(tb.y0, tb.y1)

            # Vectorized overlap check against all obstacles simultaneously
            overlaps = (tx0 < obs_x1) & (tx1 > obs_x0) & (ty0 < obs_y1) & (ty1 > obs_y0)
            if not np.any(overlaps):
                continue

            ox0, oy0, ox1, oy1 = (
                obs_x0[overlaps],
                obs_y0[overlaps],
                obs_x1[overlaps],
                obs_y1[overlaps],
            )

            # Pixel moves that would separate the text from each overlapping obstacle side
            candidate_moves = (
                [(dx, 0.0) for dx in ox0 - tx1 - padding_px]  # left
                + [(dx, 0.0) for dx in ox1 - tx0 + padding_px]  # right
                + [(0.0, dy) for dy in oy0 - ty1 - padding_px]  # down
                + [(0.0, dy) for dy in oy1 - ty0 + padding_px]  # up
            )

            # Pick the shortest escape route
            dx_px, dy_px = min(candidate_moves, key=lambda move: abs(move[0]) + abs(move[1]))

            old_display = ax.transData.transform(text.get_position())
            new_display = old_display + np.array([dx_px, dy_px])
            text.set_position(ax.transData.inverted().transform(new_display))
            moved_any = True

        if not moved_any:
            break


def draw_final_connectors(ax, x, y, texts, min_distance_px=2, lw=0.5, alpha=0.5, color="black"):
    """Draw connector lines once, after final label positions are known."""
    arrows = []
    for px, py, text in zip(x, y, texts):
        if not text.get_visible():
            arrows.append(None)
            continue

        tx, ty = text.get_position()

        # Calculate pixel distance to decide if we need a line at all
        anchor_disp = ax.transData.transform((px, py))
        text_disp = ax.transData.transform((tx, ty))
        distance_px = np.hypot(text_disp[0] - anchor_disp[0], text_disp[1] - anchor_disp[1])

        if distance_px <= min_distance_px:
            arrows.append(None)
            continue

        arrow = ax.annotate(
            "",
            xy=(px, py),
            xytext=(tx, ty),
            arrowprops=dict(
                arrowstyle="-",
                color=color,
                lw=lw,
                alpha=alpha,
                patchA=text,
                shrinkA=0.25,
                shrinkB=0,
            ),
            annotation_clip=False,
            zorder=4,
        )
        arrows.append(arrow)

    return arrows


# --- MAIN LABEL ADJUSTMENT FUNCTION ---


def adjust_all_labels(ax, texts, x, y, p):
    """Adjust peak labels while avoiding peak ellipses using iterative cooling."""
    if not texts:
        return [], []

    fig = ax.figure
    fig.canvas.draw()

    n = min(len(texts), len(x), len(y))
    if n == 0:
        return [], []

    texts = texts[:n]
    x = np.asarray(x[:n], dtype=float)
    y = np.asarray(y[:n], dtype=float)

    # User Settings
    ellipse_width = p.get("peak_ellipse_x", 0.06)
    ellipse_height = p.get("peak_ellipse_y", 0.6)
    show_ellipse = p.get("show_ellipse", False)
    random_seed = p.get("random_seed", 42)
    jitter = p.get("jitter", 0.3)

    # Internal parameters for placement
    bbox_expand = 1.05
    initial_offset_scale = 1.05
    expand = (1.25, 1.25)
    obstacle_padding_px = 1

    # Iteration Control ("Cooling" parameters)
    cycles = p.get("cycles", 10)  # Number of push/relax cycles
    base_max_move = 30.0  # Starting movement limit in pixels
    base_pull_force = 0.1  # Starting pull back to the peak
    iterations = p.get("iterations", 10)  # Number of max iterations

    # 1. Build obstacle boxes
    obstacle_bboxes = _ellipse_bboxes_from_data(
        ax=ax, x=x, y=y, width=ellipse_width, height=ellipse_height, expand=bbox_expand
    )

    if show_ellipse:
        for px, py in zip(x, y):
            ax.add_patch(
                Ellipse(
                    (px, py),
                    width=ellipse_width,
                    height=ellipse_height,
                    fill=True,
                    facecolor="black",
                    edgecolor="black",
                    alpha=0.25,
                    zorder=1,
                )
            )

    # 2. Initial deterministic label placement (golden-angle spiral + reproducible jitter)
    golden_angle = np.pi * (3.0 - np.sqrt(5.0))
    starting_phase = np.radians(p.get("starting_angle", 0.0))
    angles = (np.arange(n, dtype=float) * golden_angle) + starting_phase

    base_dx = np.cos(angles) * ellipse_width * initial_offset_scale
    base_dy = np.sin(angles) * ellipse_height * initial_offset_scale

    rng = np.random.default_rng(random_seed)
    dx = base_dx + rng.uniform(-jitter, jitter, n) * ellipse_width
    dy = base_dy + rng.uniform(-jitter, jitter, n) * ellipse_height

    for i, text in enumerate(texts):
        text.set_position((x[i] + dx[i], y[i] + dy[i]))

    # 3. Iterative Solver Loop
    for cycle in range(cycles):

        # A. Hard constraint: Push entirely out of peaks
        push_texts_out_of_obstacle_bboxes_fast(
            fig=fig,
            ax=ax,
            texts=texts,
            obstacle_bboxes=obstacle_bboxes,
            padding_px=obstacle_padding_px,
            max_iter=iterations,
        )

        # B. Soft constraint: Resolve text overlaps
        # Decrease movement allowance and pull force as cycles progress
        current_max_move = base_max_move / (cycle + 1)
        current_pull = base_pull_force / (cycle + 1) if cycle < (cycles - 1) else 0.0

        with _quiet():
            adjust_text(
                texts,
                target_x=x,
                target_y=y,
                avoid_self=False,
                ax=ax,
                prevent_crossings=True,
                ensure_inside_axes=True,
                expand_axes=False,
                expand=expand,
                force_text=(0.6, 0.6),
                force_static=(0.0, 0.0),
                force_pull=(current_pull, current_pull),
                max_move=(current_max_move, current_max_move),
                iter_lim=iterations,
            )

    # 4. One final hard push to guarantee no peak is covered after the final adjust_text
    push_texts_out_of_obstacle_bboxes_fast(
        fig=fig,
        ax=ax,
        texts=texts,
        obstacle_bboxes=obstacle_bboxes,
        padding_px=obstacle_padding_px,
        max_iter=iterations,
    )

    # 5. Draw final connector lines
    arrows = draw_final_connectors(
        ax=ax,
        x=x,
        y=y,
        texts=texts,
        min_distance_px=4,
        lw=p.get("connector_width", 0.5),
        alpha=0.5,
        color="black",
    )

    hide_texts_outside_axes(ax, texts, arrows)

    return texts, arrows


# =====================================================================
# --- 2D PLOTTING FUNCTIONS ---
# =====================================================================


def plot_everything(p):
    """Plots every loaded 2D spectrum individually."""
    _check_style(p, _KEYS_2D, "plot_everything")
    folder = p["out_single"]
    result_folder(folder)

    for idx in _valid(p, range(len(p["files"])), 2, verbose=False):
        fig, ax = _new_figure(p)
        apply_formatting(ax, p, title=p["file_names"][idx], add_labels=True)
        _draw_2d(ax, p, [idx])
        adjust_all_labels(ax, *_collect_labels(ax, p, [idx]), p)
        save_and_clear(fig, folder, f"{p['file_names'][idx]}", p)


def overlay(p, get, name):
    """Overlays the 2D spectra 'get' (list of indices); the last one is on top."""
    _check_style(p, _KEYS_2D, "overlay")
    indices = _valid(p, get, 2)
    if not indices:
        print("  -> Error: No valid 2D spectra found.")
        return

    folder = p["out_overlay"]
    result_folder(folder)

    fig, ax = _new_figure(p)
    apply_formatting(ax, p, title=name, add_labels=True)

    handles, names = _draw_2d(ax, p, indices)
    # Pool the labels of all spectra so they avoid each other
    adjust_all_labels(ax, *_collect_labels(ax, p, indices), p)
    _add_legend(ax, handles, names, p)

    save_and_clear(fig, folder, f"overlay_{name}", p)


def grid_plot(p, row=2, col=2, name="grid"):
    """Plots all 2D spectra in a row x col grid."""
    _check_style(p, _KEYS_2D, "grid_plot", grid=True)
    folder = p["out_grid"]
    result_folder(folder)

    indices = _valid(p, range(len(p["files"])), 2, verbose=False)
    fig, axes = _new_grid(p, row, col)

    for ax, idx in _zip_axes(axes, indices):
        if idx is None:
            ax.axis("off")
            continue
        apply_formatting(ax, p, title=p["file_names"][idx], is_grid=True)
        _draw_2d(ax, p, [idx])
        adjust_all_labels(ax, *_collect_labels(ax, p, [idx]), p)

    save_and_clear(fig, folder, name, p)


def grid_plot_over(p, over, row=2, col=2, reverse=False, name="grid_over"):
    """Grid of all 2D spectra, each overlaid with spectrum 'over' (on top unless reverse)."""
    _check_style(p, _KEYS_2D, "grid_plot_over", grid=True)
    if _ndim(p, over) != 2:
        print(f"  -> Error: Overlay spectrum [{over}] is missing or not a 2D spectrum.")
        return

    folder = p["out_grid"]
    result_folder(folder)

    indices = [i for i in _valid(p, range(len(p["files"])), 2, verbose=False) if i != over]
    fig, axes = _new_grid(p, row, col)

    for ax, idx in _zip_axes(axes, indices):
        if idx is None:
            ax.axis("off")
            continue
        apply_formatting(ax, p, title=p["file_names"][idx], is_grid=True)
        _draw_2d(ax, p, [over, idx] if reverse else [idx, over])
        adjust_all_labels(ax, *_collect_labels(ax, p, [idx, over]), p)

    save_and_clear(fig, folder, name, p)


def grid_plot_over_xp(p, overlay_groups, row=2, col=2, name="grid_over_xp"):
    """Grid in which every panel shows its own overlay group (list of index lists)."""
    _check_style(p, _KEYS_2D, "grid_plot_over_xp", grid=True)
    folder = p["out_grid"]
    result_folder(folder)

    fig, axes = _new_grid(p, row, col)

    for ax, group in _zip_axes(axes, overlay_groups):
        if group is None:
            ax.axis("off")
            continue
        apply_formatting(ax, p, is_grid=True)
        group = _valid(p, group, 2)
        handles, names = _draw_2d(ax, p, group)
        adjust_all_labels(ax, *_collect_labels(ax, p, group), p)
        _add_legend(ax, handles, names, p)

    save_and_clear(fig, folder, name, p)


# --- 2D WITH 1D PROJECTIONS ---


def _skyline(p, idx, axis):
    """
    Maximum ("skyline") projection of the shown region of a 2D spectrum.
    axis="x": onto the direct dimension (top panel), axis="y": onto the indirect one (left).
    """
    dic, data = p["dic_all"][idx], p["data_all"][idx]
    x_ppm, y_ppm = read_hsqc_bruker_pdata(dic, data)
    if axis == "x":
        rows = _in_range(y_ppm, p["ylim"])
        return x_ppm, (data[rows] if rows.any() else data).max(axis=0)
    cols = _in_range(x_ppm, p["xlim"])
    return y_ppm, (data[:, cols] if cols.any() else data).max(axis=1)


def _projection_traces(p, spec, indices_2d, axis):
    """Returns [(ppm, intensity, colour)] for one projection panel."""
    region = p["xlim"] if axis == "x" else p["ylim"]
    traces = []
    if isinstance(spec, str):
        if spec != "skyline":
            raise ValueError(f"Unknown projection '{spec}'. Use 'skyline', a 1D index or None.")
        for idx in indices_2d:
            ppm, y = _skyline(p, idx, axis)
            traces.append((ppm, y, _split_colors(p["colors"][idx])[0]))
    else:
        for idx in _valid(p, spec, 1):
            ppm, y = get_1d(p, idx, region=region)
            traces.append((ppm, y, _split_colors(p["colors"][idx])[0]))
    return traces


def _draw_projection(ax, p, traces, axis):
    region = p["xlim"] if axis == "x" else p["ylim"]
    lw = p.get("proj_line_width", 0.5)
    for ppm, y, color in traces:
        ppm_c, y_c = _crop(ppm, y, region)
        if axis == "x":
            ax.plot(ppm_c, y_c, color=color, lw=lw)
        else:
            ax.plot(y_c, ppm_c, color=color, lw=lw)

    limits = _auto_range([(ppm, y) for ppm, y, _ in traces], region)
    if limits is not None:
        if axis == "x":
            ax.set_ylim(limits)
        else:
            ax.set_xlim(limits[1], limits[0])  # peaks point to the left, away from the 2D
    ax.axis("off")


def plot_with_projections(p, indices, top="skyline", left="skyline", name=None):
    """
    2D spectrum (or overlay of several) with 1D projections above and to the left.

    top / left:
        "skyline"   maximum projection of the shown region of the 2D spectra
        int / list  index (or indices) of 1D spectra from the data list
        None        no projection on this side
    Saved in the 'single' folder as <name>_projections.
    """
    _check_style(p, _KEYS_2D, "plot_with_projections")
    indices = _valid(p, indices, 2)
    if not indices:
        print("  -> Error: No valid 2D spectra found.")
        return

    folder = p["out_single"]
    result_folder(folder)

    size = p.get("proj_size", 0.2)
    widths = ([size] if left is not None else []) + [1.0]
    heights = ([size] if top is not None else []) + [1.0]

    fig = plt.figure(dpi=p["dpi"], figsize=(p["xsize"] * sum(widths), p["ysize"] * sum(heights)))
    gs = fig.add_gridspec(
        len(heights),
        len(widths),
        width_ratios=widths,
        height_ratios=heights,
        wspace=0.02,
        hspace=0.02,
    )
    ax = fig.add_subplot(gs[-1, -1])
    ax_top = fig.add_subplot(gs[0, -1], sharex=ax) if top is not None else None
    ax_left = fig.add_subplot(gs[-1, 0], sharey=ax) if left is not None else None

    file_names = [p["file_names"][i] for i in indices]
    apply_formatting(ax, p, add_labels=True)
    _set_title(ax_top or ax, p, name or " / ".join(file_names))

    handles, names = _draw_2d(ax, p, indices)
    adjust_all_labels(ax, *_collect_labels(ax, p, indices), p)
    if len(indices) > 1:
        _add_legend(ax, handles, names, p)

    if ax_top is not None:
        _draw_projection(ax_top, p, _projection_traces(p, top, indices, "x"), "x")
    if ax_left is not None:
        _draw_projection(ax_left, p, _projection_traces(p, left, indices, "y"), "y")

    save_and_clear(fig, folder, f"{name or '_'.join(file_names)}_projections", p)


# =====================================================================
# --- 3D STRIP FUNCTIONS ---
# =====================================================================

# Bruker processing parameter blocks for data axis 0, 1, 2 of a 3D spectrum
_PROC_KEYS_3D = ("proc3s", "proc2s", "procs")


def get_dim_ppm(dic_3d, data_3d, dim):
    """Calculates the PPM axis for a 3D dimension (0: F1/proc3s, 1: F2/proc2s, 2: F3/procs)."""
    if dim not in (0, 1, 2):
        raise ValueError("dim must be 0, 1, or 2")
    return _proc_ppm(dic_3d[_PROC_KEYS_3D[dim]], data_3d.shape[dim])


def extract_2d_strip(dic_3d, data_3d, z_slice=None, slice_index=None, slice_axis=1):
    """
    Extracts a 2D plane from a 3D Bruker dataset by PPM or point index.
    Remaps the dictionary headers to match standard 2D processing functions.

    slice_axis:
      0: Fix F1 (proc3s) -> extracts F2-F3 plane
      1: Fix F2 (proc2s) -> extracts F1-F3 plane (most common, e.g. 15N slice in HNCACB)
      2: Fix F3 (procs)  -> extracts F1-F2 plane
    """
    ppm_axis = get_dim_ppm(dic_3d, data_3d, slice_axis)

    # If PPM is provided, find the closest point index automatically
    if slice_index is None:
        if z_slice is None:
            raise ValueError("Provide either z_slice or slice_index.")
        slice_index = int(np.argmin(np.abs(ppm_axis - z_slice)))

    actual_ppm = ppm_axis[slice_index]
    # Shallow copy is enough: only top-level entries are replaced, never modified
    dic_2d = dict(dic_3d)

    if slice_axis == 0:
        data_2d = data_3d[slice_index, :, :]
    elif slice_axis == 1:
        data_2d = data_3d[:, slice_index, :]
        dic_2d["proc2s"] = dic_3d["proc3s"]  # Map F1 (proc3s) to Y-axis (proc2s)
    else:
        data_2d = data_3d[:, :, slice_index]
        dic_2d["procs"] = dic_3d["proc2s"]
        dic_2d["proc2s"] = dic_3d["proc3s"]

    return dic_2d, data_2d, actual_ppm, slice_index


def _extract_slices(p, indices, z_slice, slice_index, slice_axis):
    """Extracts the same plane from several 3D spectra. Returns slices and the ppm value."""
    slices, actual_ppm_base = [], None
    for idx in indices:
        data_3d = p["data_all"][idx]
        dic_2d, data_2d, actual_ppm, s_idx = extract_2d_strip(
            p["dic_all"][idx], data_3d, z_slice, slice_index, slice_axis
        )
        if actual_ppm_base is None:
            actual_ppm_base = actual_ppm
        # Stable cache id: spectrum + plane (id() of the temporary slice would be reused!)
        cache_id = ("strip", idx, id(data_3d), slice_axis, s_idx)
        slices.append((idx, dic_2d, data_2d, cache_id))
    return slices, actual_ppm_base


def _strip_params(p, x_slice, strip_width):
    p_temp = p.copy()
    if x_slice is not None:
        half_width = strip_width / 2.0
        p_temp["xlim"] = (x_slice + half_width, x_slice - half_width)
    return p_temp


def _draw_strip(ax, p, slices, reverse=False, x_slice=None):
    handles, names = [], []
    for idx, dic_2d, data_2d, cache_id in reversed(slices) if reverse else slices:
        handles.append(
            draw_contours(
                ax, dic_2d, data_2d, p, p["cont"][idx], p["colors"][idx], cache_id=cache_id
            )
        )
        names.append(p["file_names"][idx])

    if x_slice is not None:
        ax.axvline(x=x_slice, color="black", linestyle="--", linewidth=0.5, alpha=0.5, zorder=0)
    return handles, names


def plot_strip(
    p,
    indices,
    z_slice=None,
    slice_index=None,
    slice_axis=1,
    name=None,
    reverse=False,
    x_slice=None,
    strip_width=0.3,
):
    """
    Plots a 2D strip from one or more 3D spectra.
    Accepts a single integer (e.g., indices=0) or a list (e.g., indices=[0, 1]).
    Can be optionally centered on an x_slice with a specific strip_width.
    """
    _check_style(p, _KEYS_2D, "plot_strip")
    indices = _valid(p, indices, 3)
    if not indices:
        print("  -> Error: No valid 3D spectra found.")
        return

    folder = p["out_strip"]
    result_folder(folder)

    slices, actual_ppm = _extract_slices(p, indices, z_slice, slice_index, slice_axis)
    file_names = [p["file_names"][i] for i in indices]

    # --- Naming ---
    is_overlay = len(slices) > 1
    style = p.get("spectra_type") or ""
    style_suffix = f"_{style}" if style else ""

    if name:
        base_name = name
        save_name = f"{base_name}{style_suffix}_{actual_ppm:.2f}ppm"
    elif is_overlay:
        base_name = f"Overlay_{'_'.join(file_names)}"
        save_name = f"{base_name}{style_suffix}_strip_overlay_{actual_ppm:.2f}ppm"
    else:
        base_name = file_names[0]
        save_name = f"{base_name}{style_suffix}_strip_{actual_ppm:.2f}ppm"

    title = f"{base_name} ({actual_ppm:.2f} ppm)"
    if x_slice is not None:
        title += f"\n{x_slice:.2f} ppm"

    p_temp = _strip_params(p, x_slice, strip_width)
    fig, ax = _new_figure(p_temp)
    apply_formatting(ax, p_temp, title=title, add_labels=True)

    handles, names = _draw_strip(ax, p_temp, slices, reverse, x_slice)
    if is_overlay:
        _add_legend(ax, handles, names, p_temp)

    # spectra_type is already part of save_name
    p_temp.pop("spectra_type", None)
    save_and_clear(fig, folder, save_name, p_temp)


def _parse_slice(slice_val):
    """A slice is either a ppm value or a dict with z_slice/slice_index/x_slice/title."""
    if isinstance(slice_val, dict):
        return (
            slice_val.get("z_slice"),
            slice_val.get("slice_index"),
            slice_val.get("x_slice"),
            slice_val.get("title"),
        )
    return float(slice_val), None, None, None


def strip_grid_plot(
    p,
    indices,
    slices,
    slice_axis=1,
    row=1,
    col=3,
    name="strip_grid",
    strip_width=0.3,
    reverse=False,
):
    """
    Plots a grid of 2D strips from one or more 3D spectra.
    Accepts a single integer (e.g., indices=0) or a list (e.g., indices=[0, 1]).
    """
    _check_style(p, _KEYS_2D, "strip_grid_plot", grid=True)
    indices = _valid(p, indices, 3)
    if not indices:
        print("  -> Error: No valid 3D spectra found.")
        return

    folder = p["out_strip"]
    result_folder(folder)

    is_overlay = len(indices) > 1
    if is_overlay and name == "strip_grid":
        name = "strip_grid_overlay"

    fig, axes = _new_grid(p, row, col, sharex=False, sharey=True)

    for i, (ax, slice_val) in enumerate(_zip_axes(axes, list(slices))):
        if slice_val is None:
            ax.axis("off")
            continue

        sl_z, sl_idx, sl_x, sl_title = _parse_slice(slice_val)
        panel, actual_ppm = _extract_slices(p, indices, sl_z, sl_idx, slice_axis)

        # Title: name + z-slice + x-slice
        title = f"{sl_title}\n({actual_ppm:.2f} ppm)" if sl_title else f"{actual_ppm:.2f} ppm"
        if sl_x is not None:
            title += f"\n{sl_x:.2f} ppm"

        p_temp = _strip_params(p, sl_x, strip_width)
        apply_formatting(ax, p_temp, title=title, is_grid=True)
        handles, names = _draw_strip(ax, p_temp, panel, reverse, sl_x)

        if not ax.get_subplotspec().is_last_row():
            ax.tick_params(labelbottom=False)

        # Legend only in the first panel for overlays
        if i == 0 and is_overlay:
            _add_legend(ax, handles, names, p_temp)

    save_and_clear(fig, folder, name, p)


# =====================================================================
# --- 1D FUNCTIONS ---
# =====================================================================


def _scale_factor(value):
    """For 1D spectra the 'contour' column of the data list is a scaling factor."""
    try:
        factor = float(value)
    except (TypeError, ValueError):
        return 1.0
    return factor if np.isfinite(factor) and factor != 0 else 1.0


def _norm_reference(p, idx, ppm, y, mode, region):
    if mode == "ns":
        acqus = p["dic_all"][idx].get("acqus", {})
        reference = float(acqus.get("NS", 1)) * float(acqus.get("RG", 1))
        return reference if reference > 0 else 1.0

    if mode == "max":
        limits = region
    elif mode == "region":
        limits = p.get("norm_region")
        if limits is None:
            raise ValueError("normalize='region' needs norm_region=(ppm_high, ppm_low).")
    else:
        raise ValueError(f"Unknown normalize='{mode}'. Use None, 'max', 'region' or 'ns'.")

    mask = _in_range(ppm, limits)
    if not mask.any():
        return 1.0
    reference = float(y[mask].max())
    if reference <= 0:
        reference = float(np.abs(y[mask]).max())
    return reference or 1.0


def get_1d(p, idx, region=None):
    """
    Returns (ppm, intensity) of a 1D spectrum after normalisation and scaling.

    p["normalize"]: None, "max" (highest point in 'region', default xlim),
                    "region" (highest point in p["norm_region"]) or "ns" (NS * RG).
    The 'contour' column of the data list is applied afterwards as scaling factor.
    """
    dic, data = p["dic_all"][idx], p["data_all"][idx]
    ppm = _proc_ppm(dic["procs"], data.shape[-1])
    y = np.real(data).astype(float)

    mode = p.get("normalize")
    if mode:
        y = y / _norm_reference(p, idx, ppm, y, mode, region if region is not None else p["xlim"])

    return ppm, y * _scale_factor(p["cont"][idx])


def apply_formatting_1d(ax, p, title=None, is_grid=False, add_labels=True):
    """Classic 1D NMR look: only the ppm axis at the bottom, no intensity axis."""
    # In grids the title sits top left inside the panel (high-ppm end is usually empty)
    _set_title(ax, p, title, is_grid, loc="left" if is_grid else "center")
    _format_x(ax, p, hide_first=is_grid and p.get("first_x", False))

    for side in ("left", "right", "top"):
        ax.spines[side].set_visible(False)
    ax.yaxis.set_visible(False)

    if add_labels and p["labels"]:
        ax.set_xlabel(p["x_label"])


def _trace(p, idx, offset=0.0, shift=0.0):
    """(idx, ppm, intensity, shift) of a 1D spectrum, optionally shifted for stacks."""
    ppm, y = get_1d(p, idx)
    return idx, ppm + shift, y + offset, shift


def _plot_traces(ax, p, traces, alpha=1.0):
    handles, names = [], []
    for idx, ppm, y, _ in traces:
        color = _split_colors(p["colors"][idx])[0]
        (line,) = ax.plot(*_crop(ppm, y, p["xlim"]), color=color, lw=p["line_width"], alpha=alpha)
        handles.append(line)
        names.append(p["file_names"][idx])
    return handles, names


def add_labels_1d(ax, csv_file, ppm, y, p, shift=0.0):
    """
    Peak labels for 1D spectra from a CSV (columns Residue, position_1).
    The label is placed above the highest point within label_window_1d ppm.
    """
    df = _read_label_csv(csv_file, p, ("Residue", "position_1"))
    if df is None:
        return [], [], []

    positions = df["position_1"].to_numpy(dtype=float) + shift
    labels = df["Residue"].astype(str).to_numpy()
    keep = _in_range(positions, p["xlim"])

    y0, y1 = ax.get_ylim()
    gap = p.get("label_offset_1d", 0.04) * abs(y1 - y0)
    half_window = p.get("label_window_1d", 0.02) / 2

    texts, xs, ys = [], [], []
    for label, pos in zip(labels[keep], positions[keep]):
        window = np.abs(ppm - pos) <= half_window
        peak = y[window].max() if window.any() else y[np.argmin(np.abs(ppm - pos))]
        texts.append(
            ax.text(
                pos,
                peak + gap,
                label,
                fontsize=p.get("label_fontsize", 6),
                rotation=p.get("label_rotation_1d", 90),
                color="black",
                ha="center",
                va="bottom",
                clip_on=False,
            )
        )
        xs.append(pos)
        ys.append(peak + 0.25 * gap)
    return texts, xs, ys


def adjust_labels_1d(ax, texts, xs, ys, p):
    """Moves 1D labels only sideways until they no longer overlap, then adds connectors."""
    if not texts:
        return []
    ax.figure.canvas.draw()
    with _quiet():
        adjust_text(
            texts,
            target_x=xs,
            target_y=ys,
            ax=ax,
            avoid_self=False,
            only_move={"text": "x", "static": "x", "explode": "x", "pull": "x"},
            expand=(1.15, 1.0),
            force_text=(0.5, 0.0),
            force_pull=(0.02, 0.0),
            prevent_crossings=False,
            ensure_inside_axes=False,
            iter_lim=20 * p.get("iterations", 10),
        )
    return draw_final_connectors(
        ax, xs, ys, texts, min_distance_px=2, lw=p.get("connector_width", 0.5), alpha=0.5
    )


def _ylim_1d(p, traces, has_labels):
    if p.get("ylim_1d") is not None:
        return p["ylim_1d"]
    headroom = p.get("label_headroom_1d", 0.3) if has_labels else 0.0
    return _auto_range([(t[1], t[2]) for t in traces], p["xlim"], top=0.05 + headroom)


def _traces_have_labels(p, traces):
    return bool(p.get("peak_labels")) and any(_has_csv(p["csv_files"][t[0]]) for t in traces)


def _finish_1d(ax, p, traces, ylim=None):
    """Sets the intensity range and adds the peak labels of all traces."""
    has_labels = _traces_have_labels(p, traces)
    ylim = ylim or _ylim_1d(p, traces, has_labels)
    if ylim is not None:
        ax.set_ylim(ylim)

    if has_labels:
        texts, xs, ys = [], [], []
        for idx, ppm, y, shift in traces:
            t, x, yy = add_labels_1d(ax, p["csv_files"][idx], ppm, y, p, shift)
            texts += t
            xs += x
            ys += yy
        adjust_labels_1d(ax, texts, xs, ys, p)


def plot_everything_1d(p):
    """Plots every loaded 1D spectrum individually (folder '1D')."""
    _check_style(p, _KEYS_1D, "plot_everything_1d", ylabel=False)
    folder = p["out_1d"]
    result_folder(folder)

    for idx in _valid(p, range(len(p["files"])), 1, verbose=False):
        fig, ax = _new_figure(p)
        apply_formatting_1d(ax, p, title=p["file_names"][idx])
        traces = [_trace(p, idx)]
        _plot_traces(ax, p, traces)
        _finish_1d(ax, p, traces)
        save_and_clear(fig, folder, f"{p['file_names'][idx]}", p)


def overlay_1d(p, get, name):
    """Overlays the 1D spectra 'get' (list of indices); the last one is on top."""
    _check_style(p, _KEYS_1D, "overlay_1d", ylabel=False)
    indices = _valid(p, get, 1)
    if not indices:
        print("  -> Error: No valid 1D spectra found.")
        return

    folder = p["out_1d"]
    result_folder(folder)

    fig, ax = _new_figure(p)
    apply_formatting_1d(ax, p, title=name)
    traces = [_trace(p, idx) for idx in indices]
    handles, names = _plot_traces(ax, p, traces, alpha=p["alpha"])
    _finish_1d(ax, p, traces)
    _add_legend(ax, handles, names, p)

    save_and_clear(fig, folder, f"overlay_{name}", p)


def stack_1d(p, get, name, offset=None, shift=0.0, names="text"):
    """
    Stacked 1D plot: the first index is at the bottom, the last one on top.

    offset: vertical distance in units of the highest visible peak
            (default: p["stack_offset"])
    shift:  horizontal shift per spectrum in ppm (pseudo-3D look)
    names:  "text" (name above each baseline), "legend" or None
    """
    _check_style(p, _KEYS_1D, "stack_1d", ylabel=False)
    indices = _valid(p, get, 1)
    if not indices:
        print("  -> Error: No valid 1D spectra found.")
        return

    folder = p["out_1d"]
    result_folder(folder)

    # Step size from the highest visible peak of all spectra
    height = 0.0
    for idx in indices:
        ppm, y = get_1d(p, idx)
        mask = _in_range(ppm, p["xlim"])
        if mask.any():
            height = max(height, float(y[mask].max()))
    step = (offset if offset is not None else p.get("stack_offset", 1.0)) * (height or 1.0)

    fig, ax = _new_figure(p)
    apply_formatting_1d(ax, p, title=name)
    traces = [_trace(p, idx, offset=k * step, shift=k * shift) for k, idx in enumerate(indices)]
    handles, labels = _plot_traces(ax, p, traces)

    if names == "text":
        for k, idx in enumerate(indices):
            ax.text(
                0.01,
                k * step + 0.05 * (height or 1.0),
                p["file_names"][idx],
                transform=ax.get_yaxis_transform(),  # x in axes fraction, y in data
                color=_split_colors(p["colors"][idx])[0],
                fontsize=plt.rcParams["legend.fontsize"],
                ha="left",
                va="bottom",
            )
    elif names == "legend":
        # Reverse so that the legend order matches the vertical order
        _add_legend(ax, handles[::-1], labels[::-1], p)

    _finish_1d(ax, p, traces)
    save_and_clear(fig, folder, f"stack_{name}", p)


def grid_plot_1d(p, indices=None, row=2, col=1, name="grid_1D", sharey=False):
    """
    Plots 1D spectra in a row x col grid with a common ppm axis.
    indices=None uses all 1D spectra. sharey=True uses one intensity scale for all panels.
    """
    _check_style(p, _KEYS_1D, "grid_plot_1d", ylabel=False)
    pool = range(len(p["files"])) if indices is None else indices
    indices = _valid(p, pool, 1, verbose=indices is not None)
    if not indices:
        print("  -> Error: No valid 1D spectra found.")
        return

    folder = p["out_1d"]
    result_folder(folder)

    fig, axes = _new_grid(p, row, col, sharex=True, sharey=sharey, suplabels=False)
    panels = list(_zip_axes(axes, indices))
    traces = {idx: _trace(p, idx) for _, idx in panels if idx is not None}

    # Shared intensity scale must be computed from all panels at once
    common_ylim = None
    if sharey:
        common_ylim = _ylim_1d(p, list(traces.values()), _traces_have_labels(p, traces.values()))

    for ax, idx in panels:
        if idx is None:
            ax.axis("off")
            continue
        # x-axis label below the bottom row (auto-positioned, unlike supxlabel)
        bottom = ax.get_subplotspec().is_last_row()
        apply_formatting_1d(ax, p, title=p["file_names"][idx], is_grid=True, add_labels=bottom)
        _plot_traces(ax, p, [traces[idx]])
        _finish_1d(ax, p, [traces[idx]], ylim=common_ylim)

    save_and_clear(fig, folder, name, p)


# =====================================================================
# --- UTILITIES ---
# =====================================================================


def set_style(p_dict, style_name, styles_dict=None):
    """
    Applies a spectra style from SPECTRA_STYLES to p (in-place).

    Parameters set by the previous style are first reset to their original value,
    so nothing leaks from one style into the next. Parameters that are not part of
    any style (and changes you made to them in p) are left untouched.
    """
    target_styles = styles_dict or p_dict.get("SPECTRA_STYLES")
    if not target_styles:
        raise ValueError("No SPECTRA_STYLES found in p or passed to set_style.")
    if style_name not in target_styles:
        raise KeyError(
            f"Style '{style_name}' not found. Available: {', '.join(map(str, target_styles))}"
        )

    defaults = p_dict.setdefault("_style_defaults", {})

    # 1. Undo the previous style
    for key in p_dict.get("_style_keys", ()):
        value = defaults.get(key, _MISSING)
        if value is _MISSING:
            p_dict.pop(key, None)
        else:
            p_dict[key] = value

    # 2. Apply the new one (remember the original value of every key it touches)
    style = target_styles[style_name]
    for key, value in style.items():
        if key not in defaults:
            defaults[key] = p_dict.get(key, _MISSING)
        p_dict[key] = value

    p_dict["_style_keys"] = list(style)
    p_dict["spectra_type"] = style_name  # Keeps track of the name for file exports


def spectrum_menu(p):
    """Prints a terminal menu of loaded spectra and their dimensionality."""
    print("\n----- Spectrum Menu -----")
    for i, name in enumerate(p["file_names"]):
        data = p["data_all"][i]
        dim_text = f"{data.ndim}D" if data is not None else "Missing/Error"
        print(f" [{i}] -> {name} ({dim_text})")
    print("-------------------------\n")


def colormap(files, maps="viridis"):
    cmap = plt.get_cmap(maps)
    colors_array = cmap(np.linspace(0, 1, len(files)))
    # Convert the array to a list of tuples
    return [tuple(c) for c in colors_array]
